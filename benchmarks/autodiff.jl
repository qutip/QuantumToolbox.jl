# Use harmonic oscillator system for both sesolve and mesolve
const N_ad = 20
const a_ad = destroy(N_ad)
const ψ0_ad = fock(N_ad, 0)
const tlist_ad = range(0, 40, 100)
const ad_a = a_ad' * a_ad

# NOTE ON WHAT IS AND IS NOT COMPARABLE HERE
#
# Forward and reverse mode necessarily measure *different* formulations of the same
# physics, and this is a library constraint rather than an oversight:
#
#   * ForwardDiff cannot differentiate through a `QobjEvo`. `QobjEvo` builds each
#     coefficient as `ScalarOperator(zero(eltype(op)), update_func)`, fixing the
#     coefficient element type to the operator's at construction time, so a `Dual`
#     parameter has nowhere type-correct to land and `update_coefficients!` throws
#     `MethodError: no method matching Float64(::ForwardDiff.Dual)`.
#   * Reverse mode cannot differentiate the constant/closure formulation. Continuous
#     adjoints form `λᵀ ∂f/∂p`, so parameters must live in the ODE's `p`. Building `H`
#     inside the function leaves `params = NullParameters()`, which either errors or —
#     with `BacksolveAdjoint(autojacvec = EnzymeVJP())` — silently returns an all-zero
#     gradient while running *faster* than the correct version.
#
# So the `Forward` entries below solve a constant `QuantumObject` problem (one fused
# sparse matvec per RHS evaluation, no forced `tstops`) while the `Reverse` entries
# solve a `QobjEvo` problem (three lazy matvecs, `tstops = tlist`). The `Primal`
# entries measure exactly the ODE the `Reverse` entries differentiate, so the tracked
# charts show AD *overhead* rather than an absolute time that moves with every
# dependency bump.

# ---- SESOLVE ----
# For direct Forward differentiation
function my_f_sesolve_direct(p)
    H = p[1] * a_ad' * a_ad + p[2] * (a_ad + a_ad')
    sol = sesolve(H, ψ0_ad, tlist_ad, progress_bar = Val(false))
    return real(expect(ad_a, sol.states[end]))
end

# For SciMLSensitivity.jl (reverse mode with Mooncake and Enzyme)
coef_Δ(p, t) = p[1]
coef_F(p, t) = p[2]
const H_evo = QobjEvo(a_ad' * a_ad, coef_Δ) + QobjEvo(a_ad + a_ad', coef_F)

function my_f_sesolve(p, sensealg)
    sol = sesolve(
        H_evo,
        ψ0_ad,
        tlist_ad,
        progress_bar = Val(false),
        params = p,
        sensealg = sensealg,
    )
    return real(expect(ad_a, sol.states[end]))
end

# `sesolve` is unitary, so the reverse solve of `BacksolveAdjoint` is stable and it is
# the only adjoint that reproduces the analytic gradient here:
#   analytic/ForwardDiff/central-difference all give [52.941064, 6.667844]
#   BacksolveAdjoint                                 [52.941063, 6.667844]  ✓
#   InterpolatingAdjoint(checkpointing = true)       [48.192356, 8.345162]  ✗ 9% off
#   GaussAdjoint                                    [-204.57,   145.85   ]  ✗
const my_f_sesolve_bsa_enzyme = Base.Fix{2}(my_f_sesolve, BacksolveAdjoint(autojacvec = EnzymeVJP()))
const my_f_sesolve_bsa_mooncake = Base.Fix{2}(my_f_sesolve, BacksolveAdjoint(autojacvec = MooncakeVJP()))

# ---- MESOLVE ----
# For direct Forward differentiation
function my_f_mesolve_direct(p)
    H = p[1] * a_ad' * a_ad + p[2] * (a_ad + a_ad')
    c_ops = [sqrt(p[3]) * a_ad]
    sol = mesolve(H, ψ0_ad, tlist_ad, c_ops, progress_bar = Val(false))
    return real(expect(ad_a, sol.states[end]))
end

# For SciMLSensitivity.jl (reverse mode with Mooncake and Enzyme)
coef_γ(p, t) = sqrt(p[3])
const c_ops_ad = [QobjEvo(a_ad, coef_γ)]
const L_ad = liouvillian(H_evo, c_ops_ad)

function my_f_mesolve(p, sensealg)
    sol = mesolve(
        L_ad,
        ψ0_ad,
        tlist_ad,
        progress_bar = Val(false),
        params = p,
        sensealg = sensealg,
    )
    return real(expect(ad_a, sol.states[end]))
end

# Lindblad dynamics is contracting, so integrating it backwards is *expanding* and
# `BacksolveAdjoint`'s reverse solve is unstable: it only stays correct here because
# `saveat = tlist` incidentally supplies 100 checkpoints. With `saveat = [tlist[end]]`
# it aborts with `dt` below floating-point epsilon and returns `NaN`.
# `InterpolatingAdjoint(checkpointing = true)` never re-integrates the state backwards,
# and on this problem it is correct and ~1.6x faster on both engines
# (Mooncake 805 -> 510 ms; Enzyme 702 -> 435 ms, 883 -> 509 MB).
const mesolve_sensealg_enzyme = InterpolatingAdjoint(autojacvec = EnzymeVJP(), checkpointing = true)
const mesolve_sensealg_mooncake = InterpolatingAdjoint(autojacvec = MooncakeVJP(), checkpointing = true)
const my_f_mesolve_bsa_enzyme = Base.Fix{2}(my_f_mesolve, mesolve_sensealg_enzyme)
const my_f_mesolve_bsa_mooncake = Base.Fix{2}(my_f_mesolve, mesolve_sensealg_mooncake)

# Parameters for benchmarks
const params_sesolve = [1.0, 1.0]
const params_mesolve = [1.0, 1.0, 1.0]

function benchmark_autodiff!(SUITE)
    # Primal references: the exact ODEs the Reverse entries differentiate, with no AD.
    # Dividing a Reverse timing by its Primal timing gives the AD overhead factor, which
    # is what we actually want to track for regressions.
    SUITE["Autodiff"]["sesolve"]["Primal"] =
        @benchmarkable my_f_sesolve($params_sesolve, nothing)
    SUITE["Autodiff"]["mesolve"]["Primal"] =
        @benchmarkable my_f_mesolve($params_mesolve, nothing)

    # Benchmark sesolve - Forward
    SUITE["Autodiff"]["sesolve"]["Forward"] = @benchmarkable ForwardDiff.gradient($my_f_sesolve_direct, $params_sesolve)

    # Benchmark sesolve - Reverse (Enzyme)
    SUITE["Autodiff"]["sesolve"]["Reverse (Enzyme)"] = @benchmarkable Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Const($my_f_sesolve_bsa_enzyme),
        Active,
        Duplicated($params_sesolve, dparams_sesolve),
    ) setup = (dparams_sesolve = Enzyme.make_zero($params_sesolve))

    # Benchmark sesolve - Reverse (Mooncake)
    cache_sesolve_mooncake = Mooncake.prepare_gradient_cache(my_f_sesolve_bsa_mooncake, params_sesolve)
    SUITE["Autodiff"]["sesolve"]["Reverse (Mooncake)"] = @benchmarkable Mooncake.value_and_gradient!!(grad_cache1, $my_f_sesolve_bsa_mooncake, $params_sesolve) setup = (grad_cache1 = $cache_sesolve_mooncake)

    # Benchmark mesolve - Forward
    SUITE["Autodiff"]["mesolve"]["Forward"] = @benchmarkable ForwardDiff.gradient($my_f_mesolve_direct, $params_mesolve)

    # Benchmark mesolve - Reverse (Enzyme)
    SUITE["Autodiff"]["mesolve"]["Reverse (Enzyme)"] = @benchmarkable Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Const($my_f_mesolve_bsa_enzyme),
        Active,
        Duplicated($params_mesolve, dparams_mesolve),
    ) setup = (dparams_mesolve = Enzyme.make_zero($params_mesolve))

    # Benchmark mesolve - Reverse (Mooncake)
    cache_mesolve_mooncake = Mooncake.prepare_gradient_cache(my_f_mesolve_bsa_mooncake, params_mesolve)
    SUITE["Autodiff"]["mesolve"]["Reverse (Mooncake)"] = @benchmarkable Mooncake.value_and_gradient!!(grad_cache2, $my_f_mesolve_bsa_mooncake, $params_mesolve) setup = (grad_cache2 = $cache_mesolve_mooncake)

    return nothing
end
