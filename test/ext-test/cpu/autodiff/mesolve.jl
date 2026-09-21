using QuantumToolbox
using ForwardDiff
using Enzyme
using Mooncake
using SciMLSensitivity
using SciMLSensitivity: MooncakeVJP

# setup problem
const N = 20
const a = destroy(N)
const ψ0_mesolve = fock(N, 0)
const tlist_mesolve = range(0, 40, 100)
const ad_a = a' * a

# For direct Forward differentiation
function my_f_mesolve_direct(p)
    H = p[1] * a' * a + p[2] * (a + a')
    c_ops = [sqrt(p[3]) * a]
    sol = mesolve(H, ψ0_mesolve, tlist_mesolve, c_ops, progress_bar = Val(false))
    return real(expect(ad_a, sol.states[end]))
end

# For SciMLSensitivity.jl
coef_Δ(p, t) = p[1]
coef_F(p, t) = p[2]
coef_γ(p, t) = sqrt(p[3])
H = QobjEvo(a' * a, coef_Δ) + QobjEvo(a + a', coef_F)
c_ops = [QobjEvo(a, coef_γ)]
const L = liouvillian(H, c_ops)
const L_assume_non_herm = liouvillian(H, c_ops, assume_hermitian = Val(false))

function my_f_mesolve(p, sensealg)
    sol = mesolve(
        L,
        ψ0_mesolve,
        tlist_mesolve,
        progress_bar = Val(false),
        params = p,
        sensealg = sensealg,
    )

    return real(expect(ad_a, sol.states[end]))
end
const my_f_mesolve_bsa_enzyme = Base.Fix{2}(my_f_mesolve, BacksolveAdjoint(autojacvec = EnzymeVJP()))
const my_f_mesolve_bsa_mooncake = Base.Fix{2}(my_f_mesolve, BacksolveAdjoint(autojacvec = MooncakeVJP()))

function my_f_mesolve_assume_non_herm(p, sensealg)
    sol = mesolve(
        L_assume_non_herm,
        ψ0_mesolve,
        tlist_mesolve,
        progress_bar = Val(false),
        params = p,
        sensealg = sensealg,
    )

    return real(expect(ad_a, sol.states[end]))
end
const my_f_mesolve_assume_non_herm_bsa_enzyme = Base.Fix{2}(my_f_mesolve_assume_non_herm, BacksolveAdjoint(autojacvec = EnzymeVJP()))
const my_f_mesolve_assume_non_herm_bsa_mooncake = Base.Fix{2}(my_f_mesolve_assume_non_herm, BacksolveAdjoint(autojacvec = MooncakeVJP()))

# Analytical solution
n_ss(Δ, F, γ) = abs2(F / (Δ + 1im * γ / 2))

@testset "mesolve" verbose = true begin
    Δ = 1.0
    F = 1.0
    γ = 1.0
    params = [Δ, F, γ]

    my_f_mesolve_direct(params)
    my_f_mesolve_bsa_enzyme(params)
    my_f_mesolve_bsa_mooncake(params)
    my_f_mesolve_assume_non_herm_bsa_enzyme(params)
    my_f_mesolve_assume_non_herm_bsa_mooncake(params)

    grad_exact_cache = Mooncake.prepare_gradient_cache(splat(n_ss), params)
    grad_exact = Mooncake.value_and_gradient!!(grad_exact_cache, splat(n_ss), params)[2][2]

    @testset "ForwardDiff.jl" begin
        grad_qt = ForwardDiff.gradient(my_f_mesolve_direct, params)
        @test grad_qt ≈ grad_exact atol = 1.0e-6
    end

    @testset "Mooncake.jl" begin
        grad_cache1 = Mooncake.prepare_gradient_cache(my_f_mesolve_bsa_mooncake, params)
        # grad_cache2 = Mooncake.prepare_gradient_cache(my_f_mesolve_assume_non_herm_bsa_mooncake, params)
        _, grad_mooncake1 = Mooncake.value_and_gradient!!(grad_cache1, my_f_mesolve_bsa_mooncake, params)
        @test grad_mooncake1[2] ≈ grad_exact atol = 1.0e-6
        # _, grad_mooncake2 = Mooncake.value_and_gradient!!(grad_cache2, my_f_mesolve_assume_non_herm_bsa_mooncake, params)
        # @test grad_mooncake2[2] ≈ grad_exact atol = 1.0e-6
    end

    @testset "Enzyme.jl" begin
        dparams1 = Enzyme.make_zero(params)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.Reverse),
            my_f_mesolve_bsa_enzyme,
            Active,
            Duplicated(params, dparams1),
        )[1]

        dparams2 = Enzyme.make_zero(params)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.Reverse),
            my_f_mesolve_assume_non_herm_bsa_enzyme,
            Active,
            Duplicated(params, dparams2),
        )[1]

        @test dparams1 ≈ grad_exact atol = 1.0e-6
        @test dparams2 ≈ grad_exact atol = 1.0e-6
    end
end
