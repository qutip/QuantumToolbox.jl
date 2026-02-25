# ---- SESOLVE ----
const ψ0_sesolve = fock(2, 1)
t_max = 10
const tlist_sesolve = range(0, t_max, 100)

# For direct Forward differentiation
function my_f_sesolve_direct(p)
    H = p[1] * sigmax()
    sol = sesolve(H, ψ0_sesolve, tlist_sesolve, progress_bar = Val(false))

    return real(expect(projection(2, 0, 0), sol.states[end]))
end

# For SciMLSensitivity.jl
coef_Ω(p, t) = p[1]
const H_evo = QobjEvo(sigmax(), coef_Ω)

function my_f_sesolve(p, sensealg)
    sol = sesolve(
        H_evo,
        ψ0_sesolve,
        tlist_sesolve,
        progress_bar = Val(false),
        params = p,
        sensealg = sensealg,
    )

    return real(expect(projection(2, 0, 0), sol.states[end]))
end
const my_f_sesolve_bsa_enzyme = Base.Fix{2}(my_f_sesolve, BacksolveAdjoint(autojacvec = EnzymeVJP()))
const my_f_sesolve_bsa_mooncake = Base.Fix{2}(my_f_sesolve, BacksolveAdjoint(autojacvec = MooncakeVJP()))

# Analytical solution
my_f_analytic(Ω) = abs2(sin(Ω * t_max))
my_f_analytic_deriv(Ω) = 2 * t_max * sin(Ω * t_max) * cos(Ω * t_max)

# ---- MESOLVE ----
const N = 20
const a = destroy(N)
const ψ0_mesolve = fock(N, 0)
const tlist_mesolve = range(0, 40, 100)

# For direct Forward differentiation
function my_f_mesolve_direct(p)
    H = p[1] * a' * a + p[2] * (a + a')
    c_ops = [sqrt(p[3]) * a]
    sol = mesolve(H, ψ0_mesolve, tlist_mesolve, c_ops, progress_bar = Val(false))
    return real(expect(a' * a, sol.states[end]))
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

    return real(expect(a' * a, sol.states[end]))
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

    return real(expect(a' * a, sol.states[end]))
end
const my_f_mesolve_assume_non_herm_bsa_enzyme = Base.Fix{2}(my_f_mesolve_assume_non_herm, BacksolveAdjoint(autojacvec = EnzymeVJP()))
const my_f_mesolve_assume_non_herm_bsa_mooncake = Base.Fix{2}(my_f_mesolve_assume_non_herm, BacksolveAdjoint(autojacvec = MooncakeVJP()))


# Analytical solution
n_ss(Δ, F, γ) = abs2(F / (Δ + 1im * γ / 2))

@testset "Autodiff" verbose = true begin
    @testset "ForwardDiff for thermal_dm" begin
        N = 100  # dimension of the system
        N_op = num(N) # number operator

        function N_expect(p) # p = [n]
            ρT = thermal_dm(N, p[1]; sparse = Val(true))
            return real(expect(N_op, ρT))
        end

        # Average photon number for thermal state
        n = rand(Float64)
        p = [n]

        # Use ForwardDiff.gradient to compute the gradient
        grad = ForwardDiff.gradient(N_expect, p)[1]

        # Compare the result
        @test isapprox(grad, 1.0; atol = 1.0e-6)
    end

    @testset "sesolve" verbose = true begin
        Ω = 1.0
        params = [Ω]

        my_f_sesolve_direct(params)
        my_f_sesolve_bsa_enzyme(params)
        my_f_sesolve_bsa_mooncake(params)

        grad_exact = [my_f_analytic_deriv(params[1])]

        @testset "ForwardDiff.jl" begin
            grad_qt = ForwardDiff.gradient(my_f_sesolve_direct, params)

            @test grad_qt ≈ grad_exact atol = 1.0e-6
        end

        @testset "Zygote.jl" begin
            grad_qt_bsa_enzyme = Zygote.gradient(my_f_sesolve_bsa_enzyme, params)[1]
            grad_qt_bsa_mooncake = Zygote.gradient(my_f_sesolve_bsa_mooncake, params)[1]

            # TODO: Fix the factor of 2 discrepancy once https://github.com/SciML/SciMLSensitivity.jl/issues/1181 is resolved.
            @test grad_qt_bsa_enzyme ≈ 2 .* grad_exact atol = 1.0e-6
            @test grad_qt_bsa_mooncake ≈ 2 .* grad_exact atol = 1.0e-6
        end

        @testset "Mooncake.jl" begin
            grad_cache = Mooncake.prepare_gradient_cache(my_f_sesolve_bsa_mooncake, params)
            _, grad_mooncake = Mooncake.value_and_gradient!!(grad_cache, my_f_sesolve_bsa_mooncake, params)
            @test grad_mooncake[2] ≈ grad_exact atol = 1.0e-6
        end

        @testset "Enzyme.jl" begin
            dparams = Enzyme.make_zero(params)
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse),
                my_f_sesolve_bsa_enzyme,
                Active,
                Duplicated(params, dparams),
            )[1]

            @test dparams ≈ grad_exact atol = 1.0e-6
        end
    end

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

        grad_exact = Zygote.gradient((p) -> n_ss(p[1], p[2], p[3]), params)[1]

        @testset "ForwardDiff.jl" begin
            grad_qt = ForwardDiff.gradient(my_f_mesolve_direct, params)
            @test grad_qt ≈ grad_exact atol = 1.0e-6
        end

        @testset "Zygote.jl" begin
            grad_qt1_bsa_enzyme = Zygote.gradient(my_f_mesolve_bsa_enzyme, params)[1]
            grad_qt1_bsa_mooncake = Zygote.gradient(my_f_mesolve_bsa_mooncake, params)[1]
            # grad_qt2_bsa_enzyme = Zygote.gradient(my_f_mesolve_assume_non_herm_bsa_enzyme, params)[1] # It doesn't work yet
            grad_qt2_bsa_mooncake = Zygote.gradient(my_f_mesolve_assume_non_herm_bsa_mooncake, params)[1]

            # TODO: Fix the factor of 2 discrepancy once

            @test grad_qt1_bsa_enzyme ≈ 2 .* grad_exact atol = 1.0e-6
            @test grad_qt1_bsa_mooncake ≈ 2 .* grad_exact atol = 1.0e-6
            # @test grad_qt2_bsa_enzyme ≈ 2 .* grad_exact atol = 1.0e-6
            @test grad_qt2_bsa_mooncake ≈ 2 .* grad_exact atol = 1.0e-6
        end

        @testset "Mooncake.jl" begin
            grad_cache1 = Mooncake.prepare_gradient_cache(my_f_mesolve_bsa_mooncake, params)
            grad_cache2 = Mooncake.prepare_gradient_cache(my_f_mesolve_assume_non_herm_bsa_mooncake, params)
            _, grad_mooncake1 = Mooncake.value_and_gradient!!(grad_cache1, my_f_mesolve_bsa_mooncake, params)
            @test grad_mooncake1[2] ≈ grad_exact atol = 1.0e-6
            _, grad_mooncake2 = Mooncake.value_and_gradient!!(grad_cache2, my_f_mesolve_assume_non_herm_bsa_mooncake, params)
            @test grad_mooncake2[2] ≈ grad_exact atol = 1.0e-6
        end

        @testset "Enzyme.jl" begin
            dparams1 = Enzyme.make_zero(params)
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse),
                my_f_mesolve_bsa_enzyme,
                Active,
                Duplicated(params, dparams1),
            )[1]

            # It doesn't work yet when assume_hermitian = Val(false)
            # dparams2 = Enzyme.make_zero(params)
            # Enzyme.autodiff(
            #     Enzyme.set_runtime_activity(Enzyme.Reverse),
            #     my_f_mesolve_assume_non_herm,
            #     Active,
            #     Duplicated(params, dparams2),
            # )[1]

            @test dparams1 ≈ grad_exact atol = 1.0e-6
            # @test dparams2 ≈ grad_exact atol = 1.0e-6
        end
    end
end
