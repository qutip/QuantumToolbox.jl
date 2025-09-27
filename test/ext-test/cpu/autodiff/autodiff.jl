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

function my_f_sesolve(p)
    sol = sesolve(
        H_evo,
        ψ0_sesolve,
        tlist_sesolve,
        progress_bar = Val(false),
        params = p,
        sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
    )

    return real(expect(projection(2, 0, 0), sol.states[end]))
end

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

function my_f_mesolve(p)
    sol = mesolve(
        L,
        ψ0_mesolve,
        tlist_mesolve,
        progress_bar = Val(false),
        params = p,
        sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
    )

    return real(expect(a' * a, sol.states[end]))
end

# Analytical solution
n_ss(Δ, F, γ) = abs2(F / (Δ + 1im * γ / 2))

@testset "Autodiff" verbose=true begin
    @testset "sesolve" verbose=true begin
        Ω = 1.0
        params = [Ω]

        my_f_sesolve_direct(params)
        my_f_sesolve(params)

        grad_exact = [my_f_analytic_deriv(params[1])]

        @testset "ForwardDiff.jl" begin
            grad_qt = ForwardDiff.gradient(my_f_sesolve_direct, params)

            @test grad_qt ≈ grad_exact atol=1e-6
        end

        @testset "Zygote.jl" begin
            grad_qt = Zygote.gradient(my_f_sesolve, params)[1]

            @test grad_qt ≈ grad_exact atol=1e-6
        end

        @testset "Enzyme.jl" begin
            dparams = Enzyme.make_zero(params)
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse),
                my_f_sesolve,
                Active,
                Duplicated(params, dparams),
            )[1]

            @test dparams ≈ grad_exact atol=1e-6
        end
    end

    @testset "mesolve" verbose=true begin
        Δ = 1.0
        F = 1.0
        γ = 1.0
        params = [Δ, F, γ]

        my_f_mesolve_direct(params)
        my_f_mesolve(params)

        grad_exact = Zygote.gradient((p) -> n_ss(p[1], p[2], p[3]), params)[1]

        @testset "ForwardDiff.jl" begin
            grad_qt = ForwardDiff.gradient(my_f_mesolve_direct, params)
            @test grad_qt ≈ grad_exact atol=1e-6
        end

        @testset "Zygote.jl" begin
            grad_qt = Zygote.gradient(my_f_mesolve, params)[1]
            @test grad_qt ≈ grad_exact atol=1e-6
        end

        @testset "Enzyme.jl" begin
            dparams = Enzyme.make_zero(params)
            Enzyme.autodiff(
                Enzyme.set_runtime_activity(Enzyme.Reverse),
                my_f_mesolve,
                Active,
                Duplicated(params, dparams),
            )[1]

            @test dparams ≈ grad_exact atol=1e-6
        end
    end
end
