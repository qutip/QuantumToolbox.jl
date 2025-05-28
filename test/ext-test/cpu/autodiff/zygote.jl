@testitem "Zygote Extension" tags=[:autodiff] default_imports=false begin
    using Test
    using QuantumToolbox
    using Zygote
    using Enzyme
    using SciMLSensitivity

    @testset "sesolve" begin
        coef_Ω(p, t) = p[1]

        H = QobjEvo(sigmax(), coef_Ω)
        ψ0 = fock(2, 1)
        t_max = 10

        function my_f_sesolve(p)
            tlist = range(0, t_max, 100)

            sol = sesolve(
                H,
                ψ0,
                tlist,
                progress_bar = Val(false),
                params = p,
                sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
            )

            return real(expect(projection(2, 0, 0), sol.states[end]))
        end

        # Analytical solution
        my_f_analytic(Ω) = abs2(sin(Ω * t_max))
        my_f_analytic_deriv(Ω) = 2 * t_max * sin(Ω * t_max) * cos(Ω * t_max)

        Ω = 1.0
        params = [Ω]

        my_f_analytic(Ω)
        my_f_sesolve(params)

        grad_qt = Zygote.gradient(my_f_sesolve, params)[1]
        grad_exact = [my_f_analytic_deriv(params[1])]

        @test grad_qt ≈ grad_exact atol=1e-6
    end

    @testset "mesolve" begin
        N = 20
        a = destroy(N)

        coef_Δ(p, t) = p[1]
        coef_F(p, t) = p[2]
        coef_γ(p, t) = sqrt(p[3])

        H = QobjEvo(a' * a, coef_Δ) + QobjEvo(a + a', coef_F)
        c_ops = [QobjEvo(a, coef_γ)]
        L = liouvillian(H, c_ops)

        ψ0 = fock(N, 0)

        function my_f_mesolve(p)
            tlist = range(0, 40, 100)

            sol = mesolve(
                L,
                ψ0,
                tlist,
                progress_bar = Val(false),
                params = p,
                sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
            )

            return real(expect(a' * a, sol.states[end]))
        end

        # Analytical solution
        n_ss(Δ, F, γ) = abs2(F / (Δ + 1im * γ / 2))

        Δ = 1.0
        F = 1.0
        γ = 1.0
        params = [Δ, F, γ]

        # The factor 2 is due to a bug
        grad_qt = Zygote.gradient(my_f_mesolve, params)[1]

        grad_exact = Zygote.gradient((p) -> n_ss(p[1], p[2], p[3]), params)[1]

        @test grad_qt ≈ grad_exact atol=1e-6
    end
end
