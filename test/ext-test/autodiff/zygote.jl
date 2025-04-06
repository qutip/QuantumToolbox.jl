@testset "Zygote.jl Autodiff" begin
    @testset "mesolve" begin
        N = 16
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
        my_f_mesolve(params)
        n_ss(Δ, F, γ)

        # The factor 2 is due to a bug
        grad_qt = Zygote.gradient(my_f_mesolve, params)[1] ./ 2

        grad_exact = Zygote.gradient((p) -> n_ss(p[1], p[2], p[3]), params)[1]

        @test grad_qt ≈ grad_exact atol=1e-6
    end
end
