using QuantumToolbox
using ForwardDiff
using Enzyme
using Mooncake
using SciMLSensitivity
using SciMLSensitivity: MooncakeVJP

# setup problem
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
