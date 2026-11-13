using Test
using QuantumToolbox
using LinearAlgebra
using SciMLBase: SciMLBase

# Reference vector-Jacobian products of the (linear) right-hand side f(du, u, p, t):
# - state: J(p, t)ᴴ λ, with J built column by column from f applied to the basis states
# - parameters: central finite differences of p ↦ Re⟨λ, f(u, p, t)⟩
function reference_vjp(f, λ, u, p, t)
    du = similar(u)
    e = zero(u)
    dλ = similar(u)
    for k in eachindex(u)
        e[k] = 1
        f(du, e, p, t)
        e[k] = 0
        dλ[k] = dot(du, λ) # (Jᴴ λ)ₖ = Σᵢ conj(Jᵢₖ) λᵢ
    end
    return dλ
end

function reference_param_vjp(f, λ, u, p, t; h = 1.0e-6)
    du = similar(u)
    g(x) = (f(du, u, x, t); real(dot(λ, du)))
    return map(eachindex(p)) do j
        e = zero(p)
        e[j] = h
        (g(p + e) - g(p - e)) / (2h)
    end
end

function test_analytic_vjp(prob; t = 0.37)
    f = prob.f
    u = randn(ComplexF64, size(prob.u0))
    λ = randn(ComplexF64, size(prob.u0))
    p = prob.p
    dλ = similar(u)
    f.vjp(dλ, λ, u, p, t)
    @test dλ ≈ reference_vjp(f, λ, u, p, t)
    if f.vjp_p !== nothing && p isa AbstractVector
        dp = zero(p)
        f.vjp_p(dp, λ, u, p, t)
        @test dp ≈ reference_param_vjp(f, λ, u, p, t) rtol = 1.0e-6 atol = 1.0e-8
    end
    return nothing
end

@testset "Analytic VJPs of the time evolution ODEs" begin
    N = 6
    a = destroy(N)
    ψ0 = basis(N, 0)
    tlist = range(0, 1, 5)
    params = [0.8, 1.3, 0.4]

    coef_drive(p, t) = p[1] * exp(1im * p[2] * t) # complex coefficient
    coef_decay(p, t) = sqrt(p[3]) * cos(t)
    H = QobjEvo(a' * a) + QobjEvo(a, coef_drive) + QobjEvo(a', (p, t) -> conj(coef_drive(p, t)))
    c_ops = [QobjEvo(a, coef_decay), sqrt(0.2) * a']

    @testset "sesolve" begin
        prob = sesolveProblem(H, ψ0, tlist; params = params, progress_bar = Val(false)).prob
        @test prob.f.vjp isa QuantumToolbox.LinearOperatorVJP
        @test prob.f.vjp_p isa QuantumToolbox.LinearOperatorParamVJP
        test_analytic_vjp(prob)

        # propagator: matrix state
        prob = sesolveProblem(H, qeye_like(a), tlist; params = params, progress_bar = Val(false)).prob
        test_analytic_vjp(prob)
    end

    @testset "mesolve" begin
        prob = mesolveProblem(H, ψ0, tlist, c_ops; params = params, progress_bar = Val(false)).prob
        @test prob.f.vjp isa QuantumToolbox.LinearOperatorVJP
        @test prob.f.vjp_p isa QuantumToolbox.LinearOperatorParamVJP
        test_analytic_vjp(prob)

        L = liouvillian(H, c_ops; assume_hermitian = Val(false))
        prob = mesolveProblem(L, ψ0, tlist; params = params, progress_bar = Val(false)).prob
        test_analytic_vjp(prob)

        # whichever VJP the matrix form gets must be correct
        prob = mesolveProblem(H, ψ0, tlist, c_ops; params = params, matrix_form = Val(true), progress_bar = Val(false)).prob
        prob.f.vjp === nothing || test_analytic_vjp(prob)
    end

    @testset "time-independent operators" begin
        prob = sesolveProblem(a' * a + a + a', ψ0, tlist; params = params, progress_bar = Val(false)).prob
        @test prob.f.vjp isa QuantumToolbox.LinearOperatorVJP
        test_analytic_vjp(prob)
        dp = ones(length(params))
        prob.f.vjp_p(dp, prob.u0, prob.u0, params, 0.0)
        @test iszero(dp) # the operator does not depend on the parameters
    end

    @testset "fallback to automatic differentiation" begin
        # structured parameters: only the state VJP is analytic
        H_nt = QobjEvo(a' * a) + QobjEvo(a + a', (p, t) -> p.F * cos(t))
        prob = sesolveProblem(H_nt, ψ0, tlist; params = (F = 0.5,), progress_bar = Val(false)).prob
        @test prob.f.vjp isa QuantumToolbox.LinearOperatorVJP
        @test prob.f.vjp_p === nothing
        test_analytic_vjp(prob)

        # out-of-place problems are left untouched
        prob = sesolveProblem(H, ψ0, tlist; params = params, inplace = Val(false), progress_bar = Val(false)).prob
        @test prob.f.vjp === nothing
        @test prob.f.vjp_p === nothing

        # unsupported operator structures are left untouched
        @test QuantumToolbox._linear_terms(QobjEvo(a, coef_drive).data * QobjEvo(a', coef_drive).data) === nothing
    end
end
