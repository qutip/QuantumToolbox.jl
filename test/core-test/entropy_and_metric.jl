@testset "entropy" begin
    base = 2
    λ = rand()
    ψ = rand_ket(10)
    ρ1 = rand_dm(10)
    ρ2 = rand_dm(10)
    σ1 = rand_dm(10)
    σ2 = rand_dm(10)

    dims = (2, 3)
    ρAB = rand_dm((dims..., dims...))
    selA = (1, 2)
    selB = (3, 4)
    ρA = ptrace(ρAB, selA)
    ρB = ptrace(ρAB, selB)
    nA = nB = prod(dims)
    IA = qeye(nA, dims = dims)
    IB = qeye(nB, dims = dims)

    # quantum relative entropy
    @test entropy_relative(ρ1, ψ) == Inf
    @test entropy_relative(ρ1, rand_dm(10, rank = 9)) == Inf
    @test entropy_relative(ψ, ψ) ≈ 0
    @test entropy_relative(λ * ρ1 + (1 - λ) * ρ2, λ * σ1 + (1 - λ) * σ2) <=
          λ * entropy_relative(ρ1, σ1) + (1 - λ) * entropy_relative(ρ2, σ2) # joint convexity

    # relations between different entropies
    @test entropy_relative(ρA, IA / nA) ≈ log(nA) - entropy_vn(ρA)
    @test entropy_relative(ρB, IB / nB; base = base) ≈ log(base, nB) - entropy_vn(ρB; base = base)
    @test entropy_relative(ρAB, tensor(ρA, ρB)) ≈ entropy_mutual(ρAB, selA, selB)
    @test entropy_relative(ρAB, tensor(ρA, ρB)) ≈ entropy_mutual(ρAB, selA, selB)
    @test entropy_relative(ρAB, tensor(ρA, IB / nB)) ≈ log(nB) - entropy_conditional(ρAB, selA)
    @test entropy_linear(ρ1) == 1 - purity(ρ1)

    ρ_wrong = Qobj(rand(ComplexF64, 10, 10))
    @test_throws ErrorException entropy_relative(ρ1, ρ_wrong)
    @test_throws ErrorException entropy_relative(ρ_wrong, ρ1)
    @test_throws ArgumentError entropy_mutual(ρAB, 1, 3)
    @test_throws ArgumentError entropy_mutual(ρAB, 1, (3, 4))
    @test_throws ArgumentError entropy_mutual(ρAB, (1, 2), 3)
    @test_throws ArgumentError entropy_mutual(ρAB, (1, 2), (1, 3))

    @testset "Type Stability (entropy)" begin
        @inferred entropy_vn(ρ1)
        @inferred entropy_vn(ρ1, base = base)
        @inferred entropy_relative(ρ1, ψ)
        @inferred entropy_relative(ρ1, σ1, base = base)
        @inferred entropy_linear(ρ1)
        @inferred entropy_mutual(ρAB, selA, selB)
        @inferred entropy_conditional(ρAB, selA)
    end
end

@testset "Entanglement" begin
    g = fock(2, 1)
    e = fock(2, 0)
    state = normalize(kron(g, e) + kron(e, g))
    rho = state * state'
    @test entanglement(state, 1) / log(2) ≈ 1
    @test entanglement(rho, 1) / log(2) ≈ 1

    @testset "Type Stability (entanglement)" begin
        @inferred entanglement(state, 1)
        @inferred entanglement(rho, 1)
    end
end

@testset "trace distance" begin
    ψz0 = basis(2, 0)
    ψz1 = basis(2, 1)
    ρz0 = to_sparse(ket2dm(ψz0))
    ρz1 = to_sparse(ket2dm(ψz1))
    ψx0 = sqrt(0.5) * (basis(2, 0) + basis(2, 1))
    @test tracedist(ψz0, ψx0) ≈ sqrt(0.5)
    @test tracedist(ρz0, ψz1) ≈ 1.0
    @test tracedist(ψz1, ρz0) ≈ 1.0
    @test tracedist(ρz0, ρz1) ≈ 1.0

    @testset "Type Inference (trace distance)" begin
        @inferred tracedist(ψz0, ψx0)
        @inferred tracedist(ρz0, ψz1)
        @inferred tracedist(ψz1, ρz0)
        @inferred tracedist(ρz0, ρz1)
    end
end

@testset "fidelity" begin
    M = sprand(ComplexF64, 5, 5, 0.5)
    M0 = Qobj(M * M')
    ψ1 = Qobj(rand(ComplexF64, 5))
    ψ2 = Qobj(rand(ComplexF64, 5))
    M1 = ψ1 * ψ1'
    @test isapprox(fidelity(M0, M1), fidelity(ψ1, M0); atol = 1e-6)
    @test isapprox(fidelity(ψ1, ψ2), fidelity(ket2dm(ψ1), ket2dm(ψ2)); atol = 1e-6)

    @testset "Type Inference (fidelity)" begin
        @inferred fidelity(M0, M1)
        @inferred fidelity(ψ1, M0)
        @inferred fidelity(ψ1, ψ2)
    end
end
