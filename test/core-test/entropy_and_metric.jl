@testitem "entropy" begin
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
    @test entropy_relative(ψ, ψ) + 1 ≈ 1
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

    if VERSION >= v"1.11" # eigen is type unstable on v1.10
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
end

@testitem "entanglement and concurrence" begin
    # bell state
    ψb = bell_state(Val(1), Val(0))
    ρb = ket2dm(ψb)
    @test entanglement(ψb, 1) / log(2) ≈ 1
    @test entanglement(ρb, 1) / log(2) ≈ 1
    @test concurrence(ψb) ≈ 1
    @test concurrence(ρb) ≈ 1

    # separable pure state
    ψs = kron(rand_ket(2), rand_ket(2))
    @test entanglement(ψs, 1) + 1 ≈ 1
    @test entanglement(ψs, 2) + 1 ≈ 1
    @test concurrence(ψs) + 1 ≈ 1

    # this only works for "pure" two-qubit states
    ψr = rand_ket((2, 2)) # might be an entangled two-qubit state
    val = concurrence(ψr)
    @test isapprox(val, sqrt(2 * entropy_linear(ptrace(ψr, 1))); atol = 1.0e-5) # √(2 * (1 - Tr(ρA^2)))
    @test isapprox(val, sqrt(2 * entropy_linear(ptrace(ψr, 2))); atol = 1.0e-5) # √(2 * (1 - Tr(ρB^2)))

    @test_throws ArgumentError entanglement(rand_dm((2, 2)), 1)
    @test_throws ArgumentError concurrence(rand_dm((2, 3)))
    @test_throws ArgumentError concurrence(rand_dm(4))

    @testset "Type Stability (entanglement)" begin
        @inferred entanglement(ψb, 1)
        @inferred entanglement(ρb, 1)
    end

    @testset "Type Stability (concurrence)" begin
        @inferred concurrence(ψb)
        @inferred concurrence(ρb)
    end
end

@testitem "trace and Hilbert-Schmidt distance" begin
    ψz0 = basis(2, 0)
    ψz1 = basis(2, 1)
    ρz0 = to_sparse(ket2dm(ψz0))
    ρz1 = to_sparse(ket2dm(ψz1))
    ψx0 = sqrt(0.5) * (basis(2, 0) + basis(2, 1))
    @test tracedist(ψz0, ψx0) ≈ sqrt(0.5)
    @test tracedist(ρz0, ψz1) ≈ 1.0
    @test tracedist(ψz1, ρz0) ≈ 1.0
    @test tracedist(ρz0, ρz1) ≈ 1.0

    ψ = rand_ket(10)
    ϕ = rand_ket(10)
    @test isapprox(tracedist(ψ, ϕ)^2, hilbert_dist(ψ, ϕ) / 2; atol = 1.0e-6)

    @testset "Type Inference (trace distance)" begin
        @inferred tracedist(ψz0, ψx0)
        @inferred tracedist(ρz0, ψz1)
        @inferred tracedist(ψz1, ρz0)
        @inferred tracedist(ρz0, ρz1)
    end

    @testset "Type Inference (Hilbert-Schmidt distance)" begin
        @inferred hilbert_dist(ψz0, ψx0)
        @inferred hilbert_dist(ρz0, ψz1)
        @inferred hilbert_dist(ψz1, ρz0)
        @inferred hilbert_dist(ρz0, ρz1)
    end
end

@testitem "fidelity, Bures metric, and Hellinger distance" begin
    M0 = rand_dm(5)
    ψ1 = rand_ket(5)
    ψ2 = rand_ket(5)
    M1 = ket2dm(ψ1)
    b00 = bell_state(Val(0), Val(0))
    b01 = bell_state(Val(0), Val(1))
    @test isapprox(fidelity(M0, M1), fidelity(ψ1, M0); atol = 1.0e-6)
    @test isapprox(fidelity(ψ1, ψ2), fidelity(ket2dm(ψ1), ket2dm(ψ2)); atol = 1.0e-6)
    @test isapprox(fidelity(b00, b00), 1; atol = 1.0e-6)
    @test isapprox(bures_dist(b00, b00) + 1, 1; atol = 1.0e-6)
    @test isapprox(bures_angle(b00, b00) + 1, 1; atol = 1.0e-6)
    @test isapprox(hellinger_dist(b00, b00) + 1, 1; atol = 1.0e-6)
    @test isapprox(fidelity(b00, b01) + 1, 1; atol = 1.0e-6)
    @test isapprox(bures_dist(b00, b01), √2; atol = 1.0e-6)
    @test isapprox(bures_angle(b00, b01), π / 2; atol = 1.0e-6)
    @test isapprox(hellinger_dist(b00, b01), √2; atol = 1.0e-6)

    # some relations between Bures and Hellinger dintances
    # together with some monotonicity under tensor products
    # [see arXiv:1611.03449 (2017); section 4.2]
    ρA = rand_dm(5)
    ρB = rand_dm(6)
    ρAB = tensor(ρA, ρB)
    σA = rand_dm(5)
    σB = rand_dm(6)
    σAB = tensor(σA, σB)
    d_Bu_A = bures_dist(ρA, σA)
    d_Bu_AB = bures_dist(ρAB, σAB)
    d_He_A = hellinger_dist(ρA, σA)
    d_He_AB = hellinger_dist(ρAB, σAB)
    @test isapprox(fidelity(ρAB, σAB), fidelity(ρA, σA) * fidelity(ρB, σB); atol = 1.0e-6)
    @test d_He_AB >= d_Bu_AB
    @test d_Bu_AB >= d_Bu_A
    @test isapprox(bures_dist(ρAB, tensor(σA, ρB)), d_Bu_A; atol = 1.0e-6)
    @test d_He_AB >= d_He_A
    @test isapprox(hellinger_dist(ρAB, tensor(σA, ρB)), d_He_A; atol = 1.0e-6)

    @testset "Type Inference (fidelity)" begin
        @inferred fidelity(M0, M1)
        @inferred fidelity(ψ1, M0)
        @inferred fidelity(ψ1, ψ2)
    end

    @testset "Type Inference (Hellinger distance)" begin
        @inferred hellinger_dist(M0, M1)
        @inferred hellinger_dist(ψ1, M0)
        @inferred hellinger_dist(ψ1, ψ2)
    end

    @testset "Type Inference (Bures metric)" begin
        @inferred bures_dist(M0, M1)
        @inferred bures_dist(ψ1, M0)
        @inferred bures_dist(ψ1, ψ2)
        @inferred bures_angle(M0, M1)
        @inferred bures_angle(ψ1, M0)
        @inferred bures_angle(ψ1, ψ2)
    end
end
