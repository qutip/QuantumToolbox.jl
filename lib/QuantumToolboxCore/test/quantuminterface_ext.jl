using Test
using QuantumToolboxCore
import QuantumInterface

const QI = QuantumInterface

@testset "QuantumInterface Extension" begin
    @test Base.get_extension(QuantumToolboxCore, :QuantumToolboxCoreQuantumInterfaceExt) isa Module

    # the functions exported by QuantumToolboxCore are not changed
    @test tensor === kron
    @test (⊗) === kron
    for f in (:expect, :variance, :fidelity, :entropy_vn, :commutator)
        @test getfield(QuantumToolboxCore, f) !== getfield(QI, f)
    end

    a = destroy(3)
    σx = sigmax()
    ψ = (basis(3, 0) + basis(3, 1)) / √2
    ϕ = basis(2, 1)
    ρ = ket2dm(ψ)
    ρm = maximally_mixed_dm(3)

    @testset "tensor" begin
        @test QI.tensor(a) == a
        @test QI.tensor(a, σx) == kron(a, σx)
        @test QI.:⊗(a, σx) == kron(a, σx)
        @test QI.tensor(a, σx, a') == kron(a, σx, a')
        @test QI.tensor(ψ, ϕ) == kron(ψ, ϕ)
        @test QI.tensor(ψ', ϕ') == kron(ψ', ϕ')
        @test QI.tensor(ψ, ϕ') == kron(ψ, ϕ')
        @test QI.tensor(a, ϕ) == kron(a, ϕ)
        Aevo = @test_logs (:warn, r"lazy tensor") QI.tensor(QobjEvo(a), σx)
        @test Aevo(0.0) == kron(a, σx)
        X_warn = @test_logs (
            :warn,
            "`tensor(A)` or `kron(A)` with `A` is a `Vector` can hurt performance. Try to use `tensor(A...)` or `kron(A...)` instead.",
        ) QI.tensor(fill(σx, 3))
        @test X_warn == kron(σx, σx, σx)
    end

    @testset "expect and variance" begin
        n = a' * a
        @test QI.expect(n, ψ) == expect(n, ψ)
        @test QI.expect(n, ψ') == expect(n, ψ')
        @test QI.expect(n, ρ) == expect(n, ρ)
        @test QI.expect(Hermitian(n), ρ) == expect(Hermitian(n), ρ)
        @test QI.expect([n, a], ψ) == expect([n, a], ψ)
        @test QI.variance(n, ψ) == variance(n, ψ)
        @test QI.variance(n, ρ) == variance(n, ρ)
        @test QI.variance(n, [ψ, ρ]) == variance(n, [ψ, ρ])
    end

    @testset "fidelity and entropy_vn" begin
        @test QI.fidelity(ρ, ρm) == fidelity(ρ, ρm)
        @test QI.fidelity(ρm, ψ) == fidelity(ρm, ψ)
        @test QI.fidelity(ψ, ρm) == fidelity(ψ, ρm)
        @test QI.fidelity(ψ, basis(3, 0)) == fidelity(ψ, basis(3, 0))
        @test QI.entropy_vn(ρm) == entropy_vn(ρm)
        @test QI.entropy_vn(ρm; base = 2) == entropy_vn(ρm; base = 2)
        @test QI.entropy_vn(ψ) == entropy_vn(ψ)
    end

    @testset "commutator" begin
        @test QI.commutator(a, a') == commutator(a, a')
        @test QI.commutator(a, a'; anti = true) == commutator(a, a'; anti = true)
    end

    @testset "Type Inference" begin
        @inferred QI.tensor(a, σx, a')
        @inferred QI.expect(a' * a, ψ)
        @inferred QI.variance(a' * a, ρ)
        @inferred QI.fidelity(ρ, ρm)
        @inferred QI.entropy_vn(ρm)
        @inferred QI.commutator(a, a')
    end

    # mimic a package which exports the functions owned by QuantumInterface.jl (e.g., QuantumOptics.jl and SecondQuantizedAlgebra.jl)
    @eval module OtherPackage
    import QuantumInterface: tensor, ⊗, expect, variance, fidelity, entropy_vn, commutator
    export tensor, ⊗, expect, variance, fidelity, entropy_vn, commutator
    end

    # the name conflicts can be resolved by explicitly importing them from that package
    @eval module UserCode
    using QuantumToolboxCore
    using ..OtherPackage
    using ..OtherPackage: tensor, ⊗, expect, variance, fidelity, entropy_vn, commutator
    f_tensor() = tensor(destroy(3), sigmax())
    f_otimes() = destroy(3) ⊗ sigmax()
    f_expect() = expect(num(3), basis(3, 1))
    f_variance() = variance(num(3), basis(3, 1))
    f_fidelity() = fidelity(basis(3, 1), basis(3, 1))
    f_entropy_vn() = entropy_vn(maximally_mixed_dm(2))
    f_commutator() = commutator(destroy(3), create(3))
    end

    @test UserCode.tensor === QI.tensor
    @test UserCode.f_tensor() == kron(a, σx)
    @test UserCode.f_otimes() == kron(a, σx)
    @test UserCode.f_expect() == expect(num(3), basis(3, 1))
    @test UserCode.f_variance() == variance(num(3), basis(3, 1))
    @test UserCode.f_fidelity() == fidelity(basis(3, 1), basis(3, 1))
    @test UserCode.f_entropy_vn() == entropy_vn(maximally_mixed_dm(2))
    @test UserCode.f_commutator() == commutator(destroy(3), create(3))
end
