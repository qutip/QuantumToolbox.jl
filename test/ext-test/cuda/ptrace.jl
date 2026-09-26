using QuantumToolbox
using CUDA

@testset "CUDA (ptrace)" begin
    g = fock(2, 1)
    e = fock(2, 0)
    α = sqrt(0.7)
    β = sqrt(0.3) * 1im
    ψ = α * kron(g, e) + β * kron(e, g) |> cu

    ρ1 = ptrace(ψ, 1)
    ρ2 = ptrace(ψ, 2)
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1.0e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1.0e-10

    ψ_d = ψ'

    ρ1 = ptrace(ψ_d, 1)
    ρ2 = ptrace(ψ_d, 2)
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1.0e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1.0e-10

    ρ = ket2dm(ψ)
    ρ1 = ptrace(ρ, 1)
    ρ2 = ptrace(ρ, 2)
    @test ρ.data isa CuArray
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1.0e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1.0e-10

    ψ1 = normalize(g + 1im * e)
    ψ2 = normalize(g + e)
    ρ1 = ket2dm(ψ1)
    ρ2 = ket2dm(ψ2)
    ρ = kron(ρ1, ρ2) |> cu
    ρ1_ptr = ptrace(ρ, 1)
    ρ2_ptr = ptrace(ρ, 2)
    @test ρ1_ptr.data isa CuArray
    @test ρ2_ptr.data isa CuArray
    @test ρ1.data ≈ Array(ρ1_ptr.data) atol = 1.0e-10
    @test ρ2.data ≈ Array(ρ2_ptr.data) atol = 1.0e-10

    ψlist = [rand_ket(2), rand_ket(3), rand_ket(4), rand_ket(5)]
    ρlist = [rand_dm(2), rand_dm(3), rand_dm(4), rand_dm(5)]
    ψtotal = tensor(ψlist...) |> cu
    ρtotal = tensor(ρlist...) |> cu
    sel_tests = (
        1,
        2,
        3,
        4,
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
        (1, 2, 3),
        (1, 2, 4),
        (1, 3, 4),
        (2, 3, 4),
        (1, 2, 3, 4),
    )
    for sel in sel_tests
        if length(sel) == 0
            @test ptrace(ψtotal, sel) ≈ 1.0
            @test ptrace(ρtotal, sel) ≈ 1.0
        else
            @test ptrace(ψtotal, sel) ≈ cu(tensor([ket2dm(ψlist[i]) for i in sel]...))
            @test ptrace(ρtotal, sel) ≈ cu(tensor([ρlist[i] for i in sel]...))
        end
    end
    @test ptrace(ψtotal, (1, 3, 4)) ≈ ptrace(ψtotal, (4, 3, 1)) # check sort of sel
    @test ptrace(ρtotal, (1, 3, 4)) ≈ ptrace(ρtotal, (3, 1, 4)) # check sort of sel

    @testset "Type Inference (ptrace)" begin
        @inferred ptrace(ρ, 1)
        @inferred ptrace(ρ, 2)
        @inferred ptrace(ψ_d, 1)
        @inferred ptrace(ψ_d, 2)
        @inferred ptrace(ψ, 1)
        @inferred ptrace(ψ, 2)
    end
end
