@testset "CUDA Extension" verbose = true begin
    # Test that scalar indexing is disallowed
    @test_throws ErrorException CUDA.rand(1)[1]

    ψdi = Qobj(Int64[1, 0])
    ψdf = Qobj(Float64[1, 0])
    ψdc = Qobj(ComplexF64[1, 0])
    ψsi = to_sparse(ψdi)
    ψsf = to_sparse(ψdf)
    ψsc = to_sparse(ψdc)

    Xdi = Qobj(Int64[0 1; 1 0])
    Xdf = Qobj(Float64[0 1; 1 0])
    Xdc = Qobj(ComplexF64[0 1; 1 0])
    Xsi = to_sparse(Xdi)
    Xsf = to_sparse(Xdf)
    Xsc = to_sparse(Xdc)

    @test_throws DomainError cu(ψdi; word_size = 16)

    # type conversion of CUDA dense arrays
    @test typeof(cu(ψdi; word_size = 64).data) == typeof(CuArray(ψdi).data) <: CuArray{Int64,1}
    @test typeof(cu(ψdi; word_size = 32).data) == typeof(CuArray{Int32}(ψdi).data) <: CuArray{Int32,1}
    @test typeof(cu(ψdf; word_size = 64).data) == typeof(CuArray(ψdf).data) <: CuArray{Float64,1}
    @test typeof(cu(ψdf; word_size = 32).data) == typeof(CuArray{Float32}(ψdf).data) <: CuArray{Float32,1}
    @test typeof(cu(ψdc; word_size = 64).data) == typeof(CuArray(ψdc).data) <: CuArray{ComplexF64,1}
    @test typeof(cu(ψdc; word_size = 32).data) == typeof(CuArray{ComplexF32}(ψdc).data) <: CuArray{ComplexF32,1}
    @test typeof(cu(Xdi; word_size = 64).data) == typeof(CuArray(Xdi).data) <: CuArray{Int64,2}
    @test typeof(cu(Xdi; word_size = 32).data) == typeof(CuArray{Int32}(Xdi).data) <: CuArray{Int32,2}
    @test typeof(cu(Xdf; word_size = 64).data) == typeof(CuArray(Xdf).data) <: CuArray{Float64,2}
    @test typeof(cu(Xdf; word_size = 32).data) == typeof(CuArray{Float32}(Xdf).data) <: CuArray{Float32,2}
    @test typeof(cu(Xdc; word_size = 64).data) == typeof(CuArray(Xdc).data) <: CuArray{ComplexF64,2}
    @test typeof(cu(Xdc; word_size = 32).data) == typeof(CuArray{ComplexF32}(Xdc).data) <: CuArray{ComplexF32,2}

    # type conversion of CUDA sparse arrays
    @test typeof(cu(ψsi; word_size = 64).data) == typeof(CuSparseVector(ψsi).data) == CuSparseVector{Int64,Int32}
    @test typeof(cu(ψsi; word_size = 32).data) == typeof(CuSparseVector{Int32}(ψsi).data) == CuSparseVector{Int32,Int32}
    @test typeof(cu(ψsf; word_size = 64).data) == typeof(CuSparseVector(ψsf).data) == CuSparseVector{Float64,Int32}
    @test typeof(cu(ψsf; word_size = 32).data) ==
          typeof(CuSparseVector{Float32}(ψsf).data) ==
          CuSparseVector{Float32,Int32}
    @test typeof(cu(ψsc; word_size = 64).data) == typeof(CuSparseVector(ψsc).data) == CuSparseVector{ComplexF64,Int32}
    @test typeof(cu(ψsc; word_size = 32).data) ==
          typeof(CuSparseVector{ComplexF32}(ψsc).data) ==
          CuSparseVector{ComplexF32,Int32}
    @test typeof(cu(Xsi; word_size = 64).data) == typeof(CuSparseMatrixCSC(Xsi).data) == CuSparseMatrixCSC{Int64,Int32}
    @test typeof(cu(Xsi; word_size = 32).data) ==
          typeof(CuSparseMatrixCSC{Int32}(Xsi).data) ==
          CuSparseMatrixCSC{Int32,Int32}
    @test typeof(cu(Xsf; word_size = 64).data) ==
          typeof(CuSparseMatrixCSC(Xsf).data) ==
          CuSparseMatrixCSC{Float64,Int32}
    @test typeof(cu(Xsf; word_size = 32).data) ==
          typeof(CuSparseMatrixCSC{Float32}(Xsf).data) ==
          CuSparseMatrixCSC{Float32,Int32}
    @test typeof(cu(Xsc; word_size = 64).data) ==
          typeof(CuSparseMatrixCSC(Xsc).data) ==
          CuSparseMatrixCSC{ComplexF64,Int32}
    @test typeof(cu(Xsc; word_size = 32).data) ==
          typeof(CuSparseMatrixCSC{ComplexF32}(Xsc).data) ==
          CuSparseMatrixCSC{ComplexF32,Int32}
    @test typeof(CuSparseMatrixCSR(Xsi).data) == CuSparseMatrixCSR{Int64,Int32}
    @test typeof(CuSparseMatrixCSR{Int32}(Xsi).data) == CuSparseMatrixCSR{Int32,Int32}
    @test typeof(CuSparseMatrixCSR(Xsf).data) == CuSparseMatrixCSR{Float64,Int32}
    @test typeof(CuSparseMatrixCSR{Float32}(Xsf).data) == CuSparseMatrixCSR{Float32,Int32}
    @test typeof(CuSparseMatrixCSR(Xsc).data) == CuSparseMatrixCSR{ComplexF64,Int32}
    @test typeof(CuSparseMatrixCSR{ComplexF32}(Xsc).data) == CuSparseMatrixCSR{ComplexF32,Int32}

    # type conversion of CUDA Diagonal arrays
    @test cu(qeye(10), word_size = Val(32)).data isa Diagonal{ComplexF32,<:CuVector{ComplexF32}}
    @test cu(qeye(10), word_size = Val(64)).data isa Diagonal{ComplexF64,<:CuVector{ComplexF64}}

    # Sparse To Dense
    # @test to_dense(cu(ψsi; word_size = 64)).data isa CuVector{Int64} # TODO: Fix this in CUDA.jl
    @test to_dense(cu(ψsf; word_size = 64)).data isa CuVector{Float64}
    @test to_dense(cu(ψsc; word_size = 64)).data isa CuVector{ComplexF64}
    # @test to_dense(cu(Xsi; word_size = 64)).data isa CuMatrix{Int64} # TODO: Fix this in CUDA.jl
    @test to_dense(cu(Xsf; word_size = 64)).data isa CuMatrix{Float64}
    @test to_dense(cu(Xsc; word_size = 64)).data isa CuMatrix{ComplexF64}

    # @test to_dense(Int32, cu(ψsf; word_size = 64)).data isa CuVector{Int32} # TODO: Fix this in CUDA.jl
    # @test to_dense(Float32, cu(ψsf; word_size = 64)).data isa CuVector{Float32} # TODO: Fix this in CUDA.jl
    # @test to_dense(ComplexF32, cu(ψsf; word_size = 64)).data isa CuVector{ComplexF32} # TODO: Fix this in CUDA.jl
    # @test to_dense(Int64, cu(Xsf; word_size = 32)).data isa CuMatrix{Int64} # TODO: Fix this in CUDA.jl
    # @test to_dense(Float64, cu(Xsf; word_size = 32)).data isa CuMatrix{Float64} # TODO: Fix this in CUDA.jl
    # @test to_dense(ComplexF64, cu(Xsf; word_size = 32)).data isa CuMatrix{ComplexF64} # TODO: Fix this in CUDA.jl

    # brief example in README and documentation
    N = 20
    ω64 = 1.0    # Float64
    ω32 = 1.0f0  # Float32
    γ64 = 0.1    # Float64
    γ32 = 0.1f0  # Float32
    tlist = range(0, 10, 100)

    ## calculate by CPU
    a_cpu = destroy(N)
    ψ0_cpu = fock(N, 3)
    H_cpu = ω64 * a_cpu' * a_cpu
    sol_cpu = mesolve(H_cpu, ψ0_cpu, tlist, [sqrt(γ64) * a_cpu], e_ops = [a_cpu' * a_cpu], progress_bar = Val(false))

    ## calculate by GPU (with 64-bit)
    a_gpu64 = cu(destroy(N))
    ψ0_gpu64 = cu(fock(N, 3))
    H_gpu64 = ω64 * a_gpu64' * a_gpu64
    sol_gpu64 = mesolve(
        H_gpu64,
        ψ0_gpu64,
        tlist,
        [sqrt(γ64) * a_gpu64],
        e_ops = [a_gpu64' * a_gpu64],
        progress_bar = Val(false),
    )

    ## calculate by GPU (with 32-bit)
    a_gpu32 = cu(destroy(N), word_size = 32)
    ψ0_gpu32 = cu(fock(N, 3), word_size = 32)
    H_gpu32 = ω32 * a_gpu32' * a_gpu32
    sol_gpu32 = mesolve(
        H_gpu32,
        ψ0_gpu32,
        tlist,
        [sqrt(γ32) * a_gpu32],
        e_ops = [a_gpu32' * a_gpu32],
        progress_bar = Val(false),
    )

    @test all([isapprox(sol_cpu.expect[i], sol_gpu64.expect[i]) for i in 1:length(tlist)])
    @test all([isapprox(sol_cpu.expect[i], sol_gpu32.expect[i]; atol = 1e-6) for i in 1:length(tlist)])
end

@testset "CUDA steadystate" begin
    N = 50
    Δ = 0.01
    F = 0.1
    γ = 0.1
    nth = 2

    a = destroy(N)
    H = Δ * a' * a + F * (a + a')
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a']

    ρ_ss_cpu = steadystate(H, c_ops)

    H_gpu_csc = cu(H)
    c_ops_gpu_csc = [cu(c_op) for c_op in c_ops]
    ρ_ss_gpu_csc = steadystate(H_gpu_csc, c_ops_gpu_csc, solver = SteadyStateLinearSolver())

    H_gpu_csr = CuSparseMatrixCSR(H_gpu_csc)
    c_ops_gpu_csr = [CuSparseMatrixCSR(c_op) for c_op in c_ops_gpu_csc]
    ρ_ss_gpu_csr = steadystate(H_gpu_csr, c_ops_gpu_csr, solver = SteadyStateLinearSolver())

    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csc.data) atol = 1e-8 * length(ρ_ss_cpu)
    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csr.data) atol = 1e-8 * length(ρ_ss_cpu)
end

@testset "CUDA ptrace" begin
    using Test
    using QuantumToolbox
    import StaticArraysCore: SVector

    using CUDA
    using CUDA.CUSPARSE

    g = fock(2, 1)
    e = fock(2, 0)
    α = sqrt(0.7)
    β = sqrt(0.3) * 1im
    ψ = α * kron(g, e) + β * kron(e, g) |> cu

    ρ1 = ptrace(ψ, 1)
    ρ2 = ptrace(ψ, 2)
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ψ_d = ψ'

    ρ1 = ptrace(ψ_d, 1)
    ρ2 = ptrace(ψ_d, 2)
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ρ = ket2dm(ψ)
    ρ1 = ptrace(ρ, 1)
    ρ2 = ptrace(ρ, 2)
    @test ρ.data isa CuArray
    @test ρ1.data isa CuArray
    @test ρ2.data isa CuArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ψ1 = normalize(g + 1im * e)
    ψ2 = normalize(g + e)
    ρ1 = ket2dm(ψ1)
    ρ2 = ket2dm(ψ2)
    ρ = kron(ρ1, ρ2) |> cu
    ρ1_ptr = ptrace(ρ, 1)
    ρ2_ptr = ptrace(ρ, 2)
    @test ρ1_ptr.data isa CuArray
    @test ρ2_ptr.data isa CuArray
    @test ρ1.data ≈ Array(ρ1_ptr.data) atol = 1e-10
    @test ρ2.data ≈ Array(ρ2_ptr.data) atol = 1e-10

    ψlist = [rand_ket(2), rand_ket(3), rand_ket(4), rand_ket(5)]
    ρlist = [rand_dm(2), rand_dm(3), rand_dm(4), rand_dm(5)]
    ψtotal = tensor(ψlist...) |> cu
    ρtotal = tensor(ρlist...) |> cu
    sel_tests = [
        SVector{0,Int}(),
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
    ]
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
