@testset "AMDGPU Extension" verbose = true begin
    # Test that scalar indexing is disallowed
    @test_throws ErrorException AMDGPU.rand(1)[1]

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

    @test_throws DomainError roc(ψdi; word_size = 16)

    # type conversion of AMDGPU dense arrays
    @test typeof(roc(ψdi; word_size = 64).data) == typeof(ROCArray(ψdi).data) <: ROCArray{Int64,1}
    @test typeof(roc(ψdi; word_size = 32).data) == typeof(ROCArray{Int32}(ψdi).data) <: ROCArray{Int32,1}
    @test typeof(roc(ψdf; word_size = 64).data) == typeof(ROCArray(ψdf).data) <: ROCArray{Float64,1}
    @test typeof(roc(ψdf; word_size = 32).data) == typeof(ROCArray{Float32}(ψdf).data) <: ROCArray{Float32,1}
    @test typeof(roc(ψdc; word_size = 64).data) == typeof(ROCArray(ψdc).data) <: ROCArray{ComplexF64,1}
    @test typeof(roc(ψdc; word_size = 32).data) == typeof(ROCArray{ComplexF32}(ψdc).data) <: ROCArray{ComplexF32,1}
    @test typeof(roc(Xdi; word_size = 64).data) == typeof(ROCArray(Xdi).data) <: ROCArray{Int64,2}
    @test typeof(roc(Xdi; word_size = 32).data) == typeof(ROCArray{Int32}(Xdi).data) <: ROCArray{Int32,2}
    @test typeof(roc(Xdf; word_size = 64).data) == typeof(ROCArray(Xdf).data) <: ROCArray{Float64,2}
    @test typeof(roc(Xdf; word_size = 32).data) == typeof(ROCArray{Float32}(Xdf).data) <: ROCArray{Float32,2}
    @test typeof(roc(Xdc; word_size = 64).data) == typeof(ROCArray(Xdc).data) <: ROCArray{ComplexF64,2}
    @test typeof(roc(Xdc; word_size = 32).data) == typeof(ROCArray{ComplexF32}(Xdc).data) <: ROCArray{ComplexF32,2}

    # type conversion of AMDGPU sparse arrays
    @test typeof(roc(ψsi; word_size = 64).data) == typeof(ROCSparseVector(ψsi).data) == ROCSparseVector{Int64,Int32}
    @test typeof(roc(ψsi; word_size = 32).data) == typeof(ROCSparseVector{Int32}(ψsi).data) == ROCSparseVector{Int32,Int32}
    @test typeof(roc(ψsf; word_size = 64).data) == typeof(ROCSparseVector(ψsf).data) == ROCSparseVector{Float64,Int32}
    @test typeof(roc(ψsf; word_size = 32).data) ==
          typeof(ROCSparseVector{Float32}(ψsf).data) ==
          ROCSparseVector{Float32,Int32}
    @test typeof(roc(ψsc; word_size = 64).data) == typeof(ROCSparseVector(ψsc).data) == ROCSparseVector{ComplexF64,Int32}
    @test typeof(roc(ψsc; word_size = 32).data) ==
          typeof(ROCSparseVector{ComplexF32}(ψsc).data) ==
          ROCSparseVector{ComplexF32,Int32}
    @test typeof(roc(Xsi; word_size = 64).data) == typeof(ROCSparseMatrixCSC(Xsi).data) == ROCSparseMatrixCSC{Int64,Int32}
    @test typeof(roc(Xsi; word_size = 32).data) ==
          typeof(ROCSparseMatrixCSC{Int32}(Xsi).data) ==
          ROCSparseMatrixCSC{Int32,Int32}
    @test typeof(roc(Xsf; word_size = 64).data) ==
          typeof(ROCSparseMatrixCSC(Xsf).data) ==
          ROCSparseMatrixCSC{Float64,Int32}
    @test typeof(roc(Xsf; word_size = 32).data) ==
          typeof(ROCSparseMatrixCSC{Float32}(Xsf).data) ==
          ROCSparseMatrixCSC{Float32,Int32}
    @test typeof(roc(Xsc; word_size = 64).data) ==
          typeof(ROCSparseMatrixCSC(Xsc).data) ==
          ROCSparseMatrixCSC{ComplexF64,Int32}
    @test typeof(roc(Xsc; word_size = 32).data) ==
          typeof(ROCSparseMatrixCSC{ComplexF32}(Xsc).data) ==
          ROCSparseMatrixCSC{ComplexF32,Int32}
    @test typeof(ROCSparseMatrixCSR(Xsi).data) == ROCSparseMatrixCSR{Int64,Int32}
    @test typeof(ROCSparseMatrixCSR{Int32}(Xsi).data) == ROCSparseMatrixCSR{Int32,Int32}
    @test typeof(ROCSparseMatrixCSR(Xsf).data) == ROCSparseMatrixCSR{Float64,Int32}
    @test typeof(ROCSparseMatrixCSR{Float32}(Xsf).data) == ROCSparseMatrixCSR{Float32,Int32}
    @test typeof(ROCSparseMatrixCSR(Xsc).data) == ROCSparseMatrixCSR{ComplexF64,Int32}
    @test typeof(ROCSparseMatrixCSR{ComplexF32}(Xsc).data) == ROCSparseMatrixCSR{ComplexF32,Int32}

    # type conversion of AMDGPU Diagonal arrays
    @test roc(qeye(10), word_size = Val(32)).data isa Diagonal{ComplexF32,<:ROCVector{ComplexF32}}
    @test roc(qeye(10), word_size = Val(64)).data isa Diagonal{ComplexF64,<:ROCVector{ComplexF64}}

    # Sparse To Dense
    # @test to_dense(roc(ψsi; word_size = 64)).data isa ROCVector{Int64} # TODO: Fix this in AMDGPU.jl
    @test to_dense(roc(ψsf; word_size = 64)).data isa ROCVector{Float64}
    @test to_dense(roc(ψsc; word_size = 64)).data isa ROCVector{ComplexF64}
    # @test to_dense(roc(Xsi; word_size = 64)).data isa ROCMatrix{Int64} # TODO: Fix this in AMDGPU.jl
    @test to_dense(roc(Xsf; word_size = 64)).data isa ROCMatrix{Float64}
    @test to_dense(roc(Xsc; word_size = 64)).data isa ROCMatrix{ComplexF64}

    # @test to_dense(Int32, roc(ψsf; word_size = 64)).data isa ROCVector{Int32} # TODO: Fix this in AMDGPU.jl
    # @test to_dense(Float32, roc(ψsf; word_size = 64)).data isa ROCVector{Float32} # TODO: Fix this in AMDGPU.jl
    # @test to_dense(ComplexF32, roc(ψsf; word_size = 64)).data isa ROCVector{ComplexF32} # TODO: Fix this in AMDGPU.jl
    # @test to_dense(Int64, roc(Xsf; word_size = 32)).data isa ROCMatrix{Int64} # TODO: Fix this in AMDGPU.jl
    # @test to_dense(Float64, roc(Xsf; word_size = 32)).data isa ROCMatrix{Float64} # TODO: Fix this in AMDGPU.jl
    # @test to_dense(ComplexF64, roc(Xsf; word_size = 32)).data isa ROCMatrix{ComplexF64} # TODO: Fix this in AMDGPU.jl

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
    a_gpu64 = roc(destroy(N))
    ψ0_gpu64 = roc(fock(N, 3))
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
    a_gpu32 = roc(destroy(N), word_size = 32)
    ψ0_gpu32 = roc(fock(N, 3), word_size = 32)
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

@testset "AMDGPU steadystate" begin
    N = 50
    Δ = 0.01
    F = 0.1
    γ = 0.1
    nth = 2

    a = destroy(N)
    H = Δ * a' * a + F * (a + a')
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a']

    ρ_ss_cpu = steadystate(H, c_ops)

    H_gpu_csc = roc(H)
    c_ops_gpu_csc = [roc(c_op) for c_op in c_ops]
    ρ_ss_gpu_csc = steadystate(H_gpu_csc, c_ops_gpu_csc, solver = SteadyStateLinearSolver())

    H_gpu_csr = ROCSparseMatrixCSR(H_gpu_csc)
    c_ops_gpu_csr = [ROCSparseMatrixCSR(c_op) for c_op in c_ops_gpu_csc]
    ρ_ss_gpu_csr = steadystate(H_gpu_csr, c_ops_gpu_csr, solver = SteadyStateLinearSolver())

    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csc.data) atol = 1e-8 * length(ρ_ss_cpu)
    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csr.data) atol = 1e-8 * length(ρ_ss_cpu)
end

@testset "AMDGPU spectrum" begin
    N = 10
    a = roc(destroy(N))
    H = a' * a
    c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']
    solver = Lanczos(steadystate_solver = SteadyStateLinearSolver())

    ω_l = range(0, 3, length = 1000)
    spec = spectrum(H, ω_l, c_ops, a', a; solver = solver)

    spec = collect(spec)
    spec = spec ./ maximum(spec)

    test_func = maximum(real.(spec)) * (0.1 / 2)^2 ./ ((ω_l .- 1) .^ 2 .+ (0.1 / 2)^2)
    idxs = test_func .> 0.05
    @test sum(abs2.(spec[idxs] .- test_func[idxs])) / sum(abs2.(test_func[idxs])) < 0.01

    # TODO: Fix this
    # @testset "Type Inference spectrum" begin
    #     @inferred spectrum(H, ω_l, c_ops, a', a; solver = solver)
    # end
end

@testset "AMDGPU ptrace" begin
    g = fock(2, 1)
    e = fock(2, 0)
    α = sqrt(0.7)
    β = sqrt(0.3) * 1im
    ψ = α * kron(g, e) + β * kron(e, g) |> roc

    ρ1 = ptrace(ψ, 1)
    ρ2 = ptrace(ψ, 2)
    @test ρ1.data isa ROCArray
    @test ρ2.data isa ROCArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ψ_d = ψ'

    ρ1 = ptrace(ψ_d, 1)
    ρ2 = ptrace(ψ_d, 2)
    @test ρ1.data isa ROCArray
    @test ρ2.data isa ROCArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ρ = ket2dm(ψ)
    ρ1 = ptrace(ρ, 1)
    ρ2 = ptrace(ρ, 2)
    @test ρ.data isa ROCArray
    @test ρ1.data isa ROCArray
    @test ρ2.data isa ROCArray
    @test Array(ρ1.data) ≈ [0.3 0.0; 0.0 0.7] atol = 1e-10
    @test Array(ρ2.data) ≈ [0.7 0.0; 0.0 0.3] atol = 1e-10

    ψ1 = normalize(g + 1im * e)
    ψ2 = normalize(g + e)
    ρ1 = ket2dm(ψ1)
    ρ2 = ket2dm(ψ2)
    ρ = kron(ρ1, ρ2) |> roc
    ρ1_ptr = ptrace(ρ, 1)
    ρ2_ptr = ptrace(ρ, 2)
    @test ρ1_ptr.data isa ROCArray
    @test ρ2_ptr.data isa ROCArray
    @test ρ1.data ≈ Array(ρ1_ptr.data) atol = 1e-10
    @test ρ2.data ≈ Array(ρ2_ptr.data) atol = 1e-10

    ψlist = [rand_ket(2), rand_ket(3), rand_ket(4), rand_ket(5)]
    ρlist = [rand_dm(2), rand_dm(3), rand_dm(4), rand_dm(5)]
    ψtotal = tensor(ψlist...) |> roc
    ρtotal = tensor(ρlist...) |> roc
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
            @test ptrace(ψtotal, sel) ≈ roc(tensor([ket2dm(ψlist[i]) for i in sel]...))
            @test ptrace(ρtotal, sel) ≈ roc(tensor([ρlist[i] for i in sel]...))
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

@testset "AMDGPU eigsolve" begin
    N = 30
    Δ = 0.5
    U = 0.1
    κ = 0.1
    F = 0.5

    a = destroy(N)
    H = Δ * a' * a + U / 2 * a' * a' * a * a + F * (a + a')

    c_ops = [sqrt(κ) * a]

    L = liouvillian(H, c_ops)
    L_gpu = ROCSparseMatrixCSR(L)

    vals_cpu, vecs_cpu = eigenstates(L; sparse = true, sigma = 0.01, eigvals = 4, krylovdim = 30)
    vals_gpu, vecs_gpu = eigenstates(
        L_gpu;
        sparse = true,
        sigma = 0.01,
        eigvals = 4,
        krylovdim = 30,
        solver = LUFactorization(),
        v0 = AMDGPU.rand(ComplexF64, size(L_gpu, 1)),
    )

    @test vals_cpu ≈ vals_gpu atol = 1e-8
    @test all(zip(vecs_cpu, vecs_gpu)) do (v_cpu, v_gpu)
        return isapprox(abs(dot(v_cpu.data, Array(v_gpu.data))), 1; atol = 1e-8)
    end
end
