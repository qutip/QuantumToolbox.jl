using Test
using QuantumToolbox
using CUDA
import SparseArrays: SparseMatrixCSC

@testset "CUDA (example in README)" verbose = true begin
    N = 20
    ω64 = 1.0    # Float64
    ω32 = 1.0f0  # Float32
    γ64 = 0.1    # Float64
    γ32 = 0.1f0  # Float32
    tlist = range(0, 10, 100)

    ## calculate by CPU (with 64-bit)
    a_cpu64 = destroy(N)
    ψ0_cpu64 = fock(N, 3)
    H_cpu64 = ω64 * a_cpu64' * a_cpu64
    c_ops_cpu64 = [sqrt(γ64) * a_cpu64]
    sol_cpu64 = mesolve(H_cpu64, ψ0_cpu64, tlist, c_ops_cpu64, e_ops = [a_cpu64' * a_cpu64], progress_bar = Val(false))

    ## calculate by CPU (with 32-bit)
    a_cpu32 = destroy(ComplexF32, N)
    ψ0_cpu32 = fock(ComplexF32, N, 3)
    H_cpu32 = ω32 * a_cpu32' * a_cpu32
    c_ops_cpu32 = [sqrt(γ32) * a_cpu32]
    sol_cpu32 = mesolve(H_cpu32, ψ0_cpu32, tlist, c_ops_cpu32, e_ops = [a_cpu32' * a_cpu32], progress_bar = Val(false))

    ## calculate by GPU (with 64-bit)
    a_gpu64 = cu(destroy(N))
    ψ0_gpu64 = cu(fock(N, 3))
    H_gpu64 = ω64 * a_gpu64' * a_gpu64
    c_ops_gpu64 = [sqrt(γ64) * a_gpu64]
    sol_gpu64 = mesolve(H_gpu64, ψ0_gpu64, tlist, c_ops_gpu64, e_ops = [a_gpu64' * a_gpu64], progress_bar = Val(false))
    # matrix_form = Val(true) case
    sol_gpu64_mat = mesolve(H_gpu64, ψ0_gpu64, tlist, c_ops_gpu64, e_ops = [a_gpu64' * a_gpu64], progress_bar = Val(false), matrix_form = Val(true))

    ## calculate by GPU (with 32-bit)
    a_gpu32 = cu(destroy(N), word_size = 32)
    ψ0_gpu32 = cu(fock(N, 3), word_size = 32)
    H_gpu32 = ω32 * a_gpu32' * a_gpu32
    c_ops_gpu32 = [sqrt(γ32) * a_gpu32]
    sol_gpu32 = mesolve(H_gpu32, ψ0_gpu32, tlist, c_ops_gpu32, e_ops = [a_gpu32' * a_gpu32], progress_bar = Val(false))
    # matrix_form = Val(true) case
    sol_gpu32_mat = mesolve(H_gpu32, ψ0_gpu32, tlist, c_ops_gpu32, e_ops = [a_gpu32' * a_gpu32], progress_bar = Val(false), matrix_form = Val(true))

    L_cpu64 = liouvillian(H_cpu64, c_ops_cpu64)
    L_gpu64 = liouvillian(H_gpu64, c_ops_gpu64)

    @test SparseMatrixCSC(L_gpu64.data) ≈ L_cpu64.data

    @test all([isapprox(sol_cpu64.expect[i], sol_gpu64.expect[i]) for i in 1:length(tlist)])
    @test all([isapprox(sol_cpu32.expect[i], sol_gpu32.expect[i]; atol = 1.0f-6) for i in 1:length(tlist)])
    @test all([isapprox(sol_cpu64.expect[i], sol_gpu64_mat.expect[i]) for i in 1:length(tlist)])
    @test all([isapprox(sol_cpu32.expect[i], sol_gpu32_mat.expect[i]; atol = 1.0f-6) for i in 1:length(tlist)])
end
