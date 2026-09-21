using Test
using QuantumToolbox
using CUDA
using CUDA.cuSPARSE
using CUDSS

@testset "CUDA (steadystate)" begin
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
    ρ_ss_gpu_csr = steadystate(H_gpu_csr, c_ops_gpu_csr) # Direct solver using CUDSS

    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csc.data) atol = 1.0e-8 * length(ρ_ss_cpu)
    @test ρ_ss_cpu.data ≈ Array(ρ_ss_gpu_csr.data) atol = 1.0e-8 * length(ρ_ss_cpu)
end