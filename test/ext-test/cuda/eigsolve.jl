using Test
using QuantumToolbox
using CUDA
using CUDA.cuSPARSE
using CUDSS

@testset "CUDA (eigsolve)" begin
    N = 30
    Δ = 0.5
    U = 0.1
    κ = 0.1
    F = 0.5

    a = destroy(N)
    H = Δ * a' * a + U / 2 * a' * a' * a * a + F * (a + a')
    H_gpu = cu(H)

    c_ops = [sqrt(κ) * a]

    L = liouvillian(H, c_ops)
    L_gpu = CuSparseMatrixCSR(L)

    # Dense eigen solver for Hamiltonian
    vals_H_cpu, vecs_H_cpu = eigenstates(H)
    vals_H_gpu, vecs_H_gpu = eigenstates(H_gpu)

    @test vals_H_cpu ≈ Array(vals_H_gpu) atol = 1.0e-8
    @test all(zip(vecs_H_cpu, vecs_H_gpu)) do (v_cpu, v_gpu)
        return isapprox(abs(dot(v_cpu.data, Array(v_gpu.data))), 1; atol = 1.0e-8)
    end

    # Sparse eigen solver for Liouvillian
    vals_cpu, vecs_cpu = eigenstates(L; sparse = Val(true), sigma = 0.01, eigvals = 4, krylovdim = 30)
    vals_gpu, vecs_gpu = eigenstates(
        L_gpu;
        sparse = Val(true),
        sigma = 0.01,
        eigvals = 4,
        krylovdim = 30,
        solver = LUFactorization(),
        v0 = normalize!(CUDA.rand(ComplexF64, size(L_gpu, 1))),
    )

    @test vals_cpu ≈ vals_gpu atol = 1.0e-8
    @test all(zip(vecs_cpu, vecs_gpu)) do (v_cpu, v_gpu)
        return isapprox(abs(dot(v_cpu.data, Array(v_gpu.data))), 1; atol = 1.0e-8)
    end

    # Arnoldi-Lindblad eigen solver for Liouvillian
    vals_al_cpu, vecs_al_cpu = eigsolve_al(L, 1 \ (30 * κ); eigvals = 4, liouvillian_eigs = Val(false))
    vals_al_gpu, vecs_al_gpu = eigsolve_al(
        L_gpu,
        1 \ (30 * κ);
        eigvals = 4,
        liouvillian_eigs = Val(false),
        ρ0 = cu(operator_to_vector(rand_dm(N))),
    )

    sort_func = x -> (round(abs(x), digits = 6), round(real(x), digits = 6), round(imag(x), digits = 6))

    idxs_al_cpu = sortperm(vals_al_cpu, by = sort_func, rev = false)
    idxs_al_gpu = sortperm(vals_al_gpu, by = sort_func, rev = false)
    vals_al_cpu = vals_al_cpu[idxs_al_cpu]
    vecs_al_cpu = vecs_al_cpu[idxs_al_cpu]
    vals_al_gpu = vals_al_gpu[idxs_al_gpu]
    vecs_al_gpu = vecs_al_gpu[idxs_al_gpu]

    @test vals_al_cpu ≈ vals_al_gpu atol = 1.0e-7
    @test all(zip(vecs_al_cpu, vecs_al_gpu)) do (v_cpu, v_gpu)
        return isapprox(abs(dot(v_cpu.data, Array(v_gpu.data))), 1; atol = 1.0e-8)
    end
end
