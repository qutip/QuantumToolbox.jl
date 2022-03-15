using QuPhys
using Test

using LinearAlgebra
using SparseArrays

@testset "QuPhys.jl" begin
    # Write your tests here.
end

@testset "Eigenvalues" begin
    A = sprand(ComplexF64, 15, 15, 0.1)
    A += A'
    eigs1 = eigen(collect(A)).values
    eigs2 = eigensystem(A, k = 8, sigma = -4, krylovdim = 10)[1]
    @test sum(abs.(eigs1[1:8] .- eigs2)) < 0.1
end

@testset "Time Evolution" begin
    N = 10
    a = destroy(N)
    a_d = a'
    H = a_d * a
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]

    psi0 = fock(N, 3)
    t_l = LinRange(0, 100, 1000)

    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops);
    sol_mc, expect_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = e_ops, ensemble_method = EnsembleSerial());

    @test sum(abs.(expect_mc .- expect_me)) / length(t_l) < 0.1
end