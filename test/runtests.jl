using QuPhys
using Test

using LinearAlgebra
using SparseArrays

@testset "QuPhys.jl" begin
    # Write your tests here.
end

@testset "Eigenvalues" begin
    A = sprand(ComplexF64, 15, 15, 0.1)
    A *= A'
    eigs1 = eigen(collect(A)).values
    eigs2 = eigensystem(A, k = 8, sigma = -5, krylovdim = 10)[1]
    eigs3 = eigensystem(A, k = 8, krylovdim = 10)[1]
    @test sum(abs.(eigs1[1:8] .- eigs2)) < 0.1
    @test sum(reverse(sort(abs.(eigs1))[end-7:end]) .- abs.(eigs3)) < 0.1
end

@testset "Time Evolution" begin
    N = 10

    a = kron(destroy(N), eye(2))
    a_d = a'
    sp = kron(eye(N), destroy(2))
    sm = sp'
    sx = sm + sp
    sy = 1im * (sm - sp)
    sz = sp * sm - sm * sp
    η = 0.01
    H = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx
    psi0 = kron(fock(N, 0), fock(2, 0))
    t_l = LinRange(0, 1000, 1000)
    e_ops = [a_d * a]
    sol, expect_se = sesolve(H, psi0, t_l, e_ops = e_ops)
    @test sum(abs.(expect_se[1:end, 1] .- sin.(η * t_l).^2)) / length(t_l) < 0.1

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

    sp1 = kron(destroy(2), eye(2))
    sm1 = sp1'
    sx1 = sm1 + sp1
    sy1 = 1im * (sm1 - sp1)
    sz1 = sp1 * sm1 - sm1 * sp1
    sp2 = kron(eye(2), destroy(2))
    sm2 = sp2'
    sx2 = sm2 + sp2
    sy2 = 1im * (sm2 - sp2)
    sz2 = sp2 * sm2 - sm2 * sp2
    ωq1, ωq2 = 1, 1
    γ1, γ2 = 0.05, 0.1
    H = 0.5 * ωq1 * sz1 + 0.5 * ωq2 * sz2
    c_ops = [sqrt(γ1) * sm1, sqrt(γ2) * sm2]
    psi0_1 = normalize(fock(2, 0) + fock(2, 1))
    psi0_2 = normalize(fock(2, 0) + fock(2, 1))
    psi0 = kron(psi0_1, psi0_2)
    t_l = LinRange(0, 10 / (γ1 + γ2), 1000)
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = [sp1 * sm1, sp2 * sm2]);
    sol_mc, expect_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = [sp1 * sm1, sp2 * sm2], ensemble_method = EnsembleSerial());
    @test sum(abs.(expect_mc[1:end, 1:2] .- expect_me[1:end, 1:2])) / length(t_l) < 0.1
end

@testset "Entanglement" begin
    g = fock(2, 1)
    e = fock(2, 0)
    state = normalize(kron(g, e) + kron(e, g))
    rho = state * state'
    @test entanglement(state, [1], (2, 2)) / log(2) ≈ 1
end