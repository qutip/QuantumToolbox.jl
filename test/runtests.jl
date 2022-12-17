using QuPhys
using Test


@testset "QuPhys.jl" begin
    # Write your tests here.
end

@testset "QuantumObjects" begin
    a = rand(ComplexF64, 10)
    @test_logs (:warn, "The norm of the input data is not one.") QuantumObject(a)
    a2 = QuantumObject(a, type=BraQuantumObject)
    a3 = QuantumObject(a, type=KetQuantumObject)
    @test isket(a2) == false
    @test isbra(a2) == true
    @test isoper(a2) == false
    @test issuper(a2) == false
    @test isket(a3) == true
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == false

    a = sprand(ComplexF64, 100, 100, 0.1)
    a2 = QuantumObject(a, type=OperatorQuantumObject)
    a3 = QuantumObject(a, type=SuperOperatorQuantumObject)

    @test isket(a2) == false
    @test isbra(a2) == false
    @test isoper(a2) == true
    @test issuper(a2) == false
    @test isket(a3) == false
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == true

    a = Array(a)
    a4 = QuantumObject(a)
    a5 = sparse(a4)
    @test isequal(a5, a2)
    @test (a5 == a3) == false
    @test a5 ≈ a2

    @test +a2 == a2
    @test -(-a2) == a2
    @test a2^3 ≈ a2 * a2 * a2
    @test a2 + 2 == 2 + a2
    @test (a2 + 2).data == a2.data + 2 * I
    @test a2 * 2 == 2 * a2

    @test transpose(transpose(a2)) == a2
    @test transpose(a2).data == transpose(a2.data)
    @test adjoint(adjoint(a2)) == a2
    @test adjoint(a2).data == adjoint(a2.data)

    N = 10
    a = fock(N, 3)
    @test ket2dm(a) ≈ projection(N, 3, 3)
    @test isket(a') == false
    @test isbra(a') == true
    @test size(a) == (N,)
    @test size(a') == (1, N)
    @test norm(a) ≈ 1
    @test norm(a') ≈ 1

    a = QuantumObject(rand(ComplexF64, N))
    @test (norm(a) ≈ 1) == false
    @test (norm(normalize(a)) ≈ 1) == true
    @test (norm(a) ≈ 1) == false # Again, to be sure that it is still non-normalized
    normalize!(a)
    @test (norm(a) ≈ 1) == true

    a = destroy(N)
    a_d = a'
    X = a + a_d
    Y = 1im * (a - a_d)
    Z = a + transpose(a)
    @test ishermitian(X) == true
    @test ishermitian(Y) == true
    @test issymmetric(Y) == false
    @test issymmetric(Z) == true

    @test Y[1, 2] == conj(Y[2, 1])

    @test triu(X) == a
    @test tril(X) == a_d

    triu!(X)
    @test X == a
    tril!(X)
    @test nnz(X) == 0

    @test eigvals(a_d * a) ≈ 0:9
end

@testset "Time Evolution and partial trace" begin
    N = 10

    a_d = kron(create(N), eye(2))
    a = a_d'
    sm = kron(eye(N), sigmam())
    sp = sm'
    sx = kron(eye(N), sigmax())
    sy = kron(eye(N), sigmay())
    sz = kron(eye(N), sigmaz())
    η = 0.01
    H = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx
    psi0 = kron(fock(N, 0), fock(2, 0))
    t_l = LinRange(0, 1000, 1000)
    e_ops = [a_d * a]
    sol = sesolve(H, psi0, t_l, e_ops=e_ops, progress=false)
    @test sum(abs.(sol.expect[1, :] .- sin.(η * t_l) .^ 2)) / length(t_l) < 0.1
    sol = sesolve(H, psi0, t_l, e_ops=e_ops, alg=Vern7(), progress=false, abstol=1e-7, reltol=1e-5)
    @test sum(abs.(sol.expect[1, :] .- sin.(η * t_l) .^ 2)) / length(t_l) < 0.1

    a = destroy(N)
    a_d = a'
    H = a_d * a
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = basis(N, 3)
    t_l = LinRange(0, 100, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops=e_ops, alg=Vern7(), progress=false)
    sol_mc = mcsolve(H, psi0, t_l, c_ops, n_traj=500, e_ops=e_ops, progress=false)
    @test sum(abs.(sol_mc.expect .- sol_me.expect)) / length(t_l) < 0.1

    sp1 = kron(sigmap(), eye(2))
    sm1 = sp1'
    sx1 = sm1 + sp1
    sy1 = 1im * (sm1 - sp1)
    sz1 = sp1 * sm1 - sm1 * sp1
    sp2 = kron(eye(2), sigmap())
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
    t_l = LinRange(0, 20 / γ1, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops=[sp1 * sm1, sp2 * sm2], progress=false)
    sol_mc = mcsolve(H, psi0, t_l, c_ops, n_traj=500, e_ops=[sp1 * sm1, sp2 * sm2], progress=false)
    @test sum(abs.(sol_mc.expect[1:2, :] .- sol_me.expect[1:2, :])) / length(t_l) < 0.1

    @test expect(sp1 * sm1, sol_me.states[300]) ≈ expect(sigmap() * sigmam(), ptrace(sol_me.states[300], [1]))


    ### DYNAMICAL FOCK DIMENSION ###
    F, Δ, κ = 5, 0.25, 1
    t_l = LinRange(0, 15, 100)

    N0 = 140
    a0 = destroy(N0)
    H0 = Δ * a0' * a0 + F * (a0 + a0')
    c_ops0 = [√κ * a0]
    e_ops0 = [a0' * a0]
    ψ00 = fock(N0, 0)
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops=e_ops0, alg=Vern7(), progress=false, saveat=[t_l[end]])

    function H_dfd(dims::AbstractVector)
        a = destroy(dims[1])
        Δ * a' * a + F * (a + a')
    end
    function c_ops_dfd(dims::AbstractVector)
        a = destroy(dims[1])
        [√κ * a]
    end
    function e_ops_dfd(dims::AbstractVector)
        a = destroy(dims[1])
        [a' * a]
    end
    maxdims = [150]
    ψ0 = fock(3, 0)
    sol = dfd_mesolve(H_dfd, ψ0, t_l, c_ops_dfd, e_ops_dfd, maxdims, progress=false,
        saveat=[t_l[end]], abstol=1e-9, reltol=1e-7)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) < 0.01

    F, Δ, κ, J = 1.5, 0.25, 1, 0.05
    N0 = 25
    N1 = 20
    a0 = kron(destroy(N0), eye(N1))
    a1 = kron(eye(N0), destroy(N1))
    H0 = Δ * a0' * a0 + F * (a0 + a0') + Δ * a1' * a1 + J * (a0' * a1 + a0 * a1')
    c_ops0 = [√κ * a0, √κ * a1]
    e_ops0 = [a0' * a0, a1' * a1]
    ψ00 = kron(fock(N0, 0), fock(N1, 15))
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops=e_ops0, alg=Vern7(), progress=false, saveat=[t_l[end]])

    function H_dfd2(dims::AbstractVector)
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        Δ * a' * a + F * (a + a') + Δ * b' * b + J * (a' * b + a * b')
    end
    function c_ops_dfd2(dims::AbstractVector)
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        [√κ * a, √κ * b]
    end
    function e_ops_dfd2(dims::AbstractVector)
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        [a' * a, b' * b]
    end
    maxdims = [50, 50]
    ψ0 = kron(fock(3, 0), fock(20, 15))
    sol = dfd_mesolve(H_dfd2, ψ0, t_l, c_ops_dfd2, e_ops_dfd2, maxdims, progress=false, saveat=[t_l[end]])

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) +
          sum(abs.((sol.expect[2, :] .- sol0.expect[2, :]) ./ (sol0.expect[2, :] .+ 1e-16))) < 0.01
end

@testset "Eigenvalues and Operators" begin
    N = 30
    a = kron(destroy(N), eye(2))
    a_d = a'

    sm = kron(eye(N), sigmam())
    sp = sm'
    sx = kron(eye(N), sigmax())
    sy = kron(eye(N), sigmay())
    sz = kron(eye(N), sigmaz())

    η = 0.2
    H_d = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx + η^2
    H_c = a_d * a + 0.5 * (sz * cosm(2 * η * (a + a_d)) + sy * sinm(2 * η * (a + a_d)))

    vals_d, vecs_d = eigen(H_d)
    vals_c, vecs_c = eigen(Hermitian(H_c))

    @test sum(vals_d[1:20] .- vals_c[1:20]) / 20 < 1e-3
end

@testset "Steadystate" begin
    N = 10

    a = destroy(N)
    a_d = a'
    H = a_d * a + 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops=e_ops, progress=false)
    ρ_ss = steadystate(H, c_ops)
    @test abs(sol_me.expect[1, end] - expect(e_ops[1], ρ_ss)) < 1e-3

    H = a_d * a
    H_t = 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops=e_ops, H_t=(t) -> sin(t) * H_t, alg=Vern7(), progress=false)
    ρ_ss = steadystate_floquet(H, c_ops, -1im * 0.5 * H_t, 1im * 0.5 * H_t, 1)
    @test abs(sum(sol_me.expect[1, end-100:end]) / 101 - expect(e_ops[1], ρ_ss)) < 1e-2
end

@testset "Entanglement" begin
    g = fock(2, 1)
    e = fock(2, 0)
    state = normalize(kron(g, e) + kron(e, g))
    rho = state * state'
    @test entanglement(state, [1]) / log(2) ≈ 1
end

@testset "Wigner" begin
    α = 0.5 + 0.8im
    ψ = coherent(30, α)
    xvec = LinRange(-3, 3, 300)
    yvec = LinRange(-3, 3, 300)

    wig = wigner(ψ, xvec, yvec)

    X, Y = meshgrid(xvec, yvec)
    wig_tmp1 = gaussian(xvec / √2, real(α), 1 / 2)
    wig_tmp2 = gaussian(yvec / √2, imag(α), 1 / 2)
    wig2 = maximum(wig) * reshape(kron(wig_tmp1, wig_tmp2), 300, 300)

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 0.1
end