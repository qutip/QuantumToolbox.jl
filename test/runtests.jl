using QuPhys
using Test


@testset "QuantumObjects" begin
    a = rand(ComplexF64, 10)
    # @test_logs (:warn, "The norm of the input data is not one.") QuantumObject(a)
    a2 = Qobj(a, type=BraQuantumObject)
    a3 = Qobj(a, type=KetQuantumObject)
    @test isket(a2) == false
    @test isbra(a2) == true
    @test isoper(a2) == false
    @test issuper(a2) == false
    @test isket(a3) == true
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == false
    @test Qobj(a3) == a3
    @test !(Qobj(a3) === a3)

    a = sprand(ComplexF64, 100, 100, 0.1)
    a2 = Qobj(a, type=OperatorQuantumObject)
    a3 = Qobj(a, type=SuperOperatorQuantumObject)

    @test isket(a2) == false
    @test isbra(a2) == false
    @test isoper(a2) == true
    @test issuper(a2) == false
    @test isket(a3) == false
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == true

    a = Array(a)
    a4 = Qobj(a)
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
    @test sparse(ket2dm(a)) ≈ projection(N, 3, 3)
    @test isket(a') == false
    @test isbra(a') == true
    @test size(a) == (N,)
    @test size(a') == (1, N)
    @test norm(a) ≈ 1
    @test norm(a') ≈ 1

    ψ = Qobj(normalize(rand(ComplexF64, N)))
    @test dot(ψ, ψ) ≈ norm(ψ)
    @test dot(ψ, ψ) ≈ ψ' * ψ

    a = Qobj(rand(ComplexF64, N))
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

    # Eigenvalues
    @test eigvals(a_d * a) ≈ 0:9

    # Random density matrix
    ρ = rand_dm(10)
    @test tr(ρ) ≈ 1
    @test isposdef(ρ) == true

    # Expectation value
    a = destroy(10)
    ψ = normalize(fock(10, 3) + 1im * fock(10, 4))
    @test expect(a, ψ) ≈ expect(a, ψ')
    ψ = fock(10, 3)
    @test norm(ψ' * a) ≈ 2
    @test expect(a' * a, ψ' * a) ≈ 16

    # REPL show
    a = destroy(N)
    ψ = fock(N, 3)

    opstring = sprint((t, s) -> show(t, "text/plain", s), a)
    datastring = sprint((t, s) -> show(t, "text/plain", s), a.data)
    a_dims = a.dims
    a_size = size(a)
    a_isherm = ishermitian(a)
    @test opstring == "Quantum Object:   type=Operator   dims=$a_dims   size=$a_size   ishermitian=$a_isherm\n$datastring"

    a = spre(a)
    opstring = sprint((t, s) -> show(t, "text/plain", s), a)
    datastring = sprint((t, s) -> show(t, "text/plain", s), a.data)
    a_dims = a.dims
    a_size = size(a)
    a_isherm = ishermitian(a)
    @test opstring == "Quantum Object:   type=SuperOperator   dims=$a_dims   size=$a_size\n$datastring"

    opstring = sprint((t, s) -> show(t, "text/plain", s), ψ)
    datastring = sprint((t, s) -> show(t, "text/plain", s), ψ.data)
    ψ_dims = ψ.dims
    ψ_size = size(ψ)
    @test opstring == "Quantum Object:   type=Ket   dims=$ψ_dims   size=$ψ_size\n$datastring"

    ψ = ψ'
    opstring = sprint((t, s) -> show(t, "text/plain", s), ψ)
    datastring = sprint((t, s) -> show(t, "text/plain", s), ψ.data)
    ψ_dims = ψ.dims
    ψ_size = size(ψ)
    @test opstring == "Quantum Object:   type=Bra   dims=$ψ_dims   size=$ψ_size\n$datastring"

    ψ = coherent(30, 3)
    α, δψ = get_coherence(ψ)
    @test isapprox(abs(α), 3, atol=1e-5) && abs2(δψ[1]) > 0.999

    # Broadcasting
    a = destroy(20)
    for op in ((+), (-), (*), (^))
        A = broadcast(op, a, a)
        @test A.data == broadcast(op, a.data, a.data) && A.type == a.type && A.dims == a.dims

        A = broadcast(op, 2.1, a)
        @test A.data == broadcast(op, 2.1, a.data) && A.type == a.type && A.dims == a.dims

        A = broadcast(op, a, 2.1)
        @test A.data == broadcast(op, a.data, 2.1) && A.type == a.type && A.dims == a.dims
    end

    # tidyup tests
    ρ1 = rand_dm(20)
    ρ2 = dense_to_sparse(ρ1)
    @test tidyup!(ρ2, 0.1) == ρ2 != ρ1
    @test dense_to_sparse(tidyup!(ρ1, 0.1)) == ρ2

    ρ1 = rand_dm(20)
    ρ2 = dense_to_sparse(ρ1)
    @test tidyup(ρ2, 0.1) != ρ2
    @test dense_to_sparse(tidyup(ρ1, 0.1)) == tidyup(ρ2, 0.1)
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
    # sol = sesolve(H, psi0, t_l, e_ops=e_ops, alg=LinearExponential(krylov=:adaptive, m=15), progress=false)
    # @test sum(abs.(sol.expect[1, :] .- sin.(η * t_l) .^ 2)) / length(t_l) < 0.1
    sol = sesolve(H, psi0, t_l, e_ops=e_ops, alg=Vern7(), progress=false)
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
    @test expect(sp1 * sm1, sol_me.states[end]) ≈ expect(sigmap() * sigmam(), ptrace(sol_me.states[end], [1]))
end

@testset "Dynamical Fock Dimension mesolve" begin
    ### DYNAMICAL FOCK DIMENSION ###
    F, Δ, κ = 5, 0.25, 1
    t_l = range(0, 15, length=100)

    N0 = 140
    a0 = destroy(N0)
    H0 = Δ * a0' * a0 + F * (a0 + a0')
    c_ops0 = [√κ * a0]
    e_ops0 = [a0' * a0]
    ψ00 = fock(N0, 0)
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops=e_ops0, progress=false)

    function H_dfd0(dims, p)
        Δ = p.Δ
        F = p.F
        a = destroy(dims[1])
        Δ * a' * a + F * (a + a')
    end
    function c_ops_dfd0(dims, p)
        κ = p.κ
        a = destroy(dims[1])
        [√κ * a]
    end
    function e_ops_dfd0(dims, p)
        a = destroy(dims[1])
        [a' * a]
    end
    maxdims = [150]
    ψ0 = fock(3, 0)
    dfd_params = (Δ=Δ, F=F, κ=κ)
    sol = dfd_mesolve(H_dfd0, ψ0, t_l, c_ops_dfd0, maxdims, dfd_params, e_ops=e_ops_dfd0, progress=false);

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) < 0.01

    ######################

    F = 0
    H0 = Δ * a0' * a0 + F * (a0 + a0')
    c_ops0 = [√κ * a0]
    e_ops0 = [a0' * a0]
    ψ00 = fock(N0, 50)
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops=e_ops0, progress=false)

    function H_dfd1(dims, p)
        Δ = p.Δ
        F = p.F
        a = destroy(dims[1])
        Δ * a' * a + F * (a + a')
    end
    function c_ops_dfd1(dims, p)
        κ = p.κ
        a = destroy(dims[1])
        [√κ * a]
    end
    function e_ops_dfd1(dims, p)
        a = destroy(dims[1])
        [a' * a]
    end
    maxdims = [150]
    ψ0 = fock(70, 50)
    dfd_params = (Δ=Δ, F=F, κ=κ)
    sol = dfd_mesolve(H_dfd1, ψ0, t_l, c_ops_dfd1, maxdims, dfd_params, e_ops=e_ops_dfd1, progress=false)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) < 0.01


    ######################

    F, Δ, κ, J = 1.5, 0.25, 1, 0.05
    N0 = 25
    N1 = 20
    a0 = kron(destroy(N0), eye(N1))
    a1 = kron(eye(N0), destroy(N1))
    H0 = Δ * a0' * a0 + F * (a0 + a0') + Δ * a1' * a1 + J * (a0' * a1 + a0 * a1')
    c_ops0 = [√κ * a0, √κ * a1]
    e_ops0 = [a0' * a0, a1' * a1]
    ψ00 = kron(fock(N0, 0), fock(N1, 15))
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops=e_ops0, progress=false)

    function H_dfd2(dims, p)
        Δ = p.Δ
        F = p.F
        J = p.J
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        Δ * a' * a + F * (a + a') + Δ * b' * b + J * (a' * b + a * b')
    end
    function c_ops_dfd2(dims, p)
        κ = p.κ
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        [√κ * a, √κ * b]
    end
    function e_ops_dfd2(dims, p)
        a = kron(destroy(dims[1]), eye(dims[2]))
        b = kron(eye(dims[1]), destroy(dims[2]))
        [a' * a, b' * b]
    end
    maxdims = [50, 50]
    ψ0 = kron(fock(3, 0), fock(20, 15))
    dfd_params = (Δ=Δ, F=F, κ=κ, J=J)
    sol = dfd_mesolve(H_dfd2, ψ0, t_l, c_ops_dfd2, maxdims, dfd_params, e_ops=e_ops_dfd2, progress=false)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) +
          sum(abs.((sol.expect[2, :] .- sol0.expect[2, :]) ./ (sol0.expect[2, :] .+ 1e-16))) < 0.01
end

@testset "Dynamical Shifted Fock" begin
    F   = 3
    Δ   = 0.25
    κ   = 1
    U = 0.01

    tlist = LinRange(0,25,300)

    # Single cavity case
    N0     = 100
    a0     = destroy(N0)
    H0     = Δ*a0'*a0 + F*(a0+a0') + U * a0'^2 * a0^2
    c_ops0 = [√(κ)*a0]
    
    α0 = 1.5
    ρ0   = coherent(N0, α0)
    sol0 = mesolve(H0, ρ0, tlist, c_ops0, e_ops=[a0'*a0, a0], progress=false)

    N = 5
    a        = destroy(N)
    function H_dsf(op_list, p)
        Δ = p.Δ
        F = p.F
        U = p.U
        a = op_list[1]
        Δ*a'*a + F*(a + a') + U * a'^2 * a^2
    end
    function c_ops_dsf(op_list, p)
        κ = p.κ
        a = op_list[1]
        [√κ * a]
    end
    function e_ops_dsf(op_list, p)
        a = op_list[1]
        [a' * a, a]
    end
    op_list = [a]
    ψ0  = fock(N, 0)
    α0_l = [α0]
    dsf_params = (Δ=Δ, F=F, κ=κ, U=U)
    
    sol_dsf_me = dsf_mesolve(H_dsf, ψ0, tlist, c_ops_dsf, op_list, α0_l, dsf_params, e_ops=e_ops_dsf, progress=false)
    sol_dsf_mc = dsf_mcsolve(H_dsf, ψ0, tlist, c_ops_dsf, op_list, α0_l, dsf_params, e_ops=e_ops_dsf, progress=false, n_traj=500)
    val_ss = abs2(sol0.expect[1,end])
    @test sum(abs2.(sol0.expect[1,:] .- sol_dsf_me.expect[1,:])) / (val_ss * length(tlist)) < 0.1
    @test sum(abs2.(sol0.expect[1,:] .- sol_dsf_mc.expect[1,:])) / (val_ss * length(tlist)) < 0.1

    # Two cavities case
    F   = 2
    Δ   = 0.25
    κ   = 1
    U = 0.01
    J = 0.5
    tlist = LinRange(0,15,300)

    N0     = 20
    a10     = kron(destroy(N0), eye(N0))
    a20     = kron(eye(N0), destroy(N0))
    H0     = Δ*a10'*a10 + Δ*a20'*a20 + U*a10'^2*a10^2 + U*a20'^2*a20^2 + F*(a10+a10') + J*(a10'*a20 + a10*a20')
    c_ops0 = [√κ*a10, √κ*a20]

    ρ0   = kron(coherent(N0, α0), coherent(N0, α0))
    sol0 = mesolve(H0, ρ0, tlist, c_ops0, e_ops=[a10'*a10, a20'*a20], progress=false)

    N = 5
    a1 = kron(destroy(N), eye(N))
    a2 = kron(eye(N), destroy(N))
    function H_dsf2(op_list, p)
        Δ = p.Δ
        F = p.F
        U = p.U
        J = p.J
        a1, a2 = op_list
        Δ*a1'*a1 + Δ*a2'*a2 + U*a1'^2*a1^2 + U*a2'^2*a2^2 + F*(a1 + a1') + J*(a1'*a2 + a1*a2')
    end
    function c_ops_dsf2(op_list, p)
        κ = p.κ
        a1, a2 = op_list
        [√κ * a1, √κ * a2]
    end
    function e_ops_dsf2(op_list, p)
        a1, a2 = op_list
        [a1' * a1, a2' * a2]
    end
    op_list = [a1, a2]
    ψ0  = kron(fock(N, 0), fock(N, 0))
    α0_l = [α0, α0]
    dsf_params = (Δ=Δ, F=F, κ=κ, U=U, J=J)

    sol_dsf_me = dsf_mesolve(H_dsf2, ψ0, tlist, c_ops_dsf2, op_list, α0_l, dsf_params, e_ops=e_ops_dsf2, progress=false)
    sol_dsf_mc = dsf_mcsolve(H_dsf2, ψ0, tlist, c_ops_dsf2, op_list, α0_l, dsf_params, e_ops=e_ops_dsf2, progress=false, n_traj=500)

    val_ss = abs2(sol0.expect[1,end])
    @test sum(abs2.(sol0.expect[1,:] .- sol_dsf_me.expect[1,:])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[1,:] .- sol_dsf_mc.expect[1,:])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[2,:] .- sol_dsf_me.expect[2,:])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[2,:] .- sol_dsf_mc.expect[2,:])) / (val_ss * length(tlist)) < 0.6
end

@testset "Generalized Master Equation" begin
    N_c = 30
    N_trunc = 10
    tol=1e-14
    
    a = kron(destroy(N_c), eye(2))
    sm = kron(eye(N_c), sigmam())
    sp = sm'
    sx = sm + sp
    sz = sp * sm - sm * sp
    
    H = 1 * a' * a + 1 * sz / 2 + 0.5 * (a + a') * sx
    
    fields = [sqrt(0.01) * (a + a'), sqrt(0.01) * sx]
    Tlist = [0, 0.0]
    
    E, U, L1 = liouvillian_generalized(H, fields, Tlist, N_trunc=N_trunc, tol=tol)
    Ω = dense_to_sparse((E' .- E)[1:N_trunc,1:N_trunc], tol)
    
    H_d = Qobj(dense_to_sparse((U' * H * U)[1:N_trunc,1:N_trunc], tol))
    Xp = Qobj( Ω .* dense_to_sparse(triu((U' * (a + a') * U).data[1:N_trunc,1:N_trunc], 1), tol))
    a2 = Qobj( dense_to_sparse((U' * a * U).data[1:N_trunc,1:N_trunc], tol))
    sm2 = Qobj( dense_to_sparse((U' * sm * U).data[1:N_trunc,1:N_trunc], tol))
    
    # Standard liouvillian case
    c_ops = [sqrt(0.01) * a2, sqrt(0.01) * sm2] 
    L2 = liouvillian(H_d, c_ops)

    @test (expect(Xp'*Xp, steadystate(L1)) < 1e-10 && expect(Xp'*Xp, steadystate(L2)) > 1e-3)

    H = 1 * a' * a + 1 * sz / 2 + 1e-5 * (a * sp + a' * sm)

    Tlist = [0.2, 0.0]

    E, U, L1 = liouvillian_generalized(H, fields, Tlist, N_trunc=N_trunc, tol=tol)
    Ω = dense_to_sparse((E' .- E)[1:N_trunc,1:N_trunc], tol)

    H_d = Qobj(dense_to_sparse((U' * H * U)[1:N_trunc,1:N_trunc], tol))
    Xp = Qobj( Ω .* dense_to_sparse(triu((U' * (a + a') * U).data[1:N_trunc,1:N_trunc], 1), tol))
    a2 = Qobj( dense_to_sparse((U' * a * U).data[1:N_trunc,1:N_trunc], tol))
    sm2 = Qobj( dense_to_sparse((U' * sm * U).data[1:N_trunc,1:N_trunc], tol))

    @test abs(expect(Xp'*Xp, steadystate(L1)) - n_th(1, Tlist[1])) / n_th(1, Tlist[1]) < 1e-4
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
    vals_c, vecs_c = eigen((H_c))
    vals2, vecs2 = eigsolve(H_d, sigma=-0.9, k=10, krylovdim=30)
    sort!(vals_c, by=real)
    sort!(vals2, by=real)

    @test sum(real.(vals_d[1:20]) .- real.(vals_c[1:20])) / 20 < 1e-3
    @test sum(real.(vals_d[1:10]) .- real.(vals2[1:10])) / 20 < 1e-3


    N = 5
    a = kron(destroy(N), eye(N))
    a_d = a'
    b = kron(eye(N), destroy(N))
    b_d = b'

    ωc = 1
    ωb = 1
    g = 0.01
    κ = 0.1
    n_th = 0.01

    H = ωc * a_d * a + ωb * b_d * b + g * (a + a_d) * (b + b_d)
    c_ops = [√((1+n_th)*κ) * a, √κ * b, √(n_th*κ) * a_d]
    L = liouvillian(H, c_ops).data

    vals, vecs = eigsolve(L, sigma=0.01, k=10, krylovdim=50)
    vals2, vecs2 = eigen(sparse_to_dense(L))
    vals3, vecs3 = eigsolve_al(liouvillian(H, c_ops), 1\(40*κ), k=10, krylovdim=50)
    idxs = sortperm(vals2, by=abs)
    vals2 = vals2[idxs][1:10]
    vecs2 = vecs2[:, idxs][:, 1:10]

    @test isapprox(sum(abs2, vals), sum(abs2, vals2), atol=1e-7)
    @test isapprox(abs2(vals2[1]), abs2(vals3[1]), atol=1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im*angle(vecs[1,1])), vec2mat(vecs2[:, 1]), atol=1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im*angle(vecs[1,1])), vec2mat(vecs3[:, 1]), atol=1e-5)
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
    function H_t_f(t,p)
        copyto!(p.L_t_cache2.data, p.L_t.data)
        lmul!(sin(t), p.L_t_cache2)
    end
    L_t = liouvillian(H_t)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops=e_ops, H_t=H_t_f, alg=Vern7(), progress=false, params=(L_t=L_t, L_t_cache2=deepcopy(L_t)))
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
    ρ = dense_to_sparse(ket2dm(ψ), 1e-6)
    xvec = LinRange(-3, 3, 300)
    yvec = LinRange(-3, 3, 300)

    wig = wigner(ψ, xvec, yvec, solver=WignerLaguerre(tol=1e-6))
    wig2 = wigner(ρ, xvec, yvec, solver=WignerLaguerre(parallel=false))
    wig3 = wigner(ρ, xvec, yvec, solver=WignerLaguerre(parallel=true))
    wig4 = wigner(ψ, xvec, yvec, solver=WignerClenshaw())

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 1e-3
    @test sqrt(sum(abs.(wig3 .- wig)) / length(wig)) < 1e-3
    @test sqrt(sum(abs.(wig4 .- wig)) / length(wig)) < 1e-3

    X, Y = meshgrid(xvec, yvec)
    wig_tmp1 = gaussian.(xvec / √2, real(α), 1 / 2)
    wig_tmp2 = gaussian.(yvec / √2, imag(α), 1 / 2)
    wig2 = maximum(wig) * reshape(kron(wig_tmp1, wig_tmp2), 300, 300)

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 0.1
end

@testset "Permutation" begin
    # Block Diagonal Form
    N = 20
    Δ = 0
    G = 5
    tg = 0
    θ  = atan(tg)
    U  = sin(θ)
    κ2 = cos(θ)
    κ1  = 0.
    κϕ  = 1e-3
    nth = 0.

    a     = destroy(N)
    ad    = create(N)
    H     = -Δ*ad*a + G/2*(ad^2 + a^2) + U/2*(ad^2*a^2)
    c_ops = [√(κ2)*a^2, √(κ1*(nth+1))*a, √(κ1*nth)*ad, √(κϕ)*ad*a]
    L     = liouvillian(H,c_ops)

    P, L_bd, block_sizes = bdf(L)
    blocks_list, block_indices = get_bdf_blocks(L_bd, block_sizes)
    @test size(L_bd) == size(L)
    @test length(block_sizes) == 4
    @test length(blocks_list) == 4
    @test length(block_indices) == 4
    @test sum(block_sizes .== 100) == 4
end

@testset "Correlations and Spectrum" begin
    a = destroy(10)
    H = a' * a
    c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']

    ω_l = range(0, 3, length=1000)
    ω_l1, spec1 = spectrum(H, ω_l, a', a, c_ops, solver=FFTCorrelation(), progress=false)
    ω_l2, spec2 = spectrum(H, ω_l, a', a, c_ops)
    spec1 = spec1 ./ maximum(spec1)
    spec2 = spec2 ./ maximum(spec2)

    test_func1 = maximum(real.(spec1)) * (0.1/2)^2 ./ ((ω_l1 .- 1).^2 .+ (0.1/2)^2)
    test_func2 = maximum(real.(spec2)) * (0.1/2)^2 ./ ((ω_l2 .- 1).^2 .+ (0.1/2)^2)
    idxs1 = test_func1 .> 0.05
    idxs2 = test_func2 .> 0.05
    @test sum(abs2.(spec1[idxs1] .- test_func1[idxs1])) / sum(abs2.(test_func1[idxs1])) < 0.01
    @test sum(abs2.(spec2[idxs2] .- test_func2[idxs2])) / sum(abs2.(test_func2[idxs2])) < 0.01
end


@testset "LowRankDynamics" begin
    # Define lattice
    Nx,Ny = 2, 3
    latt  = Lattice(Nx=Nx, Ny=Ny)
    N_cut = 2         
    N_modes = latt.N  
    N = N_cut^N_modes 
    M = Nx*Ny+1       

    # Define initial state
    ϕ = Vector{QuantumObject{Vector{ComplexF64}, KetQuantumObject}}(undef, M)
    ϕ[1] = kron(repeat([basis(2,0)],N_modes)...)
    global i=1
    for j in 1:N_modes
        global i+=1
        i<=M && (ϕ[i] = mb(sp, j, latt) * ϕ[1])
    end
    for k in 1:N_modes-1
        for l=k+1:N_modes
            global i+=1
            i<=M && (ϕ[i] = mb(sp, k, latt) * mb(sp, l, latt) * ϕ[1])
        end
    end
    for i in i+1:M
        ϕ[i] = Qobj(rand(ComplexF64,size(ϕ[1])[1]), dims=ϕ[1].dims)
        normalize!(ϕ[i])
    end
    z  = hcat(broadcast(x->x.data, ϕ)...)
    B  = Matrix(Diagonal([1+0im; zeros(M-1)]))
    S  = z'*z 
    B = B / tr(S*B)
    ρ = Qobj(z*B*z', dims=ones(Int,N_modes)*N_cut)

    # Define Hamiltonian and collapse operators
    Jx = 0.9 
    Jy = 1.02
    Jz = 1. 
    hx = 0.
    γ  = 1
    Sz = sum([mb(sz, i, latt) for i in 1:latt.N])
    tl = LinRange(0,10,100)

    H, c_ops = TFIM(Jx, Jy, Jz, hx, γ, latt; bc=pbc, order=1)
    e_ops = (Sz,)

    # Full solution
    mesol = mesolve(H, ρ, tl, c_ops; e_ops=[e_ops...]);
    A = Matrix(mesol.states[end].data)
    λ = eigvals(Hermitian(A))
    Strue  = -sum(λ.*log2.(λ))

    # Low rank solution
    function f_entropy(p,z,B)
        C = p.A0
        σ = p.Bi
        mul!(C, z, sqrt(B))
        mul!(σ, C', C)
        λ = eigvals(Hermitian(σ))
        λ = λ[λ.>1e-10]
        return -sum(λ .* log2.(λ))
    end

    opt = LRMesolveOptions( 
        alg        = Tsit5(),
        err_max    = 1e-3,
        p0         = 0.,
        atol_inv   = 1e-6,
        adj_condition="variational",
        Δt = 0.2, )
    lrsol = lr_mesolve(H, z, B, tl, c_ops; e_ops=e_ops, f_ops=(f_entropy,), opt=opt)

    # Test
    m_me = real(mesol.expect[1,:])
    m_lr = real(lrsol.expvals[1,:])
    @test all(abs.((m_me .- m_lr)./m_me).<0.1)

    S_lr = real(lrsol.funvals[1,end])
    @test abs((S_lr - Strue)/Strue) < 0.5
end