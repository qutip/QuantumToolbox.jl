@testset "Quantum Objects" begin
    # unsupported size of array
    for a in [rand(ComplexF64, 3, 2), rand(ComplexF64, 2, 2, 2)]
        for t in [Nothing, KetQuantumObject, BraQuantumObject, OperatorQuantumObject, SuperOperatorQuantumObject]
            @test_throws DomainError Qobj(a, type=t)
        end
    end

    N = 10
    a = rand(ComplexF64, 10)
    # @test_logs (:warn, "The norm of the input data is not one.") QuantumObject(a)
    @test_throws DomainError Qobj(a,  type=BraQuantumObject)
    @test_throws DomainError Qobj(a,  type=OperatorQuantumObject)
    @test_throws DomainError Qobj(a,  type=SuperOperatorQuantumObject)
    @test_throws DomainError Qobj(a', type=KetQuantumObject)
    @test_throws DomainError Qobj(a', type=OperatorQuantumObject)
    @test_throws DomainError Qobj(a', type=SuperOperatorQuantumObject)
    @test_throws DimensionMismatch Qobj(a, dims=[2])
    @test_throws DimensionMismatch Qobj(a', dims=[2])
    a2 = Qobj(a')
    a3 = Qobj(a)
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
    @test isket(Qobj(Matrix([2 3])')) == true

    a = sprand(ComplexF64, 100, 100, 0.1)
    a2 = Qobj(a)
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