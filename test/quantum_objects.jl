@testset "Quantum Objects" begin
    # unsupported size of array
    for a in [rand(ComplexF64, 3, 2), rand(ComplexF64, 2, 2, 2)]
        for t in [nothing, Ket, Bra, Operator, SuperOperator, OperatorBra, OperatorKet]
            @test_throws DomainError Qobj(a, type = t)
        end
    end

    N = 10
    a = rand(ComplexF64, 10)
    # @test_logs (:warn, "The norm of the input data is not one.") QuantumObject(a)
    @test_throws DomainError Qobj(a, type = Bra)
    @test_throws DomainError Qobj(a, type = Operator)
    @test_throws DomainError Qobj(a, type = SuperOperator)
    @test_throws DomainError Qobj(a, type = OperatorBra)
    @test_throws DomainError Qobj(a', type = Ket)
    @test_throws DomainError Qobj(a', type = Operator)
    @test_throws DomainError Qobj(a', type = SuperOperator)
    @test_throws DomainError Qobj(a', type = OperatorKet)
    @test_throws DimensionMismatch Qobj(a, dims = [2])
    @test_throws DimensionMismatch Qobj(a', dims = [2])
    a2 = Qobj(a')
    a3 = Qobj(a)
    @test isket(a2) == false
    @test isbra(a2) == true
    @test isoper(a2) == false
    @test issuper(a2) == false
    @test isoperket(a2) == false
    @test isoperbra(a2) == false
    @test isket(a3) == true
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == false
    @test isoperket(a3) == false
    @test isoperbra(a3) == false
    @test Qobj(a3) == a3
    @test !(Qobj(a3) === a3)
    @test isket(Qobj(Matrix([2 3])')) == true

    a = sprand(ComplexF64, 100, 100, 0.1)
    a2 = Qobj(a)
    a3 = Qobj(a, type = SuperOperator)

    @test isket(a2) == false
    @test isbra(a2) == false
    @test isoper(a2) == true
    @test issuper(a2) == false
    @test isoperket(a2) == false
    @test isoperbra(a2) == false
    @test isket(a3) == false
    @test isbra(a3) == false
    @test isoper(a3) == false
    @test issuper(a3) == true
    @test isoperket(a3) == false
    @test isoperbra(a3) == false
    @test_throws DimensionMismatch Qobj(a, dims = [2])
    @test_throws DimensionMismatch Qobj(a, dims = [2])

    # Operator-Ket, Operator-Bra tests
    H = 0.3 * sigmax() + 0.7 * sigmaz()
    L = liouvillian(H)
    ρ = Qobj(rand(ComplexF64, 2, 2))
    ρ_ket = mat2vec(ρ)
    ρ_bra = ρ_ket'
    @test ρ_bra == Qobj(mat2vec(ρ.data)', type = OperatorBra)
    @test ρ == vec2mat(ρ_ket)
    @test isket(ρ_ket) == false
    @test isbra(ρ_ket) == false
    @test isoper(ρ_ket) == false
    @test issuper(ρ_ket) == false
    @test isoperket(ρ_ket) == true
    @test isoperbra(ρ_ket) == false
    @test isket(ρ_bra) == false
    @test isbra(ρ_bra) == false
    @test isoper(ρ_bra) == false
    @test issuper(ρ_bra) == false
    @test isoperket(ρ_bra) == false
    @test isoperbra(ρ_bra) == true
    @test ρ_bra.dims == [2]
    @test ρ_ket.dims == [2]
    @test H * ρ ≈ spre(H) * ρ
    @test ρ * H ≈ spost(H) * ρ
    @test H * ρ * H ≈ sprepost(H, H) * ρ
    @test (L * ρ_ket).dims == [2]
    @test L * ρ_ket ≈ -1im * (+(spre(H) * ρ_ket) - spost(H) * ρ_ket)
    @test (ρ_bra * L')' == L * ρ_ket
    @test sum((conj(ρ) .* ρ).data) ≈ dot(ρ_ket, ρ_ket) ≈ ρ_bra * ρ_ket
    @test_throws DimensionMismatch Qobj(ρ_ket.data, type = OperatorKet, dims = [4])
    @test_throws DimensionMismatch Qobj(ρ_bra.data, type = OperatorBra, dims = [4])

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
    @test adjoint(a2) ≈ transpose(conj(a2))
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
    @test eigenenergies(a_d * a) ≈ 0:9

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
    @test opstring ==
          "Quantum Object:   type=Operator   dims=$a_dims   size=$a_size   ishermitian=$a_isherm\n$datastring"

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

    ψ2 = Qobj(rand(ComplexF64, 4), type = OperatorKet)
    opstring = sprint((t, s) -> show(t, "text/plain", s), ψ2)
    datastring = sprint((t, s) -> show(t, "text/plain", s), ψ2.data)
    ψ2_dims = ψ2.dims
    ψ2_size = size(ψ2)
    @test opstring == "Quantum Object:   type=OperatorKet   dims=$ψ2_dims   size=$ψ2_size\n$datastring"

    ψ2 = ψ2'
    opstring = sprint((t, s) -> show(t, "text/plain", s), ψ2)
    datastring = sprint((t, s) -> show(t, "text/plain", s), ψ2.data)
    ψ2_dims = ψ2.dims
    ψ2_size = size(ψ2)
    @test opstring == "Quantum Object:   type=OperatorBra   dims=$ψ2_dims   size=$ψ2_size\n$datastring"

    # get coherence
    ψ = coherent(30, 3)
    α, δψ = get_coherence(ψ)
    @test isapprox(abs(α), 3, atol = 1e-5) && abs2(δψ[1]) > 0.999

    # svdvals, Schatten p-norm
    vd = Qobj(rand(ComplexF64, 10))
    vs = Qobj(sprand(ComplexF64, 100, 0.1))
    Md = Qobj(rand(ComplexF64, 10, 10))
    Ms = Qobj(sprand(ComplexF64, 10, 10, 0.5))
    @test svdvals(vd)[1] ≈ √(vd' * vd)
    @test svdvals(vs)[1] ≈ √(vs' * vs)
    @test norm(Md, 1) ≈ sum(sqrt, abs.(eigenenergies(Md' * Md))) atol = 1e-6
    @test norm(Ms, 1) ≈ sum(sqrt, abs.(eigenenergies(Ms' * Ms))) atol = 1e-6

    # trace distance
    ψz0 = basis(2, 0)
    ψz1 = basis(2, 1)
    ρz0 = dense_to_sparse(ket2dm(ψz0))
    ρz1 = dense_to_sparse(ket2dm(ψz1))
    ψx0 = sqrt(0.5) * (basis(2, 0) + basis(2, 1))
    @test tracedist(ψz0, ψx0) ≈ sqrt(0.5)
    @test tracedist(ρz0, ψz1) ≈ 1.0
    @test tracedist(ψz1, ρz0) ≈ 1.0
    @test tracedist(ρz0, ρz1) ≈ 1.0

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

    # data element type conversion
    vd = Qobj(Int64[0, 0])
    vs = Qobj(dense_to_sparse(vd))
    Md = Qobj(Int64[0 0; 0 0])
    Ms = Qobj(dense_to_sparse(Md))
    @test typeof(Vector(vd).data) == Vector{Int64}
    @test typeof(Vector(vs).data) == Vector{Int64}
    @test typeof(Vector{ComplexF64}(vd).data) == Vector{ComplexF64}
    @test typeof(Vector{ComplexF64}(vs).data) == Vector{ComplexF64}
    @test typeof(SparseVector(vd).data) == SparseVector{Int64,Int64}
    @test typeof(SparseVector(vs).data) == SparseVector{Int64,Int64}
    @test typeof(SparseVector{ComplexF64}(vs).data) == SparseVector{ComplexF64,Int64}
    @test typeof(Matrix(Md).data) == Matrix{Int64}
    @test typeof(Matrix(Ms).data) == Matrix{Int64}
    @test typeof(Matrix{ComplexF64}(Ms).data) == Matrix{ComplexF64}
    @test typeof(Matrix{ComplexF64}(Md).data) == Matrix{ComplexF64}
    @test typeof(SparseMatrixCSC(Md).data) == SparseMatrixCSC{Int64,Int64}
    @test typeof(SparseMatrixCSC(Ms).data) == SparseMatrixCSC{Int64,Int64}
    @test typeof(SparseMatrixCSC{ComplexF64}(Ms).data) == SparseMatrixCSC{ComplexF64,Int64}
end
