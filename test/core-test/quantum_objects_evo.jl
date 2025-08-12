@testitem "Quantum Objects Evolution" begin
    using LinearAlgebra
    using SparseArrays
    using StaticArraysCore
    using SciMLOperators

    # DomainError: incompatible between size of array and type
    @testset "Thrown Errors" begin
        a = MatrixOperator(rand(ComplexF64, 3, 2))
        @test_throws DomainError QobjEvo(a, type = SuperOperator())

        a = MatrixOperator(rand(ComplexF64, 4, 4))
        @test_throws DomainError QobjEvo(a, type = SuperOperator(), dims = ((2,), (2,)))

        a = MatrixOperator(rand(ComplexF64, 3, 2))
        for t in (Ket(), Bra(), OperatorKet(), OperatorBra())
            @test_throws ArgumentError QobjEvo(a, type = t)
        end

        a = QobjEvo(destroy(20))
        @test_throws ArgumentError QobjEvo(a, type = SuperOperator())

        a = MatrixOperator(rand(ComplexF64, 5, 5))
        @test_throws DimensionMismatch QobjEvo(a, type = SuperOperator())

        ψ = fock(10, 3)
        @test_throws MethodError QobjEvo(ψ)
    end

    # unsupported type of dims
    @testset "unsupported dims" begin
        a = MatrixOperator(rand(2, 2))
        @test_throws ArgumentError QobjEvo(a, dims = 2.0)
        @test_throws ArgumentError QobjEvo(a, dims = 2.0 + 0.0im)
        @test_throws DomainError QobjEvo(a, dims = 0)
        @test_throws DomainError QobjEvo(a, dims = (2, -2))
        @test_logs (
            :warn,
            "The argument dims should be a Tuple or a StaticVector for better performance. Try to use `dims = (2, 2)` instead of `dims = [2, 2]`. Alternatively, you can do `import QuantumToolbox: SVector` and use `dims = SVector(2, 2)`.",
        ) QobjEvo(MatrixOperator(rand(4, 4)), dims = [2, 2])
    end

    @testset "Operator and SuperOperator" begin
        a = MatrixOperator(sprand(ComplexF64, 100, 100, 0.1))
        a2 = QobjEvo(a)
        a3 = QobjEvo(a, type = SuperOperator())

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
        @test_throws DimensionMismatch QobjEvo(a, dims = 2)
    end

    @testset "Promote Operators Type" begin
        a = destroy(20)
        A = QobjEvo(a)
        @test QuantumToolbox.promote_op_type(a, A) == QuantumObjectEvolution
        @test QuantumToolbox.promote_op_type(A, a) == QuantumObjectEvolution
        @test QuantumToolbox.promote_op_type(A, A) == QuantumObjectEvolution
        @test QuantumToolbox.promote_op_type(a, a) == QuantumObject
    end

    @testset "arithmetic" begin
        m = sprand(ComplexF64, 100, 100, 0.1)
        a = MatrixOperator(m)
        a2 = QobjEvo(a)
        a3 = QobjEvo(a, type = SuperOperator())
        a4 = QobjEvo(Qobj(m), (p, t) -> 1)

        @test +a2 == a2
        @test -(-a2) == a2
        @test a2 + 2 == 2 + a2
        @test (a2 + 2).data == a2.data + 2 * I
        @test a2 * 2 == 2 * a2
        @test 1*a*a*1 == 1*(a*a) == (a*a)*1 == 1*(a*a)*1 # to check QobjEvo multiplication with Number is correct

        zero_like = zero(a2)
        iden_like = one(a3)
        zero_array = NullOperator(100)
        iden_array = IdentityOperator(100)
        @test zero_like == QobjEvo(zero_array, type = a2.type, dims = a2.dims)
        @test typeof(zero_like.data) == typeof(zero_array)
        @test iden_like == QobjEvo(iden_array, type = a3.type, dims = a3.dims)
        @test typeof(iden_like.data) == typeof(iden_array)
        @test trans(trans(a2)) == a2
        @test trans(a2).data == transpose(a2.data)
        # @test adjoint(a2) ≈ trans(conj(a2)) # Currently doesn't work
        @test adjoint(adjoint(a2)) == a2
        @test adjoint(a2).data == adjoint(a2.data)

        N = 10
        # We use MatrixOperator instead of directly using a Qobj to increase coverage
        a = QobjEvo(MatrixOperator(sprand(ComplexF64, N, N, 5 / N)), Operator(), N)
        a_d = a'
        X = a + a_d
        # Y = 1im * (a - a_d) # Currently doesn't work. Fix in SciMLOperators.jl
        Z = a + trans(a)
        @test isherm(X) == true
        # @test isherm(Y) == true
        # @test issymmetric(Y) == false
        @test issymmetric(Z) == true
    end

    @testset "REPL show" begin
        N = 10
        a = destroy(N)
        coef(p, t) = exp(-1im * t)
        H = QobjEvo((a' * a, (a, coef)))

        opstring = sprint((t, s) -> show(t, "text/plain", s), H)
        datastring = sprint((t, s) -> show(t, "text/plain", s), H.data)
        H_dims = H.dims
        H_size = size(H)
        H_isherm = isherm(H)
        H_isconst = isconstant(H)
        @test opstring ==
              "\nQuantum Object Evo.:   type=Operator()   dims=$H_dims   size=$H_size   ishermitian=$H_isherm   isconstant=$H_isconst\n$datastring"

        L = QobjEvo(spre(a))
        opstring = sprint((t, s) -> show(t, "text/plain", s), L)
        datastring = sprint((t, s) -> show(t, "text/plain", s), L.data)
        L_dims = L.dims
        L_size = size(L)
        L_isherm = isherm(L)
        L_isconst = isconstant(L)
        @test opstring ==
              "\nQuantum Object Evo.:   type=SuperOperator()   dims=$L_dims   size=$L_size   ishermitian=$L_isherm   isconstant=$L_isconst\n$datastring"
    end

    @testset "Type Inference (QobjEvo)" begin
        N = 4
        for T in [ComplexF32, ComplexF64]
            a = MatrixOperator(rand(T, N, N))
            UnionType = Union{
                QuantumObjectEvolution{Operator,GeneralDimensions{1,Tuple{Space},Tuple{Space}},typeof(a)},
                QuantumObjectEvolution{Operator,Dimensions{1,Tuple{Space}},typeof(a)},
            }
            @inferred UnionType QobjEvo(a)
            @inferred UnionType QobjEvo(a, type = Operator())
            @inferred QobjEvo(a, type = SuperOperator())
        end

        a = destroy(N)
        coef1(p, t) = exp(-t)
        coef2(p::Vector, t) = sin(p[1] * t)
        coef3(p::NamedTuple, t) = cos(p.ω * t)
        @inferred QobjEvo(a, coef1)
        @inferred QobjEvo((a', coef2))
        @inferred QobjEvo((a' * a, (a, coef1), (a', coef2), (a + a', coef3)))

        @testset "Math Operation" begin
            a = QobjEvo(destroy(20))
            σx = QobjEvo(sigmax())
            @inferred a + a
            @inferred a + a'
            # @inferred a + 2 # TODO fix in SciMLOperators.jl
            @inferred 2 * a
            @inferred a / 2
            @inferred a * a
            @inferred a * a'

            @inferred kron(a)
            @test_logs (:warn,) @inferred kron(a, σx)
            @test_logs (:warn,) @inferred kron(a, eye(2))
            @test_logs (:warn,) (:warn,) @inferred kron(a, eye(2), eye(2))
        end
    end

    @testset "tensor" begin
        σx = QobjEvo(sigmax())
        X3 = @test_logs (:warn,) (:warn,) tensor(σx, σx, σx)
        X_warn = @test_logs (:warn,) (:warn,) (:warn,) tensor(fill(σx, 3))
        @test X_warn(0) == X3(0) == tensor(sigmax(), sigmax(), sigmax())
    end

    @testset "Time Dependent Operators and SuperOperators" begin
        N = 10
        a = destroy(N)
        coef1(p, t) = exp(-1im * p.ω1 * t)
        coef2(p, t) = sin(p.ω2 * t)
        coef3(p, t) = sin(p.ω3 * t)
        t = rand()
        p = (ω1 = rand(), ω2 = rand(), ω3 = rand())

        # Operator
        H_td = QobjEvo(((a, coef1), a' * a, (a', coef2)))
        H_ti = coef1(p, t) * a + a' * a + coef2(p, t) * a'
        ψ = rand_ket(N)
        @test H_td(p, t) ≈ H_ti
        @test iscached(H_td) == true
        H_td = cache_operator(H_td, ψ)
        @test iscached(H_td) == true
        @test H_td(ψ, p, t) ≈ H_ti * ψ
        @test isconstant(a) == true
        @test isconstant(H_td) == false
        @test isconstant(QobjEvo(a)) == true
        @test isoper(H_td) == true
        @test QobjEvo(a, coef1) == QobjEvo((a, coef1))

        # SuperOperator
        X = a * a'
        c_op1 = QobjEvo(a', coef1)
        c_op2 = QobjEvo(((a, coef2), (X, coef3)))
        c_ops = [c_op1, c_op2]
        D1_ti = abs2(coef1(p, t)) * lindblad_dissipator(a')
        D2_ti =
            abs2(coef2(p, t)) * lindblad_dissipator(a) + # normal dissipator for first  element in c_op2
            abs2(coef3(p, t)) * lindblad_dissipator(X) + # normal dissipator for second element in c_op2
            coef2(p, t) * conj(coef3(p, t)) * (spre(a) * spost(X') - 0.5 * spre(X' * a) - 0.5 * spost(X' * a)) + # cross terms
            conj(coef2(p, t)) * coef3(p, t) * (spre(X) * spost(a') - 0.5 * spre(a' * X) - 0.5 * spost(a' * X))   # cross terms
        L_ti = liouvillian(H_ti) + D1_ti + D2_ti
        L_td = @test_logs (:warn,) (:warn,) liouvillian(H_td, c_ops) # warnings from lazy tensor in `lindblad_dissipator(c_op2)`
        ρvec = mat2vec(rand_dm(N))
        @test L_td(p, t) ≈ L_ti
        @test iscached(L_td) == false
        L_td = cache_operator(L_td, ρvec)
        @test iscached(L_td) == true
        @test L_td(ρvec, p, t) ≈ L_ti * ρvec
        @test isconstant(L_td) == false
        @test issuper(L_td) == true

        coef_wrong1(t) = nothing
        coef_wrong2(p, t::ComplexF64) = nothing
        @test_logs (:warn,) (:warn,) liouvillian(H_td * H_td) # warnings from lazy tensor
        @test_throws ArgumentError QobjEvo(a, coef_wrong1)
        @test_throws ArgumentError QobjEvo(a, coef_wrong2)
        @test_throws MethodError QobjEvo([[a, coef1], a' * a, [a', coef2]])
        @test_throws ArgumentError H_td(ρvec, p, t)
        @test_throws ArgumentError cache_operator(H_td, ρvec)
        @test_throws ArgumentError L_td(ψ, p, t)
        @test_throws ArgumentError cache_operator(L_td, ψ)

        @testset "Type Inference" begin
            # we use destroy and create here because they somehow causes type instability before
            H_td2 = H_td + QobjEvo(destroy(N) + create(N), coef3)
            c_ops1 = (destroy(N), create(N))
            c_ops2 = (destroy(N), QobjEvo(create(N), coef1))

            @inferred liouvillian(H_td, c_ops1)
            @inferred liouvillian(H_td, c_ops2)
            @inferred liouvillian(H_td2, c_ops1)
            @inferred liouvillian(H_td2, c_ops2)
        end
    end
end
