@testset "Quantum Objects Evolution" verbose = true begin
    # DomainError: incompatible between size of array and type
    @testset "Thrown Errors" begin
        a = MatrixOperator(rand(ComplexF64, 3, 2))
        for t in [Operator, SuperOperator]
            @test_throws DomainError QobjEvo(a, type = t)
        end

        a = MatrixOperator(rand(ComplexF64, 3, 2))
        for t in (Ket, Bra, OperatorKet, OperatorBra)
            @test_throws ArgumentError QobjEvo(a, type = t)
        end

        a = QobjEvo(destroy(20))
        @test_throws ArgumentError QobjEvo(a, type = SuperOperator)

        a = MatrixOperator(rand(ComplexF64, 5, 5))
        @test_throws DimensionMismatch QobjEvo(a, type = SuperOperator)

        ψ = fock(10, 3)
        @test_throws TypeError QobjEvo(ψ)
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
            "The argument dims should be a Tuple or a StaticVector for better performance. Try to use `dims = (2, 2)` or `dims = SVector(2, 2)` instead of `dims = [2, 2]`.",
        ) QobjEvo(MatrixOperator(rand(4, 4)), dims = [2, 2])
    end

    @testset "Operator and SuperOperator" begin
        a = MatrixOperator(sprand(ComplexF64, 100, 100, 0.1))
        a2 = QobjEvo(a)
        a3 = QobjEvo(a, type = SuperOperator)

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
        a = MatrixOperator(sprand(ComplexF64, 100, 100, 0.1))
        a2 = QobjEvo(a)
        a3 = QobjEvo(a, type = SuperOperator)

        @test +a2 == a2
        @test -(-a2) == a2
        @test a2 + 2 == 2 + a2
        @test (a2 + 2).data == a2.data + 2 * I
        @test a2 * 2 == 2 * a2

        @test trans(trans(a2)) == a2
        @test trans(a2).data == transpose(a2.data)
        # @test adjoint(a2) ≈ trans(conj(a2)) # Currently doesn't work
        @test adjoint(adjoint(a2)) == a2
        @test adjoint(a2).data == adjoint(a2.data)

        N = 10
        # We use MatrixOperator instead of directly using a Qobj to increase coverage
        a = QobjEvo(MatrixOperator(sprand(ComplexF64, N, N, 5 / N)), Operator, N)
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
              "Quantum Object Evo.:   type=Operator   dims=$H_dims   size=$H_size   ishermitian=$H_isherm   isconstant=$H_isconst\n$datastring"

        L = QobjEvo(spre(a))
        opstring = sprint((t, s) -> show(t, "text/plain", s), L)
        datastring = sprint((t, s) -> show(t, "text/plain", s), L.data)
        L_dims = L.dims
        L_size = size(L)
        L_isherm = isherm(L)
        L_isconst = isconstant(L)
        @test opstring ==
              "Quantum Object Evo.:   type=SuperOperator   dims=$L_dims   size=$L_size   ishermitian=$L_isherm   isconstant=$L_isconst\n$datastring"
    end

    @testset "Type Inference (QuantumObject)" begin
        for T in [ComplexF32, ComplexF64]
            N = 4
            a = MatrixOperator(rand(T, N, N))
            @inferred QobjEvo(a)
            for type in [Operator, SuperOperator]
                @inferred QobjEvo(a, type = type)
            end
        end

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

            # TODO: kron is currently not supported
            # @inferred kron(a)
            # @inferred kron(a, σx)
            # @inferred kron(a, eye(2))
        end
    end

    # TODO: tensor is currently not supported
    # @testset "tensor" begin
    #     σx = sigmax()
    #     X3 = kron(σx, σx, σx)
    #     @test tensor(σx) == kron(σx)
    #     @test tensor(fill(σx, 3)...) == X3
    #     X_warn = @test_logs (
    #         :warn,
    #         "`tensor(A)` or `kron(A)` with `A` is a `Vector` can hurt performance. Try to use `tensor(A...)` or `kron(A...)` instead.",
    #     ) tensor(fill(σx, 3))
    #     @test X_warn == X3
    # end

    @testset "Time Dependent Operators" begin
        N = 10
        a = destroy(N)
        coef1(p, t) = exp(-1im * p.ω1 * t)
        coef2(p, t) = sin(p.ω2 * t)

        @test_throws MethodError QobjEvo([[a, coef1], a' * a, [a', coef2]])

        op1 = QobjEvo(((a, coef1), a' * a, (a', coef2)))
        op1 = QobjEvo(((a, coef1), a' * a, (a', coef2)))

        p = (ω1 = 1, ω2 = 2)
        @test op1(p, 0.1) ≈ coef1(p, 0.1) * a + a' * a + coef2(p, 0.1) * a'

        ψ = fock(N, 1)
        @test op1(ψ, p, 0.1) ≈ (coef1(p, 0.1) * a + a' * a + coef2(p, 0.1) * a') * ψ

        @test isconstant(a) == true
        @test isconstant(op1) == false
        @test isconstant(Qobj(a)) == true
    end
end
