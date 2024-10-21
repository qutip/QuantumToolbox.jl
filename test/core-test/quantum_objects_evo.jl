@testset "Quantum Objects Evolution" verbose = true begin
    # DomainError: incompatible between size of array and type
    @testset "DomainError" begin
        a = MatrixOperator(rand(ComplexF64, 3, 2))
        for t in [Operator, SuperOperator]
            @test_throws DomainError QobjEvo(a, type = t)
        end
    end

    @testset "ArgumentErro" begin
        a = MatrixOperator(rand(ComplexF64, 3, 2))
        for t in (Ket, Bra, OperatorKet, OperatorBra)
            @test_throws ArgumentError QobjEvo(a, type = t)
        end
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
        a = QobjEvo(destroy(N))
        a_d = a'
        X = a + a_d
        # Y = 1im * (a - a_d) # Currently doesn't work. Fix in SciMLOperators.jl
        Z = a + trans(a)
        @test isherm(X) == true
        # @test isherm(Y) == true
        # @test issymmetric(Y) == false
        @test issymmetric(Z) == true
    end

    # TODO: Implement a new show method for QuantumObjectEvolution
    # @testset "REPL show" begin
    #     N = 10
    #     a = QobjEvo(destroy(N))

    #     opstring = sprint((t, s) -> show(t, "text/plain", s), a)
    #     datastring = sprint((t, s) -> show(t, "text/plain", s), a.data)
    #     a_dims = a.dims
    #     a_size = size(a)
    #     a_isherm = isherm(a)
    #     @test opstring ==
    #           "Quantum Object:   type=Operator   dims=$a_dims   size=$a_size   ishermitian=$a_isherm\n$datastring"

    #     a = spre(a)
    #     opstring = sprint((t, s) -> show(t, "text/plain", s), a)
    #     datastring = sprint((t, s) -> show(t, "text/plain", s), a.data)
    #     a_dims = a.dims
    #     a_size = size(a)
    #     a_isherm = isherm(a)
    #     @test opstring == "Quantum Object:   type=SuperOperator   dims=$a_dims   size=$a_size\n$datastring"
    # end

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
    end
end
