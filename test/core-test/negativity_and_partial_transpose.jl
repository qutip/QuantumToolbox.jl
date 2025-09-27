@testitem "Negativity and Partial Transpose" begin
    @testset "negativity" begin
        rho1 = (1 / 40) * Qobj(
            [
                15 1 1 15
                1 5 -3 1
                1 -3 5 1
                15 1 1 15
            ];
            dims = (2, 2),
        )
        Neg1 = negativity(rho1, 1)
        @test Neg1 ≈ 0.25
        @test negativity(rho1, 2) ≈ Neg1
        @test negativity(rho1, 1; logarithmic = true) ≈ log2(2 * Neg1 + 1)
        @test_throws ArgumentError negativity(rho1, 3)

        # a maximally entanglment state with subsystem dimension (3, 2):
        # (|1,0⟩ - i|2,1⟩) / √2
        rho2_d = ket2dm((tensor(basis(3, 1), basis(2, 0)) - 1im * tensor(basis(3, 2), basis(2, 1))) / √2)
        rho2_s = to_sparse(rho2_d)
        @test negativity(rho2_d, 1) ≈ 0.5
        @test negativity(rho2_d, 2) ≈ 0.5
        @test negativity(rho2_s, 1) ≈ 0.5
        @test negativity(rho2_s, 2) ≈ 0.5

        # a separable state with subsystem dimension (3, 2)
        rho3_d = tensor(rand_dm(3), rand_dm(2))
        rho3_s = to_sparse(rho3_d)
        @test abs(negativity(rho3_d, 1)) < 1e-10
        @test abs(negativity(rho3_d, 2)) < 1e-10
        @test abs(negativity(rho3_s, 1)) < 1e-10
        @test abs(negativity(rho3_s, 2)) < 1e-10

        @testset "Type Inference (negativity)" begin
            @inferred negativity(rho1, 1)
            @inferred negativity(rho1, 1; logarithmic = true)
        end
    end

    @testset "partial_transpose" begin
        # A (24 * 24)-matrix which contains number 1 ~ 576
        A_dense = Qobj(reshape(1:(24^2), (24, 24)), dims = (2, 3, 4))
        A_sparse = to_sparse(A_dense)
        PT = (true, false)
        for s1 in PT
            for s2 in PT
                for s3 in PT
                    mask = [s1, s2, s3]
                    @test partial_transpose(A_dense, mask) == partial_transpose(A_sparse, mask)
                end
            end
        end
        @test_throws ArgumentError partial_transpose(A_dense, [true])
        @test_throws ArgumentError partial_transpose(Qobj(zeros(ComplexF64, 3, 2)), [true]) # invalid GeneralDimensions

        @testset "Type Inference (partial_transpose)" begin
            @inferred partial_transpose(A_dense, [true, false, true])
            @inferred partial_transpose(A_sparse, [true, false, true])
        end
    end
end
