@testset "Negativity and Partial Transpose" verbose = true begin
    @testset "negativity" begin
        rho = (1 / 40) * Qobj(
            [
                15 1 1 15
                1 5 -3 1
                1 -3 5 1
                15 1 1 15
            ];
            dims = (2, 2),
        )
        Neg = negativity(rho, 1)
        @test Neg ≈ 0.25
        @test negativity(rho, 2) ≈ Neg
        @test negativity(rho, 1; logarithmic = true) ≈ log2(2 * Neg + 1)
        @test_throws ArgumentError negativity(rho, 3)

        @testset "Type Inference (negativity)" begin
            @inferred negativity(rho, 1)
            @inferred negativity(rho, 1; logarithmic = true)
        end
    end

    @testset "partial_transpose" begin
        # A (24 * 24)-matrix which contains number 1 ~ 576
        A_dense = Qobj(reshape(1:(24^2), (24, 24)), dims = (2, 3, 4))
        A_sparse = dense_to_sparse(A_dense)
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

        @testset "Type Inference (partial_transpose)" begin
            @inferred partial_transpose(A_dense, [true, false, true])
            @inferred partial_transpose(A_sparse, [true, false, true])
        end
    end
end
