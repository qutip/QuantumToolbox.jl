@testitem "linalg helpers" begin
    # meshgrid
    x = [1, 2, 3]
    y = [10, 20]
    X, Y = meshgrid(x, y)
    @test X == [1 2 3; 1 2 3]
    @test Y == [10 10 10; 20 20 20]
end

@testitem "Arnoldi and expv" begin
    import LinearAlgebra: norm

    # Use n = 3 and m = 2 to avoid exact Arnoldi breakdown (and NaN in AS).
    A = 0.01 * [1.0 2.0 0.0; 0.0 3.0 4.0; 5.0 0.0 6.0]
    b = [1.0, 2.0, 3.0]
    m = 2

    AS = arnoldi(A, b, m)
    @test AS.m == m
    @test size(AS.V) == (3, 3)
    @test size(AS.H) == (3, 2)
    @test norm(AS.V[:, 1]) ≈ 1.0

    x = expv(A, 0.2, b; m = m)
    @test x ≈ exp(0.2 * A) * b atol = 1.0e-4

    x2 = zeros(3)
    expv!(x2, AS, 0.2, b)
    @test x2 ≈ exp(0.2 * A) * b atol = 1.0e-4

    # test dimension mismatch
    A_bad = [1.0 0.0; 0.0 1.0]
    b_bad = [1.0, 0.0]
    AS_bad_V = ArnoldiSpace(zeros(3, 3), zeros(3, 2), zeros(3, 2), 2)
    AS_bad_b = ArnoldiSpace(zeros(2, 3), zeros(3, 2), zeros(3, 2), 2)
    @test_throws DimensionMismatch arnoldi!(AS_bad_V, A_bad, b_bad)
    @test_throws DimensionMismatch arnoldi!(AS_bad_b, A_bad, [1.0, 0.0, 0.0])
end
