setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox
using LinearAlgebra
using SparseArrays
using GenericSchur

include("setup.jl") # parameters and operators are defined in this file

@testset "Arbitrary Precision (eigsolve_al)" begin
    L = liouvillian(H, c_ops)
    L_big = liouvillian(H_big, c_ops_big)

    vals_al, vecs_al = eigsolve_al(L, 1 \ (50 * κ), eigvals = 7, krylovdim = 30)
    vals_big_al, vecs_big_al = eigsolve_al(
        L_big,
        1 \ (50 * BigFloat(κ)),
        eigvals = 7,
        krylovdim = 30,
        maxiter = 3,
    )

    # Align eigenvalues
    idxs_al = [findmin(abs.(vals_big_al .- val))[2] for val in vals_al]

    @test eltype(vals_big_al) == Complex{BigFloat}

    @test vals_al[1:(end - 1)] ≈ vals_big_al[idxs_al][1:(end - 1)] atol = 1.0e-6
    @test all(zip(vecs_al[1:(end - 1)], vecs_big_al[idxs_al][1:(end - 1)])) do (v, v_big)
        return isapprox(abs(dot(v.data, v_big.data)), 1; atol = 1.0e-6)
    end
end
