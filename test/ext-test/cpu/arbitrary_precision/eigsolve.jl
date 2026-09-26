setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox
using LinearAlgebra
using SparseArrays
using GenericSchur

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "Arbitrary Precision (eigsolve)" begin
    L = liouvillian(TESetup.H, TESetup.c_ops)
    L_big = liouvillian(TESetup.H_big, TESetup.c_ops_big)

    # eigenstates(..., sparse=Val(true), ...) directly passes to eigsolve(...)
    vals, vecs = eigenstates(L; sparse = Val(true), sigma = 0.01, eigvals = 7, krylovdim = 30)
    vals_big, vecs_big = eigenstates(L_big; sparse = Val(true), sigma = 0.01, eigvals = 7, krylovdim = 30)

    # Align eigenvalues
    idxs = [findmin(abs.(vals_big .- val))[2] for val in vals]

    @test eltype(vals_big) == Complex{BigFloat}

    @test vals ≈ vals_big[idxs] atol = 1.0e-7
    @test all(zip(vecs, vecs_big[idxs])) do (v, v_big)
        return isapprox(abs(dot(v.data, v_big.data)), 1; atol = 1.0e-7)
    end
end
