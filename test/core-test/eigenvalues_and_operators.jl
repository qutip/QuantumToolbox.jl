@testitem "Eigenvalues" begin
    σx = sigmax()
    result = eigenstates(σx, sparse = Val(false))
    vals_pauli_dense, vecs_pauli_dense, mat_pauli_dense = result
    resstring = sprint((t, s) -> show(t, "text/plain", s), result)
    valstring = sprint((t, s) -> show(t, "text/plain", s), result.values)
    vecsstring = sprint((t, s) -> show(t, "text/plain", s), result.vectors)
    vals_pauli_sparse, vecs_pauli_sparse, mat_pauli_sparse =
        eigenstates(σx, sparse = Val(true), eigvals = 2)

    @test all([ψ.type isa Ket for ψ in vecs_pauli_dense])
    @test mat_pauli_dense isa AbstractMatrix
    @test mat_pauli_sparse isa AbstractMatrix
    @test sort(eigenenergies(σx, sparse = Val(false)); by = real) ≈ sort(vals_pauli_dense; by = real)
    @test sort(eigenenergies(σx, sparse = Val(true), eigvals = 2); by = real) ≈ sort(vals_pauli_sparse; by = real)
    @test resstring ==
        "EigsolveResult:   type=$(Operator())   dims=$(result.dims)\nvalues:\n$(valstring)\nvectors:\n$vecsstring"

    N = 30
    a = kron(destroy(N), qeye(2))
    a_d = a'

    sm = kron(qeye(N), sigmam())
    sp = sm'
    sx = kron(qeye(N), sigmax())
    sy = kron(qeye(N), sigmay())
    sz = kron(qeye(N), sigmaz())

    η = 0.2
    H_d = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx + η^2
    H_c = a_d * a + 0.5 * (sz * cosm(2 * η * (a + a_d)) + sy * sinm(2 * η * (a + a_d)))

    vals_disp_dense, vecs_disp_dense, _ = eigenstates(H_d)
    vals_disp_cosine, vecs_disp_cosine, _ = eigenstates(H_c)
    vals_disp_sparse, vecs_disp_sparse, _ =
        eigenstates(H_d, sparse = Val(true), sigma = -0.9, eigvals = 10, krylovdim = 30, by = real)

    @test real.(vals_disp_dense[1:20]) ≈ real.(vals_disp_cosine[1:20])
    @test real.(vals_disp_dense[1:10]) ≈ real.(vals_disp_sparse[1:10])

    N = 5
    a = kron(destroy(N), qeye(N))
    a_d = a'
    b = kron(qeye(N), destroy(N))
    b_d = b'

    ωc = 1
    ωb = 1
    g = 0.01
    κ = 0.1
    n_th = 0.01

    H = ωc * a_d * a + ωb * b_d * b + g * (a + a_d) * (b + b_d)
    c_ops = [√((1 + n_th) * κ) * a, √κ * b, √(n_th * κ) * a_d]
    L = liouvillian(H, c_ops)

    # eigen solve for general matrices
    vals_liou_dense, _, vecs_liou_dense = eigsolve(L.data, sigma = 0.01, eigvals = 10, krylovdim = 50)
    vals_liou_sparse, _, vecs_liou_sparse = eigenstates(L; sortby = abs)
    vals_liou_sparse = vals_liou_sparse[1:10]
    vecs_liou_sparse = vecs_liou_sparse[:, 1:10]
    vals_liou_al, states_liou_al, vecs_liou_al =
        eigsolve_al(L, 1 \ (10 * κ), eigvals = 10, krylovdim = 50, liouvillian_eigs = Val(true))

    # We now sort the eigenvalues and eigenvectors to make sure they are in the same order.
    # The Arnoldi-Lindblad takes the values closest to the unit circle, not to the origin,
    # so, we only pick those also captured by the standard shift-invert method
    sort_func = x -> (round(abs(x), digits = 6), round(real(x), digits = 6), round(imag(x), digits = 6))
    idxs_dense = sortperm(vals_liou_dense, by = sort_func, rev = false)
    idxs_sparse = sortperm(vals_liou_sparse, by = sort_func, rev = false)
    idxs_al = sortperm(vals_liou_al, by = sort_func, rev = false)
    vals_liou_dense = vals_liou_dense[idxs_dense]
    vecs_liou_dense = vecs_liou_dense[:, idxs_dense]
    vals_liou_sparse = vals_liou_sparse[idxs_sparse]
    vecs_liou_sparse = vecs_liou_sparse[:, idxs_sparse]
    vals_liou_al = vals_liou_al[idxs_al]
    states_liou_al = states_liou_al[idxs_al]
    vecs_liou_al = vecs_liou_al[:, idxs_al]

    @test vals_liou_dense ≈ vals_liou_sparse
    @test vals_liou_sparse[1:2] ≈ vals_liou_al[1:2] atol = 1.0e-6
    @test vec2mat(vecs_liou_dense[:, 1]) * exp(-1im * angle(vecs_liou_dense[1, 1])) ≈
        vec2mat(vecs_liou_sparse[:, 1]) * exp(-1im * angle(vecs_liou_sparse[1, 1])) atol = 1.0e-7
    @test vec2mat(vecs_liou_dense[:, 1]) * exp(-1im * angle(vecs_liou_dense[1, 1])) ≈
        vec2mat(vecs_liou_al[:, 1]) * exp(-1im * angle(vecs_liou_al[1, 1])) atol = 1.0e-5

    # eigen solve for QuantumObject
    result_sparse = eigenstates(L, sparse = Val(true), sigma = 0.01, eigvals = 10, krylovdim = 50)
    vals_sparse_qobj, vecs_sparse_qobj = result_sparse
    resstring = sprint((t, s) -> show(t, "text/plain", s), result_sparse)
    valstring = sprint((t, s) -> show(t, "text/plain", s), result_sparse.values)
    vecsstring = sprint((t, s) -> show(t, "text/plain", s), result_sparse.vectors)
    @test resstring ==
        "EigsolveResult:   type=$(SuperOperator())   dims=$(result_sparse.dims)\nvalues:\n$(valstring)\nvectors:\n$vecsstring"

    vals_qobj_sorted, vecs_qobj_sorted = eigenstates(L, sortby = abs)
    vals_qobj_sorted = vals_qobj_sorted[1:10]
    vecs_qobj_sorted = vecs_qobj_sorted[1:10]

    ρss_dense = vector_to_operator(vecs_sparse_qobj[1])
    ρss_qobj = vector_to_operator(vecs_qobj_sorted[1])
    ρss_al = vector_to_operator(states_liou_al[1])
    ρss_dense /= tr(ρss_dense)
    ρss_qobj /= tr(ρss_qobj)
    ρss_al /= tr(ρss_al)

    @test result_sparse.type isa SuperOperator
    @test result_sparse.dims == L.dims
    @test all([v.type isa OperatorKet for v in vecs_sparse_qobj])
    @test result_sparse.vectors isa AbstractMatrix
    @test sum(abs2, vals_sparse_qobj) ≈ sum(abs2, vals_qobj_sorted)
    @test vals_qobj_sorted[1:2] ≈ vals_liou_al[1:2]
    @test fidelity(ρss_dense, ρss_qobj) ≈ 1
    @test fidelity(ρss_dense, ρss_al) ≈ 1

    @testset "Type Inference (eigen)" begin
        UnionType = Union{
            QuantumToolbox.EigsolveResult{
                Vector{ComplexF64},
                Matrix{ComplexF64},
                QuantumToolbox.Operator,
                QuantumToolbox.Dimensions{2, Tuple{QuantumToolbox.Space, QuantumToolbox.Space}},
            },
            QuantumToolbox.EigsolveResult{
                Vector{Float64},
                Matrix{ComplexF64},
                QuantumToolbox.Operator,
                QuantumToolbox.Dimensions{2, Tuple{QuantumToolbox.Space, QuantumToolbox.Space}},
            },
        }

        @inferred UnionType eigenstates(H, sparse = Val(false))
        @inferred eigenstates(H, sparse = Val(true))
        @inferred eigenstates(L, sparse = Val(true))
        @inferred eigsolve_al(L, 1 \ (10 * κ), eigvals = 10, liouvillian_eigs = Val(true))
    end
end
