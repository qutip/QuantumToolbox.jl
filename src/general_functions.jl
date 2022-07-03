# row_major_reshape(X::AbstractArray, size...) = permutedims(reshape(X, reverse([size...])...), length(size):-1:1)

"""
    row_major_reshape(Q::AbstractArray, shapes)

Reshaping arrays in a different order.
"""
function row_major_reshape(Q::AbstractArray, shapes)
    shapes = reverse(shapes)
    perm = collect(length(shapes):-1:1)
    return permutedims(reshape(Q, shapes), perm)
end

function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]
    X, Y
end

function chop_op(O::AbstractArray, tol = 1e-8)
    tmp_r = (abs.(real.(O)) .> tol) .* real.(O)
    tmp_i = (abs.(imag.(O)) .> tol) .* imag.(O)
    return tmp_r + 1im .* tmp_i
end

"""
    gaussian(x, mu, sig)

Gaussian function.
"""
function gaussian(x, mu, sig)
    return exp.(- 0.5 * (x .- mu).^2 / sig^2)
end

function gaussian_derivative(x, mu, sig)
    return - (x .- mu) ./ sig^2 .* exp.(- 0.5 * (x .- mu).^2 / sig^2)
end

function trunc_op(op::AbstractArray, states)
    N_trunc = size(states)[1]
    # qstates = [qtp.Qobj(states[i], dims = [[N_s, N_s], [1, 1]]) for i in range(N_trunc)]
    res = spzeros(N_trunc, N_trunc)
    for i in range(1, N_trunc, step = 1)
        for j in range(1, N_trunc, step = 1)
            res += (adjoint(states[i, :]) * op * states[j, :]) * projection(N_trunc, i, j, true)
        end
    end
    return chop_op(res)
end

function eigensystem(A::AbstractArray; k = 6, v0 = nothing, sigma = nothing, krylovdim = 30)
    is_A_hermitian = ishermitian(A)

    if v0 === nothing
        v0 = rand(eltype(A), size(A, 1))
    end
    v0 /= norm(v0)

    if sigma === nothing
        vals, vecs, info = eigsolve(A, v0, k, ishermitian = is_A_hermitian)
    else
        # fac = factorize(A - sigma * I)
        # vals, vecs, info = eigsolve(x -> fac \ x, v0, k, ishermitian = is_A_hermitian)
        # vals = (1 .+ sigma * vals) ./ vals
        
        A_s = A - sigma * I

        P_cpu = ilu(A_s, τ = 0.01)
        vals, vecs, info = eigsolve(x -> cg(A_s, x; Pl = P_cpu, maxiter = 500), v0, k, ishermitian = true, krylovdim = krylovdim)
        vals = (1 .+ sigma * vals) ./ vals
    end
    if is_A_hermitian
        vals = real.(vals)
    end
    vals = vals[1:k]
    vecs = hcat(vecs...)[:, 1:k]
    return vals, vecs
end

function ptrace(Q::AbstractArray, sel, rd)
    nd = length(rd)
    dkeep = rd[sel]
    qtrace = filter(e->e∉sel,1:nd)
    dtrace = rd[qtrace]

    if length(size(Q)) == 1 ## Is Ket or Bra
        vmat = row_major_reshape(Q, (prod(rd), 1))
        vmat = row_major_reshape(vmat, rd)
        topermute = []
        append!(topermute, sel)
        append!(topermute, qtrace)
        vmat = permutedims(vmat, topermute)
        vmat = permutedims(vmat, nd:-1:1)
        vmat = row_major_reshape(vmat, (prod(dkeep), prod(dtrace)))
        return vmat * adjoint(vmat)
    elseif length(size(Q)) == 2 ## Is matrix
        ρmat = row_major_reshape(Q, (rd..., rd...))
        topermute = []
        append!(topermute, qtrace)
        append!(topermute, [nd + q for q in qtrace])
        append!(topermute, sel)
        append!(topermute, [nd + q for q in sel])
        ρmat = permutedims(ρmat, topermute)
        ρmat = permutedims(ρmat, 2*nd:-1:1)
        ρmat = row_major_reshape(ρmat, (prod(dtrace), prod(dtrace), prod(dkeep), prod(dkeep)))
        dims = size(ρmat)
        res = [tr(ρmat[:, :, i, j]) for i in 1:dims[3] for j in 1:dims[4]]
        return row_major_reshape(res, dims[3:length(dims)])
    end
end

function entropy_vn(rho::AbstractArray, base = 0, tol = 1e-15)
    vals = eigvals(rho)
    indexes = abs.(vals) .> tol
    if 1 ∈ indexes
        nzvals = vals[indexes]
    #     nzvals = vals[vals .!= 0]
        if base != 0
            logvals = log.(base, Complex.(nzvals))
        else
            logvals = log.(Complex.(nzvals))
        end
        return -real(sum(nzvals .* logvals))
    else
        return 0
    end
end

function entanglement(psi, subspaces, size)
    norm = sqrt(sum(abs.(psi).^2))
    ψ = deepcopy(psi)
    if norm > 1e-5
        ψ = ψ ./ norm
    else
        return 0
    end
    rho_rabi = ptrace(ψ, subspaces, size)
    entropy = entropy_vn(rho_rabi)
    return (entropy > 0) * entropy
end