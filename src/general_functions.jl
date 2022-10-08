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

function meshgrid(x::Union{Vector{T1}, LinRange{T2}}, y::Union{Vector{T1}, LinRange{T2}}) where {T1,T2}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    X, Y
end

"""
    gaussian(x, μ, σ)

Gaussian function.
"""
function gaussian(x::Union{Vector{T1}, LinRange{T2}}, μ::Real, σ::Real) where {T1,T2}
    return exp.(- 0.5 * (x .- μ).^2 / σ::Real^2)
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

function coherent(N::Real, α::T) where T <: Number
    a = destroy(N)
    ad = a'
    return exp(collect(α * ad - α' * a)) * fock(N, 0)
end