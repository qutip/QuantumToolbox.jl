row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} = PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    X, Y
end

function gaussian(x::AbstractVector{T}, μ::Real, σ::Real) where {T}
    return exp.(- 0.5 * (x .- μ).^2 / σ^2)
end

function ptrace(QO::QuantumObject{<:AbstractArray{T}, OpType}, sel) where 
        {T,OpType<:Union{BraQuantumObject, KetQuantumObject, OperatorQuantumObject}}
        
    rd = QO.dims
    nd = length(rd)
    dkeep = rd[sel]
    qtrace = filter(e->e∉sel,1:nd)
    dtrace = rd[qtrace]

    nd == 1 && return QO

    if isket(QO) || isbra(QO)
        vmat = row_major_reshape(QO.data, prod(rd), 1)
        vmat = row_major_reshape(vmat, rd...)
        topermute = Int64[]
        append!(topermute, sel)
        append!(topermute, qtrace)
        reverse!(topermute)
        vmat = PermutedDimsArray(vmat, topermute)
        vmat = row_major_reshape(vmat, prod(dkeep), prod(dtrace))
        return QuantumObject(vmat * vmat', OperatorQuantumObject, dkeep)
    elseif isoper(QO)
        ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
        topermute = Int64[]
        append!(topermute, qtrace)
        append!(topermute, [nd + q for q in qtrace])
        append!(topermute, sel)
        append!(topermute, [nd + q for q in sel])
        reverse!(topermute)
        ρmat = PermutedDimsArray(ρmat, topermute)
        ρmat = row_major_reshape(ρmat, prod(dtrace), prod(dtrace), prod(dkeep), prod(dkeep))
        dims = size(ρmat)
        res = [tr(@view(ρmat[:, :, i, j])) for i in 1:dims[3] for j in 1:dims[4]]
        return QuantumObject(row_major_reshape(res, dims[3:length(dims)]...), OperatorQuantumObject, dkeep)
    end
end

function entropy_vn(ρ::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, base = 0, tol = 1e-15) where {T}
    vals = eigvals(ρ)
    indexes = abs.(vals) .> tol
    if 1 ∈ indexes
        nzvals = vals[indexes]
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

function entanglement(psi::QuantumObject{<:AbstractArray{T}, OpType}, sel) where 
        {T,OpType<:Union{BraQuantumObject, KetQuantumObject, OperatorQuantumObject}}

    ψ = normalize(psi)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end

function expect(op::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T}, KetQuantumObject}) where {T}
    ishermitian(op) && return real(ψ' * op * ψ)
    return ψ' * op * ψ
end
function expect(op::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T}, BraQuantumObject}) where {T}
    ishermitian(op) && return real(ψ * op * ψ')
    return ψ * op * ψ'
end
function expect(op::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, ρ::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T}
    ishermitian(op) && return real(tr(op * ρ))
    return tr(op * ρ)
end

function wigner(state::QuantumObject{<:AbstractArray{T1}, OpType}, xvec::AbstractVector{T2}, 
    yvec::AbstractVector{T2}; g::Real = √2) where {T1,T2,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    if isket(state)
        ρ = (state * state').data
    elseif isbra(state)
        ρ = (state' * state).data
    else
        ρ = state.data
    end
    M = size(ρ, 1)
    X,Y = meshgrid(xvec, yvec)
    A2 = g * (X + 1im * Y)

    B = abs.(A2)
    B .*= B
    w0 = (2*ρ[1, end]) .* ones(eltype(A2), size(A2)...)
    L = M-1

    ρ = ρ .* (2 * ones(M,M) - diagm(ones(M)))
    while L > 0
        L -= 1
        w0 = _wig_laguerre_val(L, B, diag(ρ, L)) .+ w0 .* A2 .* (L+1)^(-0.5)
    end

    return @. real(w0) * exp(-B * 0.5) * (g*g*0.5 / π)
end

function _wig_laguerre_val(L, x, c)
    if length(c) == 1
        y0 = c[1]
        y1 = 0
    elseif  length(c) == 2
        y0 = c[1]
        y1 = c[2]
    else
        k = length(c)
        y0 = c[end-1]
        y1 = c[end]
        for i in range(3, length(c), step = 1)
            k -= 1
            y0, y1 = @. c[end+1-i] - y1 * ( (k - 1)*(L+k-1)/((L + k)*k) )^0.5, y0 - y1 * ((L + 2*k -1) - x) * ((L+k)*k)^(-0.5)
        end
    end
    return @. y0 - y1 * ((L + 1) - x) * (L + 1)^(-0.5)
end