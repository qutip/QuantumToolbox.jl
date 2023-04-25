"""
    row_major_reshape(Q::AbstractArray, shapes...)

Reshapes `Q` in the row-major order, as numpy. 
"""
row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} = PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner=length(y)), length(y), length(x))
    Y = repeat(y, outer=(1, length(x)))
    X, Y
end

@doc raw"""
    gaussian(x, μ::Real, σ::Real)

Returns the gaussian function ``\exp \left[- \frac{(x - \mu)^2}{2 \sigma^2} \right]``,
where ``\mu`` and ``\sigma^2`` are the mean and the variance respectively.
"""
gaussian(x::AbstractVector{T}, μ::Real, σ::Real) where {T} = @. exp(-0.5 * (x - μ)^2 / σ^2)

gaussian(x::Real, μ::Real, σ::Real) = exp(-0.5 * (x - μ)^2 / σ^2)

@doc raw"""
    ptrace(QO::QuantumObject, sel::Vector{Int})

[Partial trace](https://en.wikipedia.org/wiki/Partial_trace) of a quantum state `QO` leaving only the dimensions
with the indices present in the `sel` vector.

# Examples
Two qubits in the state ``\ket{\psi} = \ket{e,g}``:
```jldoctest; setup=(using QuPhys)
julia> ψ = kron(fock(2,0), fock(2,1))
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> ptrace(ψ, [2])
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im
```

or in an entangled state ``\ket{\psi} = \frac{1}{\sqrt{2}} \left( \ket{e,e} + \ket{g,g} \right)``:
```jldoctest; setup=(using QuPhys)
julia> ψ = 1 / √2 * (kron(fock(2,0), fock(2,0)) + kron(fock(2,1), fock(2,1)))
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> ptrace(ψ, [1])
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im
```
"""
function ptrace(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::Vector{Int}) where
{T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    rd = QO.dims
    nd = length(rd)
    dkeep = rd[sel]
    qtrace = filter(e -> e ∉ sel, 1:nd)
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

@doc raw"""
    entropy_vn(ρ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy)
``S = - \Tr \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right]`` where ``\hat{\rho}``
is the density matrix of the system.

The `base` parameter specifies the base of the logarithm to use, and when using the default value 0,
the natural logarithm is used. The `tol` parameter
describes the absolute tolerance for detecting the zero-valued eigenvalues of the density
matrix ``\hat{\rho}``.

# Examples

Pure state:
```jldoctest; setup=(using QuPhys)
julia> ψ = fock(2,0)
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im

julia> ρ = ket2dm(ψ)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> entropy_vn(ρ, base=2)
-0.0
```

Mixed state:
```jldoctest; setup=(using QuPhys)
julia> ρ = 1 / 2 * ( ket2dm(fock(2,0)) + ket2dm(fock(2,1)) )
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im

julia> entropy_vn(ρ, base=2)
1.0
```
"""
function entropy_vn(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}; base::Int=0, tol::Real=1e-15) where {T}
    vals = eigvals(ρ)
    indexes = abs.(vals) .> tol
    1 ∉ indexes && return 0 
    nzvals = vals[indexes]
    logvals = base != 0 ? log.(base, Complex.(nzvals)) : log.(Complex.(nzvals))
    return -real(sum(nzvals .* logvals))
end

"""
    entanglement(QO::QuantumObject, sel::Vector)

Calculates the entanglement by doing the partial trace of `QO`, selecting only the dimensions
with the indices contained in the `sel` vector, and then using the Von Neumann entropy [`entropy_vn`](@ref). 
"""
function entanglement(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::Vector{Int}) where
{T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    ψ = normalize(QO)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end

@doc raw"""
    expect(O::QuantumObject, ψ::QuantumObject)

Expectation value of the operator `O` with the state `ψ`. The latter
can be a [`KetQuantumObject`](@ref), [`BraQuantumObject`](@ref) or a [`OperatorQuantumObject`](@ref).
If `ψ` is a density matrix, the function calculates ``\Tr \left[ \hat{O} \hat{\psi} \right]``, while if `ψ`
is a state, the function calculates ``\mel{\psi}{\hat{O}}{\psi}``.

The function returns a real number if the operator is hermitian, and returns a complex number otherwise.

# Examples

```jldoctest; setup=(using QuPhys)
julia> ψ = 1 / √2 * (fock(10,2) + fock(10,4));

julia> a = destroy(10);

julia> expect(a' * a, ψ) ≈ 3
true
```
"""
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T2},KetQuantumObject}) where {T1,T2}
    ψd = ψ.data
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od * ψd)) : dot(ψd, Od * ψd)
end
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T2},BraQuantumObject}) where {T1,T2}
    ψd = ψ.data'
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od * ψd)) : dot(ψd, Od * ψd)
end
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ρ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject}) where {T1,T2}
    return ishermitian(O) ? real(tr(O * ρ)) : tr(O * ρ)
end

@doc raw"""
    get_coherence(ψ::QuantumObject)

Get the coherence value ``\alpha`` by measuring the expectation value of the destruction
operator ``\hat{a}`` on the state ``\ket{\psi}``.

It returns both ``\alpha`` and the state 
``\ket{\delta_\psi} = \exp ( \bar{\alpha} \hat{a} - \alpha \hat{a}^\dagger )``. The
latter corresponds to the quantum fulctuations around the coherent state ``\ket{\alpha}``.
"""
function get_coherence(ψ::QuantumObject{<:AbstractArray{T}, StateOpType}) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    a = destroy(size(ψ,1))
    α = expect(a, ψ)
    D = exp(α*a' - conj(α)*a)

    α, D' * ψ
end


@doc raw"""
    wigner(state::QuantumObject, xvec::AbstractVector, yvec::AbstractVector; g::Real=√2)

Generates the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution)
of `state` at points `xvec + 1im * yvec`. The `g` parameter is a scaling factor related to the value of ``\hbar`` in the
commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2`` giving the default value ``\hbar=1``.
"""
function wigner(state::QuantumObject{<:AbstractArray{T1},OpType}, xvec::AbstractVector{T2},
    yvec::AbstractVector{T2}; g::Real=√2) where {T1,T2,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    if isket(state)
        ρ = (state * state').data
    elseif isbra(state)
        ρ = (state' * state).data
    else
        ρ = state.data
    end
    M = size(ρ, 1)
    X, Y = meshgrid(xvec, yvec)
    A2 = g * (X + 1im * Y)

    B = abs.(A2)
    B .*= B
    w0 = (2 * ρ[1, end]) .* ones(eltype(A2), size(A2)...)
    L = M - 1

    ρ = ρ .* (2 * ones(M, M) - diagm(ones(M)))
    while L > 0
        L -= 1
        w0 = _wig_laguerre_val(L, B, diag(ρ, L)) .+ w0 .* A2 .* (L + 1)^(-0.5)
    end

    return @. real(w0) * exp(-B * 0.5) * (g * g * 0.5 / π)
end

function _wig_laguerre_val(L, x, c)
    if length(c) == 1
        y0 = c[1]
        y1 = 0
    elseif length(c) == 2
        y0 = c[1]
        y1 = c[2]
    else
        k = length(c)
        y0 = c[end-1]
        y1 = c[end]
        for i in range(3, length(c), step=1)
            k -= 1
            y0, y1 = @. c[end+1-i] - y1 * ((k - 1) * (L + k - 1) / ((L + k) * k))^0.5, y0 - y1 * ((L + 2 * k - 1) - x) * ((L + k) * k)^(-0.5)
        end
    end
    return @. y0 - y1 * ((L + 1) - x) * (L + 1)^(-0.5)
end