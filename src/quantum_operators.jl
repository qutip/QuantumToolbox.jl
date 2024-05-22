export spre, spost, sprepost, lindblad_dissipator
export fock, basis, coherent
export sigmam, sigmap, sigmax, sigmay, sigmaz
export destroy, create, eye, qeye, projection, rand_dm
export sinm, cosm

@doc raw"""
    spre(O::QuantumObject, Id_cache=I(size(O,1)))

Returns the super-operator form of `O` acting on the left
of the density matrix operator, ``\mathcal{O} \left(\hat{O}\right) \left[ \hat{\rho} \right] = \hat{O} \hat{\rho}``.

Since the density matrix is vectorized, this super-operator is always
a matrix, obtained from ``\mathcal{O} \left(\hat{O}\right) \boldsymbol{\cdot} = \hat{\mathbb{1}} \otimes \hat{O}``.

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
spre(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, Id_cache = I(size(O, 1))) where {T} =
    QuantumObject(kron(Id_cache, O.data), SuperOperator, O.dims)

@doc raw"""
    spost(O::QuantumObject)

Returns the super-operator form of `O` acting on the right
of the density matrix operator, ``\mathcal{O} \left(\hat{O}\right) \left[ \hat{\rho} \right] = \hat{\rho} \hat{O}``.

Since the density matrix is vectorized, this super-operator is always
a matrix, obtained from ``\mathcal{O} \left(\hat{O}\right) \boldsymbol{\cdot} = \hat{O}^T \otimes \hat{\mathbb{1}}``.

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
spost(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, Id_cache = I(size(O, 1))) where {T} =
    QuantumObject(kron(sparse(transpose(sparse(O.data))), Id_cache), SuperOperator, O.dims) # TODO: fix the sparse conversion

@doc raw"""
    sprepost(A::QuantumObject, B::QuantumObject)

Returns the super-operator form of `A` and `B` acting on the left and the right
of the density matrix operator respectively, ``\mathcal{O} \left( \hat{A}, \hat{B} \right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho} \hat{B}``.

Since the density matrix is vectorized, this super-operator is always
a matrix, obtained from ``\mathcal{O} \left(\hat{A}, \hat{B}\right) \boldsymbol{\cdot} = \text{spre}(A) * \text{spost}(B)``.
"""
sprepost(
    A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2} = QuantumObject(kron(sparse(transpose(sparse(B.data))), A.data), SuperOperator, A.dims) # TODO: fix the sparse conversion

@doc raw"""
    lindblad_dissipator(O::QuantumObject, Id_cache=I(size(O,1))

Returns the Lindblad super-operator defined as
``
\mathcal{D} \left( \hat{O} \right) \left[ \hat{\rho} \right] = \frac{1}{2} \left( 2 \hat{O} \hat{\rho} \hat{O}^\dagger - 
\hat{O}^\dagger \hat{O} \hat{\rho} - \hat{\rho} \hat{O}^\dagger \hat{O} \right)
``
considering the density matrix ``\hat{\rho}`` in the vectorized form.

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
function lindblad_dissipator(
    O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject},
    Id_cache = I(size(O, 1)),
) where {T}
    Od_O = O' * O
    return sprepost(O, O') - spre(Od_O, Id_cache) / 2 - spost(Od_O, Id_cache) / 2
end

# It is already a SuperOperator
lindblad_dissipator(O::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject}, Id_cache) where {T} = O

@doc raw"""
    destroy(N::Int)

Bosonic annihilation operator with Hilbert space cutoff `N`. This operator
acts on a fock state as ``\hat{a} \ket{n} = \sqrt{n} \ket{n-1}``.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> fock(20, 3)' * a * fock(20, 4)
2.0 + 0.0im
```
"""
destroy(N::Int) = QuantumObject(spdiagm(1 => Array{ComplexF64}(sqrt.(1:N-1))), Operator, [N])

@doc raw"""
    create(N::Int)

Bosonic creation operator with Hilbert space cutoff `N`. This operator
acts on a fock state as ``\hat{a}^\dagger \ket{n} = \sqrt{n+1} \ket{n+1}``.

# Examples

```
julia> a_d = create(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠈⠢⡀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠈⠢⡀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⠢⡀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀

julia> fock(20, 4)' * a_d * fock(20, 3)
2.0 + 0.0im
```
"""
create(N::Int) = QuantumObject(spdiagm(-1 => Array{ComplexF64}(sqrt.(1:N-1))), Operator, [N])

@doc raw"""
    sigmap()

Pauli ladder operator ``\hat{\sigma}_+ = \hat{\sigma}_x + i \hat{\sigma}_y``.
"""
sigmap() = destroy(2)

@doc raw"""
    sigmam()

Pauli ladder operator ``\hat{\sigma}_- = \hat{\sigma}_x - i \hat{\sigma}_y``.
"""
sigmam() = create(2)

@doc raw"""
    sigmax()

Pauli operator ``\hat{\sigma}_x = \hat{\sigma}_- + \hat{\sigma}_+``.
"""
sigmax() = sigmam() + sigmap()

@doc raw"""
    sigmay()

Pauli operator ``\hat{\sigma}_y = i \left( \hat{\sigma}_- - \hat{\sigma}_+ \right)``.
"""
sigmay() = 1im * (sigmam() - sigmap())

@doc raw"""
    sigmaz()

Pauli operator ``\hat{\sigma}_z = \comm{\hat{\sigma}_+}{\hat{\sigma}_-}``.
"""
sigmaz() = sigmap() * sigmam() - sigmam() * sigmap()

@doc raw"""
    eye(N::Int; type=OperatorQuantumObject, dims=[N])

Identity operator ``\hat{\mathbb{1}}`` with Hilbert dimension `N`.
"""
eye(
    N::Int;
    type::ObjType = Operator,
    dims::Vector{Int} = [N],
) where {ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Diagonal(ones(ComplexF64, N)), type, dims)

@doc raw"""
    qeye(N::Int; type=OperatorQuantumObject, dims=[N])

Identity operator ``\hat{\mathbb{1}}`` with Hilbert dimension `N`.
"""
qeye(
    N::Int;
    type::ObjType = Operator,
    dims::Vector{Int} = [N],
) where {ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = eye(N, type = type, dims = dims)

@doc raw"""
    fock(N::Int, pos::Int; dims::Vector{Int}=[N], sparse::Bool=false)

Generates a fock state ``\ket{\psi}`` of dimension `N`. It is also possible
to specify the list of dimensions `dims` if different subsystems are present.
"""
function fock(N::Int, pos::Int; dims::Vector{Int} = [N], sparse::Bool = false)
    if sparse
        return QuantumObject(sparsevec([pos + 1], [1.0 + 0im], N), Ket, dims)
    else
        array = zeros(ComplexF64, N)
        array[pos+1] = 1
        return QuantumObject(array, Ket, dims)
    end
end

"""
    basis(N::Int, pos::Int; dims::Vector{Int}=[N])

Generates a fock state like [`fock`](@ref).
"""
basis(N::Int, pos::Int; dims::Vector{Int} = [N]) = fock(N, pos, dims = dims)

@doc raw"""
    coherent(N::Real, α::T)

Generates a coherent state ``\ket{\alpha}``, which is defined as an eigenvector of the
bosonic annihilation operator ``\hat{a} \ket{\alpha} = \alpha \ket{\alpha}``.
"""
function coherent(N::Real, α::T) where {T<:Number}
    a = destroy(N)
    return exp(α * a' - α' * a) * fock(N, 0)
end

@doc raw"""
    rand_dm(N::Integer; kwargs...)

Generates a random density matrix ``\hat{\rho}``, with the property to be positive definite,
and that ``\Tr \left[ \hat{\rho} \right] = 1``.
"""
function rand_dm(N::Integer; kwargs...)
    ρ = rand(ComplexF64, N, N)
    ρ *= ρ'
    ρ /= tr(ρ)
    return QuantumObject(ρ; kwargs...)
end

@doc raw"""
    projection(N::Int, i::Int, j::Int)

Generates the projection operator ``\hat{O} = \dyad{i}{j}`` with Hilbert space dimension `N`.
"""
projection(N::Int, i::Int, j::Int) = QuantumObject(sparse([i + 1], [j + 1], [1.0 + 0.0im], N, N))

@doc raw"""
    sinm(O::QuantumObject)

Generates the sine of the operator `O`, defined as

``\sin \left( \hat{O} \right) = \frac{e^{i \hat{O}} - e^{-i \hat{O}}}{2 i}``
"""
sinm(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = -0.5im * (exp(1im * O) - exp(-1im * O))

@doc raw"""
    cosm(O::QuantumObject)

Generates the cosine of the operator `O`, defined as

``\cos \left( \hat{O} \right) = \frac{e^{i \hat{O}} + e^{-i \hat{O}}}{2}``
"""
cosm(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = 0.5 * (exp(1im * O) + exp(-1im * O))
