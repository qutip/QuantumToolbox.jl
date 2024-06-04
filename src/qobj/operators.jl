#=
Functions for generating (common) quantum operators.
=#

export sigmam, sigmap, sigmax, sigmay, sigmaz
export destroy, create, eye, qeye, projection
export commutator
export spre, spost, sprepost, lindblad_dissipator

@doc raw"""
    commutator(A::QuantumObject, B::QuantumObject; anti::Bool=false)

Return the commutator (or `anti`-commutator) of the two [`QuantumObject`](@ref):
- commutator (`anti=false`): ``AB-BA``
- anticommutator (`anti=true`): ``AB+BA``

Note that `A` and `B` must be [`Operator`](@ref)
"""
commutator(A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject}; anti::Bool=false) = A * B - (-1)^anti * B * A

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
    projection(N::Int, i::Int, j::Int)

Generates the projection operator ``\hat{O} = \dyad{i}{j}`` with Hilbert space dimension `N`.
"""
projection(N::Int, i::Int, j::Int) = QuantumObject(sparse([i + 1], [j + 1], [1.0 + 0.0im], N, N))
