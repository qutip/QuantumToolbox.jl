#=
Functions for generating (common) quantum states.
=#

export zero_ket, fock, coherent, rand_ket
export fock_dm, coherent_dm, thermal_dm, maximally_mixed_dm, rand_dm
export spin_state, spin_coherent
export bell_state, singlet_state, triplet_states, w_state, ghz_state

@doc raw"""
    zero_ket(dimensions)

Returns a zero [`Ket`](@ref) vector with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{Dimensions,AbstractVector{Int}, Tuple}`: list of dimensions representing the each number of basis in the subsystems.

!!! warning "Beware of type-stability!"
    It is highly recommended to use `zero_ket(dimensions)` with `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
zero_ket(dimensions::Int) = QuantumObject(zeros(ComplexF64, dimensions), Ket(), dimensions)
zero_ket(dimensions::Union{Dimensions, AbstractVector{Int}, Tuple}) =
    QuantumObject(zeros(ComplexF64, prod(dimensions)), Ket(), dimensions)

@doc raw"""
    fock(N::Int, j::Int=0; dims::Union{Int,AbstractVector{Int},Tuple}=N, sparse::Union{Bool,Val}=Val(false))
    basis(N::Int, j::Int=0; dims::Union{Int,AbstractVector{Int},Tuple}=N, sparse::Union{Bool,Val}=Val(false))

Generates a fock state ``\ket{\psi}`` of dimension `N`. 

It is also possible to specify the list of dimensions `dims` if different subsystems are present.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `fock(N, j, dims=dims, sparse=Val(sparse))` instead of `fock(N, j, dims=dims, sparse=sparse)`. Consider also to use `dims` as a `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) instead of `Vector`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.

!!! note
    `basis(N, j; dims = dims, sparse = sparse)` is a synonym of `fock(N, j; dims = dims, sparse = sparse)`.
"""
function fock(N::Int, j::Int = 0; dims::Union{Int, AbstractVector{Int}, Tuple} = N, sparse::Union{Bool, Val} = Val(false))
    (0 <= j < N) || throw(ArgumentError("Invalid argument j, must satisfy: 0 ≤ j ≤ N-1"))
    if getVal(sparse)
        array = sparsevec([j + 1], [1.0 + 0im], N)
    else
        array = [i == (j + 1) ? 1.0 + 0im : 0.0 + 0im for i in 1:N]
    end
    return QuantumObject(array; type = Ket(), dims = dims)
end

@doc raw"""
    coherent(N::Int, α::Number)

Generates a [coherent state](https://en.wikipedia.org/wiki/Coherent_state) ``|\alpha\rangle``, which is defined as an eigenvector of the bosonic annihilation operator ``\hat{a} |\alpha\rangle = \alpha |\alpha\rangle``.

This state is constructed via the displacement operator [`displace`](@ref) and zero-fock state [`fock`](@ref): ``|\alpha\rangle = \hat{D}(\alpha) |0\rangle``
"""
coherent(N::Int, α::T) where {T <: Number} = displace(N, α) * fock(N, 0)

@doc raw"""
    rand_ket(dimensions)

Generate a random normalized [`Ket`](@ref) vector with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{Dimensions,AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `rand_ket(dimensions)` with `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
rand_ket(dimensions::Int) = rand_ket(SVector(dimensions))
function rand_ket(dimensions::Union{Dimensions, AbstractVector{Int}, Tuple})
    N = prod(dimensions)
    ψ = rand(ComplexF64, N) .- (0.5 + 0.5im)
    return QuantumObject(normalize!(ψ); type = Ket(), dims = dimensions)
end

@doc raw"""
    fock_dm(N::Int, j::Int=0; dims::Union{Int,AbstractVector{Int},Tuple}=N, sparse::Union{Bool,Val}=Val(false))

Density matrix representation of a Fock state.

Constructed via outer product of [`fock`](@ref).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `fock_dm(N, j, dims=dims, sparse=Val(sparse))` instead of `fock_dm(N, j, dims=dims, sparse=sparse)`. Consider also to use `dims` as a `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) instead of `Vector`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function fock_dm(
        N::Int,
        j::Int = 0;
        dims::Union{Int, AbstractVector{Int}, Tuple} = N,
        sparse::Union{Bool, Val} = Val(false),
    )
    ψ = fock(N, j; dims = dims, sparse = sparse)
    return ket2dm(ψ)
end

@doc raw"""
    coherent_dm(N::Int, α::Number)

Density matrix representation of a [coherent state](https://en.wikipedia.org/wiki/Coherent_state).

Constructed via outer product of [`coherent`](@ref).
"""
function coherent_dm(N::Int, α::T) where {T <: Number}
    ψ = coherent(N, α)
    return ket2dm(ψ)
end

@doc raw"""
    thermal_dm(N::Int, n::Real; sparse::Union{Bool,Val}=Val(false))

Density matrix for a thermal state (generating thermal state probabilities) with the following arguments:
- `N::Int`: Number of basis states in the Hilbert space
- `n::Real`: Expectation value for number of particles in the thermal state.
- `sparse::Union{Bool,Val}`: If `true`, return a sparse matrix representation.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `thermal_dm(N, n, sparse=Val(sparse))` instead of `thermal_dm(N, n, sparse=sparse)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function thermal_dm(N::Int, n::T; sparse::Union{Bool, Val} = Val(false)) where {T <: Real}
    β = log(1 + 1 / n)
    P = _complex_float_type(T)[_Boltzmann_weight(β, j) for j in 0:(N - 1)]
    P /= sum(P)
    if getVal(sparse)
        return QuantumObject(spdiagm(0 => P), Operator(), N)
    else
        return QuantumObject(diagm(0 => P), Operator(), N)
    end
end

@doc raw"""
    maximally_mixed_dm(dimensions)

Returns the maximally mixed density matrix with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{Dimensions,AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `maximally_mixed_dm(dimensions)` with `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
maximally_mixed_dm(dimensions::Int) =
    QuantumObject(diagm(0 => fill(ComplexF64(1 / dimensions), dimensions)), Operator(), SVector(dimensions))
function maximally_mixed_dm(dimensions::Union{Dimensions, AbstractVector{Int}, Tuple})
    N = prod(dimensions)
    return QuantumObject(diagm(0 => fill(ComplexF64(1 / N), N)), Operator(), dimensions)
end

@doc raw"""
    rand_dm(dimensions; rank::Int=prod(dimensions))

Generate a random density matrix from Ginibre ensemble with given argument `dimensions` and `rank`, ensuring that it is positive semi-definite and trace equals to `1`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{Dimensions,AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

The default keyword argument `rank = prod(dimensions)` (full rank).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `rand_dm(dimensions; rank=rank)` with `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) instead of `Vector`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.

# References
- [J. Ginibre, Statistical ensembles of complex, quaternion, and real matrices, Journal of Mathematical Physics 6.3 (1965): 440-449](https://doi.org/10.1063/1.1704292)
- [K. Życzkowski, et al., Generating random density matrices, Journal of Mathematical Physics 52, 062201 (2011)](http://dx.doi.org/10.1063/1.3595693)
"""
rand_dm(dimensions::Int; rank::Int = prod(dimensions)) = rand_dm(SVector(dimensions), rank = rank)
function rand_dm(dimensions::Union{Dimensions, AbstractVector{Int}, Tuple}; rank::Int = prod(dimensions))
    N = prod(dimensions)
    (rank < 1) && throw(DomainError(rank, "The argument rank must be larger than 1."))
    (rank > N) && throw(DomainError(rank, "The argument rank cannot exceed dimensions."))

    X = _Ginibre_ensemble(N, rank)
    ρ = X * X'
    ρ /= tr(ρ)
    return QuantumObject(ρ; type = Operator(), dims = dimensions)
end

@doc raw"""
    spin_state(j::Real, m::Real)

Generate the spin state: ``|j, m\rangle``

The eigenstate of the Spin-`j` ``\hat{S}_z`` operator with eigenvalue `m`, where where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
function spin_state(j::Real, m::Real)
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    Δ = j - m
    ((floor(Δ) != Δ) || (Δ < 0)) &&
        throw(ArgumentError("Invalid eigenvalue m: (j - m) must be a non-negative integer."))
    (m < (-j)) && throw(ArgumentError("Invalid eigenvalue m, must satisfy: -j ≤ m ≤ j"))

    return fock(Int(J), Int(Δ))
end

@doc raw"""
    spin_coherent(j::Real, θ::Real, ϕ::Real)

Generate the coherent spin state (rotation of the ``|j, j\rangle`` state), namely

```math
|\theta, \phi \rangle = \hat{R}(\theta, \phi) |j, j\rangle
```

where the rotation operator is defined as

```math
\hat{R}(\theta, \phi) = \exp \left( \frac{\theta}{2} (\hat{S}_- e^{i\phi} - \hat{S}_+ e^{-i\phi}) \right)
```

and ``\hat{S}_\pm`` are plus and minus Spin-`j` operators, respectively.

# Arguments
- `j::Real`: The spin quantum number and can be a non-negative integer or half-integer
- `θ::Real`: rotation angle from z-axis
- `ϕ::Real`: rotation angle from x-axis

See also [`jmat`](@ref) and [`spin_state`](@ref).

# Reference
- [Robert Jones, Spin Coherent States and Statistical Physics](https://web.mit.edu/8.334/www/grades/projects/projects19/JonesRobert.pdf)
"""
function spin_coherent(j::Real, θ::Real, ϕ::Real)
    Sm = jmat(j, Val(:-))
    return exp(0.5 * θ * (Sm * exp(1im * ϕ) - Sm' * exp(-1im * ϕ))) * spin_state(j, j)
end

@doc raw"""
    bell_state(x::Union{Int}, z::Union{Int})

Return the [Bell state](https://en.wikipedia.org/wiki/Bell_state) depending on the arguments `(x, z)`:
- `(0, 0)`: ``| \Phi^+ \rangle = ( |00\rangle + |11\rangle ) / \sqrt{2}``
- `(0, 1)`: ``| \Phi^- \rangle = ( |00\rangle - |11\rangle ) / \sqrt{2}``
- `(1, 0)`: ``| \Psi^+ \rangle = ( |01\rangle + |10\rangle ) / \sqrt{2}``
- `(1, 1)`: ``| \Psi^- \rangle = ( |01\rangle - |10\rangle ) / \sqrt{2}``

Here, `x = 1` (`z = 1`) means applying Pauli-``X`` ( Pauli-``Z``) unitary transformation on ``| \Phi^+ \rangle``.

# Example

```jldoctest
julia> bell_state(0, 0)

Quantum Object:   type=Ket()   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> bell_state(Val(1), Val(0))

Quantum Object:   type=Ket()   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
                0.0 + 0.0im
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
```

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `bell_state(Val(x), Val(z))` instead of `bell_state(x, z)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
bell_state(x::Int, z::Int) = bell_state(Val(x), Val(z))
bell_state(::Val{0}, ::Val{0}) = QuantumObject(ComplexF64[1, 0, 0, 1] / sqrt(2), Ket(), (2, 2))
bell_state(::Val{0}, ::Val{1}) = QuantumObject(ComplexF64[1, 0, 0, -1] / sqrt(2), Ket(), (2, 2))
bell_state(::Val{1}, ::Val{0}) = QuantumObject(ComplexF64[0, 1, 1, 0] / sqrt(2), Ket(), (2, 2))
bell_state(::Val{1}, ::Val{1}) = QuantumObject(ComplexF64[0, 1, -1, 0] / sqrt(2), Ket(), (2, 2))
bell_state(::Val{T1}, ::Val{T2}) where {T1, T2} = throw(ArgumentError("Invalid Bell state: $(T1), $(T2)"))

@doc raw"""
    singlet_state()

Return the two particle singlet state: ``\frac{1}{\sqrt{2}} ( |01\rangle - |10\rangle )``
"""
singlet_state() = QuantumObject(ComplexF64[0, 1, -1, 0] / sqrt(2), Ket(), (2, 2))

@doc raw"""
    triplet_states()

Return a list of the two particle triplet states: 

- ``|11\rangle``
- ``( |01\rangle + |10\rangle ) / \sqrt{2}``
- ``|00\rangle``
"""
function triplet_states()
    return QuantumObject[
        QuantumObject(ComplexF64[0, 0, 0, 1], Ket(), (2, 2)),
        QuantumObject(ComplexF64[0, 1, 1, 0] / sqrt(2), Ket(), (2, 2)),
        QuantumObject(ComplexF64[1, 0, 0, 0], Ket(), (2, 2)),
    ]
end

@doc raw"""
    w_state(n::Union{Int,Val})

Returns the `n`-qubit [W-state](https://en.wikipedia.org/wiki/W_state):

```math
\frac{1}{\sqrt{n}} \left( |100...0\rangle + |010...0\rangle + \cdots + |00...01\rangle \right)
```

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `w_state(Val(n))` instead of `w_state(n)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
function w_state(::Val{n}) where {n}
    (n >= 2) || throw(ArgumentError("Invalid argument n, must satisfy: n ≥ 2"))

    nzind = 2 .^ (0:(n - 1)) .+ 1
    nzval = fill(ComplexF64(1 / sqrt(n)), n)
    data = zeros(ComplexF64, 2^n)
    @inbounds data[nzind] .= nzval
    return QuantumObject(data, Ket(), ntuple(x -> 2, Val(n)))
end
w_state(n::Int) = w_state(Val(n))

@doc raw"""
    ghz_state(n::Union{Int,Val}; d::Int=2)

Returns the generalized `n`-qudit [Greenberger–Horne–Zeilinger (GHZ) state](https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state):

```math
\frac{1}{\sqrt{d}} \sum_{i=0}^{d-1} | i \rangle \otimes \cdots \otimes | i \rangle
```

Here, `d` specifies the dimension of each qudit. Default to `d=2` (qubit).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `ghz_state(Val(n))` instead of `ghz_state(n)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
function ghz_state(::Val{n}; d::Int = 2) where {n}
    (n >= 2) || throw(ArgumentError("Invalid argument n, must satisfy: n ≥ 2"))
    (d >= 2) || throw(ArgumentError("Invalid argument d, must satisfy: d ≥ 2"))

    nzind = collect((0:(d - 1)) .* Int((d^n - 1) / (d - 1)) .+ 1)
    nzval = fill(ComplexF64(1 / sqrt(d)), d)
    data = zeros(ComplexF64, d^n)
    @inbounds data[nzind] .= nzval
    return QuantumObject(data, Ket(), ntuple(x -> d, Val(n)))
end
ghz_state(n::Int; d::Int = 2) = ghz_state(Val(n), d = d)
