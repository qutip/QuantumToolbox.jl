export Lattice, MultiSiteOperator, DissipativeIsing

@doc raw"""
    Lattice

A Julia constructor for a lattice object. The lattice object is used to define the geometry of the lattice. `Nx` and `Ny` are the number of sites in the x and y directions, respectively. `N` is the total number of sites. `lin_idx` is a `LinearIndices` object and `car_idx` is a `CartesianIndices` object, and they are used to efficiently select sites on the lattice.
"""
Base.@kwdef struct Lattice{TN<:Integer,TLI<:LinearIndices,TCI<:CartesianIndices}
    Nx::TN
    Ny::TN
    N::TN = Nx * Ny
    lin_idx::TLI = LinearIndices((Nx, Ny))
    car_idx::TCI = CartesianIndices((Nx, Ny))
end

#Definition of many-body operators
@doc raw"""
    MultiSiteOperator(dims::Union{AbstractVector, Tuple}, pairs::Pair{<:Integer,<:QuantumObject}...)

A Julia function for generating a multi-site operator ``\\hat{O} = \\hat{O}_i \\hat{O}_j \\cdots \\hat{O}_k``. The function takes a vector of dimensions `dims` and a list of pairs `pairs` where the first element of the pair is the site index and the second element is the operator acting on that site.

# Arguments
- `dims::Union{AbstractVector, Tuple}`: A vector of dimensions of the lattice.
- `pairs::Pair{<:Integer,<:QuantumObject}...`: A list of pairs where the first element of the pair is the site index and the second element is the operator acting on that site.
    
# Returns
`QuantumObject`: A `QuantumObject` representing the multi-site operator.

# Example
```jldoctest
julia> op = MultiSiteOperator(Val(8), 5=>sigmax(), 7=>sigmaz());

julia> op.dims
8-element SVector{8, Int64} with indices SOneTo(8):
 2
 2
 2
 2
 2
 2
 2
 2
```
"""
function MultiSiteOperator(dims::Union{AbstractVector,Tuple}, pairs::Pair{<:Integer,<:QuantumObject}...)
    sites_unsorted = collect(first.(pairs))
    idxs = sortperm(sites_unsorted)
    _sites = sites_unsorted[idxs]
    _ops = collect(last.(pairs))[idxs]
    _dims = collect(dims) # Use this instead of a Tuple, to avoid type instability when indexing on a slice

    sites, ops = _get_unique_sites_ops(_sites, _ops)

    _dims[sites] == [op.to[1].size for op in ops] || throw(ArgumentError("The dimensions of the operators do not match the dimensions of the lattice."))

    data = kron(I(prod(_dims[1:sites[1]-1])), ops[1].data)
    for i in 2:length(sites)
        data = kron(data, I(prod(_dims[sites[i-1]+1:sites[i]-1])), ops[i].data)
    end
    data = kron(data, I(prod(_dims[sites[end]+1:end])))

    return QuantumObject(data; type = Operator, dims = dims)
end
function MultiSiteOperator(N::Union{Integer,Val}, pairs::Pair{<:Integer,<:QuantumObject}...)
    dims = ntuple(j -> 2, makeVal(N))

    return MultiSiteOperator(dims, pairs...)
end
function MultiSiteOperator(latt::Lattice, pairs::Pair{<:Integer,<:QuantumObject}...)
    return MultiSiteOperator(makeVal(latt.N), pairs...)
end

#Definition of nearest-neighbour sites on lattice
periodic_boundary_conditions(i::Integer, N::Integer) = 1 + (i - 1 + N) % N
open_boundary_conditions(i::Integer, N::Integer) = (i >= 1 && i <= N)
periodic_boundary_conditions(i::Vector{Int}, N::Integer) = periodic_boundary_conditions.(i, N)
open_boundary_conditions(i::Vector{Int}, N::Integer) = filter(x -> open_boundary_conditions(x, N), i)

function nearest_neighbor(i::CartesianIndex, latt::Lattice, ::Val{:periodic_bc}; order::Integer = 1)
    row = periodic_boundary_conditions([i[1] + order, i[1] - order], latt.Nx)
    col = periodic_boundary_conditions([i[2] + order, i[2] - order], latt.Ny)
    return vcat([CartesianIndex(r, i[2]) for r in row], [CartesianIndex(i[1], c) for c in col])
end

function nearest_neighbor(i::CartesianIndex, latt::Lattice, ::Val{:open_bc}; order::Integer = 1)
    row = periodic_boundary_conditions([i[1] + order, i[1] - order], latt.Nx)
    col = periodic_boundary_conditions([i[2] + order, i[2] - order], latt.Ny)
    return vcat([CartesianIndex(r, i[2]) for r in row], [CartesianIndex(i[1], c) for c in col])
end

@doc """
    DissipativeIsing(Jx::Real, Jy::Real, Jz::Real, hx::Real, hy::Real, hz::Real, γ::Real, latt::Lattice; boundary_condition::Union{Symbol, Val} = Val(:periodic_bc), order::Integer = 1)

A Julia constructor for a dissipative Ising model. The function returns the Hamiltonian

```math
\\hat{H} = \\frac{J_x}{2} \\sum_{\\langle i, j \\rangle} \\hat{\\sigma}_i^x \\hat{\\sigma}_j^x + \\frac{J_y}{2} \\sum_{\\langle i, j \\rangle} \\hat{\\sigma}_i^y \\hat{\\sigma}_j^y + \\frac{J_z}{2} \\sum_{\\langle i, j \\rangle} \\hat{\\sigma}_i^z \\hat{\\sigma}_j^z + h_x \\sum_i \\hat{\\sigma}_i^x
```

and the collapse operators

```math
\\hat{c}_i = \\sqrt{\\gamma} \\hat{\\sigma}_i^-
```

# Arguments
- `Jx::Real`: The coupling constant in the x-direction.
- `Jy::Real`: The coupling constant in the y-direction.
- `Jz::Real`: The coupling constant in the z-direction.
- `hx::Real`: The magnetic field in the x-direction.
- `hy::Real`: The magnetic field in the y-direction.
- `hz::Real`: The magnetic field in the z-direction.
- `γ::Real`: The local dissipation rate.
- `latt::Lattice`: A [`Lattice`](@ref) object that defines the geometry of the lattice.
- `boundary_condition::Union{Symbol, Val}`: The boundary conditions of the lattice. The possible inputs are `periodic_bc` and `open_bc`, for periodic or open boundary conditions, respectively. The default value is `Val(:periodic_bc)`.
- `order::Integer`: The order of the nearest-neighbour sites. The default value is 1.
"""
function DissipativeIsing(
    Jx::Real,
    Jy::Real,
    Jz::Real,
    hx::Real,
    hy::Real,
    hz::Real,
    γ::Real,
    latt::Lattice;
    boundary_condition::Union{Symbol,Val} = Val(:periodic_bc),
    order::Integer = 1,
)
    S = [MultiSiteOperator(latt, i => sigmam()) for i in 1:latt.N]
    c_ops = sqrt(γ) .* S

    op_sum(S, i::CartesianIndex) =
        S[latt.lin_idx[i]] * sum(S[latt.lin_idx[nearest_neighbor(i, latt, makeVal(boundary_condition); order = order)]])

    H = 0
    if (Jx != 0 || hx != 0)
        S = [MultiSiteOperator(latt, i => sigmax()) for i in 1:latt.N]
        H += Jx / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx) #/2 because we are double counting
        H += hx * sum(S)
    end
    if (Jy != 0 || hy != 0)
        S = [MultiSiteOperator(latt, i => sigmay()) for i in 1:latt.N]
        H += Jy / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
        H += hy * sum(S)
    end
    if (Jz != 0 || hz != 0)
        S = [MultiSiteOperator(latt, i => sigmaz()) for i in 1:latt.N]
        H += Jz / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
        H += hz * sum(S)
    end
    return H, c_ops
end

function _get_unique_sites_ops(sites, ops)
    unique_sites = unique(sites)
    unique_ops = map(i -> prod(ops[findall(==(i), sites)]), unique_sites)

    return unique_sites, unique_ops
end
