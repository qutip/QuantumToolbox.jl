export Lattice, SingleSiteOperator, DissipativeIsing

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
    SingleSiteOperator(O::QuantumObject, i::Integer, N::Integer)

A Julia constructor for a single-site operator. `s` is the operator acting on the site. `i` is the site index, and `N` is the total number of sites. The function returns a `QuantumObject` given by ``\\mathbb{1}^{\\otimes (i - 1)} \\otimes \hat{O} \\otimes \\mathbb{1}^{\\otimes (N - i)}``.
"""
function SingleSiteOperator(O::QuantumObject{DT,OperatorQuantumObject}, i::Integer, N::Integer) where {DT}
    T = O.dims[1]
    return QuantumObject(kron(eye(T^(i - 1)), O, eye(T^(N - i))); dims = ntuple(j -> 2, Val(N)))
end
SingleSiteOperator(O::QuantumObject{DT,OperatorQuantumObject}, i::Integer, latt::Lattice) where {DT} =
    SingleSiteOperator(O, i, latt.N)
SingleSiteOperator(O::QuantumObject{DT,OperatorQuantumObject}, row::Integer, col::Integer, latt::Lattice) where {DT} =
    SingleSiteOperator(O, latt.idx[row, col], latt.N)
SingleSiteOperator(O::QuantumObject{DT,OperatorQuantumObject}, x::CartesianIndex, latt::Lattice) where {DT} =
    SingleSiteOperator(O, latt.idx[x], latt.N)

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
    S = [SingleSiteOperator(sigmam(), i, latt) for i in 1:latt.N]
    c_ops = sqrt(γ) .* S

    op_sum(S, i::CartesianIndex) =
        S[latt.lin_idx[i]] * sum(S[latt.lin_idx[nearest_neighbor(i, latt, makeVal(boundary_condition); order = order)]])

    H = 0
    if (Jx != 0 || hx != 0)
        S = [SingleSiteOperator(sigmax(), i, latt) for i in 1:latt.N]
        H += Jx / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx) #/2 because we are double counting
        H += hx * sum(S)
    end
    if (Jy != 0 || hy != 0)
        S = [SingleSiteOperator(sigmay(), i, latt) for i in 1:latt.N]
        H += Jy / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
        H += hy * sum(S)
    end
    if (Jz != 0 || hz != 0)
        S = [SingleSiteOperator(sigmaz(), i, latt) for i in 1:latt.N]
        H += Jz / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
        H += hz * sum(S)
    end
    return H, c_ops
end
