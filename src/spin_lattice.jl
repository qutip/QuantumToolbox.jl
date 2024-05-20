export Lattice, mb, TFIM, nn, sx, sy, sz, sm, sp, pbc, obc

sx = sigmax()
sy = -sigmay()
sz = -sigmaz()
sm = (sx - 1im * sy) / 2
sp = (sx + 1im * sy) / 2

#Lattice structure
Base.@kwdef struct Lattice{TN<:Integer,TLI<:LinearIndices,TCI<:CartesianIndices}
    Nx::TN
    Ny::TN
    N::TN = Nx * Ny
    lin_idx::TLI = LinearIndices((Nx, Ny))
    car_idx::TCI = CartesianIndices((Nx, Ny))
end

#Definition of many-body operators
function mb(
    s::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    i::Integer,
    N::Integer,
) where {T1}
    T = s.dims[1]
    QuantumObject(kron(eye(T^(i - 1)), s, eye(T^(N - i))); dims = fill(2, N))
end
mb(
    s::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    i::Integer,
    latt::Lattice,
) where {T1} = mb(s, i, latt.N)
mb(
    s::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    row::Integer,
    col::Integer,
    latt::Lattice,
) where {T1} = mb(s, latt.idx[row, col], latt.N)
mb(
    s::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    x::CartesianIndex,
    latt::Lattice,
) where {T1} = mb(s, latt.idx[x], latt.N)

#Definition of nearest-neighbour sites on lattice
pbc(i::Integer, N::Integer) = 1 + (i - 1 + N) % N
obc(i::Integer, N::Integer) = (i >= 1 && i <= N)
pbc(i::Vector{Int}, N::Integer) = pbc.(i, N)
obc(i::Vector{Int}, N::Integer) = filter(x -> obc(x, N), i)

function nn(i::CartesianIndex, latt::Lattice, bc::Function; order::Integer = 1)
    row = bc([i[1] + order, i[1] - order], latt.Nx)
    col = bc([i[2] + order, i[2] - order], latt.Ny)
    vcat([CartesianIndex(r, i[2]) for r in row], [CartesianIndex(i[1], c) for c in col])
end

function TFIM(
    Jx::Real,
    Jy::Real,
    Jz::Real,
    hx::Real,
    γ::Real,
    latt::Lattice;
    bc::Function = pbc,
    order::Integer = 1,
)
    S = [mb(sm, i, latt) for i in 1:latt.N]
    c_ops = sqrt(γ) .* S

    op_sum(
        S::Vector{QuantumObject{SparseMatrixCSC{ComplexF64,Int64},OperatorQuantumObject}},
        i::CartesianIndex,
    ) = S[latt.lin_idx[i]] * sum(S[latt.lin_idx[nn(i, latt, bc; order = order)]])

    H = 0
    if (Jx != 0 || hx != 0)
        S .= [mb(sx, i, latt) for i in 1:latt.N]
        H += Jx / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx) #/2 because we are double counting
        H += hx * sum(S)
    end
    if Jy != 0
        S .= [mb(sy, i, latt) for i in 1:latt.N]
        H += Jy / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
    end
    if Jz != 0
        S .= [mb(sz, i, latt) for i in 1:latt.N]
        H += Jz / 2 * mapreduce(i -> op_sum(S, i), +, latt.car_idx)
    end
    H, c_ops
end;
