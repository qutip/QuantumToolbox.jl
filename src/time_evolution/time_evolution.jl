abstract type LiouvillianSolver end
struct LiouvillianDirectSolver{T<:Real} <: LiouvillianSolver 
    tol::T
end

abstract type SteadyStateSolver end
abstract type SteadyStateDirectSolver <: SteadyStateSolver end

struct TimeEvolutionSol{TT<:Vector{<:Real}, TS<:AbstractVector, TE<:Matrix{ComplexF64}}
    times::TT
    states::TS
    expect::TE
end

struct TimeEvolutionMCSol{TT<:Vector{<:Vector{<:Real}}, TS<:AbstractVector, TE<:Matrix{ComplexF64}, 
                TEA<:Array{ComplexF64, 3}, TJT<:Vector{<:Vector{<:Real}}, TJW<:Vector{<:Vector{<:Integer}}}
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    jump_times::TJT
    jump_which::TJW
end

LiouvillianDirectSolver(;tol=1e-16) = LiouvillianDirectSolver(tol)

# It is needed to keep track of the index for saving the expectation values
mutable struct ODEProgress{T<:Integer}
    counter::T
end

function next!(p::ODEProgress)
    p.counter += 1
end

abstract type LindbladJumpCallbackType end

struct ContinuousLindbladJumpCallback <: LindbladJumpCallbackType
    interp_points::Int
end

struct DiscreteLindbladJumpCallback <: LindbladJumpCallbackType
end

ContinuousLindbladJumpCallback(;interp_points::Int=10) = ContinuousLindbladJumpCallback(interp_points)

#######################################
    

### LIOUVILLIAN AND STEADYSTATE ###
@doc raw"""
    liouvillian(H::QuantumObject, c_ops::AbstractVector, Id_cache=I(prod(H.dims))

Construct the Liouvillian superoperator for a system Hamiltonian and a set of collapse operators:
``\mathcal{L} \cdot = -i[\hat{H}, \cdot] + \sum_i \mathcal{D}[\hat{O}_i] \cdot``,
where ``\mathcal{D}[\hat{O}_i] \cdot = \hat{O}_i \cdot \hat{O}_i^\dagger - \frac{1}{2} \hat{O}_i^\dagger \hat{O}_i \cdot - \frac{1}{2} \cdot \hat{O}_i^\dagger \hat{O}_i``.

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
function liouvillian(H::QuantumObject{MT,OpType},
    c_ops::Vector{QuantumObject{MT,OperatorQuantumObject}}=Vector{QuantumObject{MT,OperatorQuantumObject}}([]),
    Id_cache=I(prod(H.dims))) where {MT<:AbstractMatrix,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}

    L = isoper(H) ? -1im * (spre(H, Id_cache) - spost(H, Id_cache)) : H
    for c_op in c_ops
        isoper(c_op) ? L += lindblad_dissipator(c_op, Id_cache) : L += c_op
    end
    L
end


# liouvillian(H::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = isoper(H) ? -1im * (spre(H) - spost(H)) : H

function liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real; n_max::Int=4, solver::LSolver=LiouvillianDirectSolver()) where {T1,T2,T3,LSolver<:LiouvillianSolver}

    ((L₀.dims == Lₚ.dims) && (L₀.dims == Lₘ.dims)) || throw(ErrorException("The operators are not of the same Hilbert dimension."))

    _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, solver, n_max=n_max)
end

function liouvillian_floquet(H::QuantumObject{<:AbstractArray{T1},OpType1},
    c_ops::AbstractVector,
    Hₚ::QuantumObject{<:AbstractArray{T2},OpType2},
    Hₘ::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, solver::LSolver=LiouvillianDirectSolver()) where {T1,T2,T3,
                                                                            OpType1<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            OpType2<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            OpType3<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            LSolver<:LiouvillianSolver}

    liouvillian_floquet(liouvillian(H, c_ops), liouvillian(Hₚ), liouvillian(Hₘ), ω, solver=solver, n_max=n_max)
end

@doc raw"""
    liouvillian_generalized(H::QuantumObject, fields::Vector, 
    T_list::Vector; N_trunc::Int=size(H,1), tol::Float64=0.0, σ_filter::Union{Nothing, Real}=nothing)

Constructs the generalized Liouvillian for a system coupled to a bath of harmonic oscillators.

See, e.g., Settineri, Alessio, et al. "Dissipation and thermal noise in hybrid quantum systems in the ultrastrong-coupling regime." Physical Review A 98.5 (2018): 053834.
"""
function liouvillian_generalized(H::QuantumObject{MT, OperatorQuantumObject}, fields::Vector, 
    T_list::Vector{<:Real}; N_trunc::Int=size(H,1), tol::Real=1e-12, σ_filter::Union{Nothing, Real}=nothing) where {MT<:AbstractMatrix}

    (length(fields) == length(T_list)) || throw(DimensionMismatch("The number of fields, ωs and Ts must be the same."))

    dims = N_trunc == size(H,1) ? H.dims : [N_trunc]
    E2, U2 = eigen(H)
    E = real.(E2[1:N_trunc])
    U = QuantumObject(U2, dims=H.dims)

    H_d = QuantumObject(spdiagm(complex(E)), dims=dims)

    Ω = E' .- E
    Ωp = triu(dense_to_sparse(Ω, tol), 1)

    # Filter in the Hilbert space
    σ = isnothing(σ_filter) ? 500 * maximum([norm(field) / length(field) for field in fields]) : σ_filter
    F1 = QuantumObject(gaussian.(Ω, 0, σ), dims=dims)
    F1 = dense_to_sparse(F1, tol)
    
    # Filter in the Liouville space
    M1 = ones(N_trunc, N_trunc)
    Ω1 = kron(Ω, M1)
    Ω2 = kron(M1, Ω)
    Ωdiff = Ω1 .- Ω2
    F2 = QuantumObject(gaussian.(Ωdiff, 0, σ), SuperOperatorQuantumObject, dims)
    F2 = dense_to_sparse(F2, tol)

    L = liouvillian(H_d)

    for i in eachindex(fields)
        # The operator that couples the system to the bath in the eigenbasis
        X_op = dense_to_sparse((U' * fields[i] * U).data[1:N_trunc, 1:N_trunc], tol)
        if ishermitian(fields[i])
            X_op = (X_op + X_op') / 2 # Make sure it's hermitian
        end

        # Ohmic reservoir
        N_th = n_th.(Ωp, T_list[i])
        Sp₀ = QuantumObject( triu(X_op, 1), dims=dims )
        Sp₁ = QuantumObject( droptol!( (@. Ωp * N_th * Sp₀.data), tol), dims=dims )
        Sp₂ = QuantumObject( droptol!( (@. Ωp * (1 + N_th) * Sp₀.data), tol), dims=dims )
        # S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

        L += 1 / 2 * ( F2 .* (sprepost(Sp₁', Sp₀) + sprepost(Sp₀', Sp₁)) - spre(F1 .* (Sp₀ * Sp₁')) - spost(F1 .* (Sp₁ * Sp₀')) )
        L += 1 / 2 * ( F2 .* (sprepost(Sp₂, Sp₀') + sprepost(Sp₀, Sp₂')) - spre(F1 .* (Sp₀' * Sp₂)) - spost(F1 .* (Sp₂' * Sp₀)) )
    end

    return E, U, L
end

function _liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real, solver::LiouvillianDirectSolver; n_max::Int=4) where {T1,T2,T3}

    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data
    L_p_dense = sparse_to_dense(Lₚ.data)
    L_m_dense = sparse_to_dense(Lₘ.data)

    S = -(L_0 - 1im * n_max * ω * I) \ L_p_dense
    T = -(L_0 + 1im * n_max * ω * I) \ L_m_dense

    for n_i in n_max-1:-1:1
        S = -(L_0 - 1im * n_i * ω * I + L_m * S) \ L_p_dense
        T = -(L_0 + 1im * n_i * ω * I + L_p * T) \ L_m_dense
    end

    solver.tol == 0 && return QuantumObject(L_0 + L_m * S + L_p * T, SuperOperatorQuantumObject, L₀.dims)
    return QuantumObject(dense_to_sparse(L_0 + L_m * S + L_p * T, solver.tol), SuperOperatorQuantumObject, L₀.dims)
end

function steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject};
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    _steadystate(L, solver)
end

function steadystate(H::QuantumObject{<:AbstractArray{T},OpType}, c_ops::Vector,
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},SSSolver<:SteadyStateSolver}

    L = liouvillian(H, c_ops)
    steadystate(L, solver=solver)
end

function _steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::Type{SteadyStateDirectSolver}) where {T}

    L_tmp = copy(L.data)
    N = prod(L.dims)
    weight = norm(L_tmp, 1) / length(L_tmp)
    v0 = zeros(ComplexF64, N^2) # This is not scalable for GPU arrays
    v0[1] = weight

    L_tmp[1, [N * (i - 1) + i for i in 1:N]] .+= weight

    ρss_vec = L_tmp \ v0
    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    QuantumObject(ρss, OperatorQuantumObject, L.dims)
end

@doc raw"""
    steadystate_floquet(H_0::QuantumObject,
        c_ops::Vector, H_p::QuantumObject,
        H_m::QuantumObject,
        ω::Real; n_max::Int=4, lf_solver::LSolver=LiouvillianDirectSolver(),
        ss_solver::Type{SSSolver}=SteadyStateDirectSolver)

Calculates the steady state of a periodically driven system.
Here `H_0` is the Hamiltonian or the Liouvillian of the undriven system.
Considering a monochromatic drive at frequency ``\\omega``, we divide it into two parts,
`H_p` and `H_m`, where `H_p` oscillates
as ``e^{i \\omega t}`` and `H_m` oscillates as ``e^{-i \\omega t}``.
`n_max` is the number of iterations used to obtain the effective Liouvillian,
`lf_solver` is the solver used to solve the effective Liouvillian,
and `ss_solver` is the solver used to solve the steady state.
"""
function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T1},OpType1},
    c_ops::AbstractVector, H_p::QuantumObject{<:AbstractArray{T2},OpType2},
    H_m::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, lf_solver::LSolver=LiouvillianDirectSolver(),
    ss_solver::Type{SSSolver}=SteadyStateDirectSolver) where {T1,T2,T3,OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    LSolver<:LiouvillianSolver,SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end

function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T1},OpType1},
    H_p::QuantumObject{<:AbstractArray{T2},OpType2},
    H_m::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, lf_solver::LSolver=LiouvillianDirectSolver(),
    ss_solver::Type{SSSolver}=SteadyStateDirectSolver) where {T1,T2,T3,OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    LSolver<:LiouvillianSolver,SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end