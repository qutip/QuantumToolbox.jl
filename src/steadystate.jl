export steadystate, steadystate_floquet
export SteadyStateSolver,
    SteadyStateLinearSolver,
    SteadyStateEigenSolver,
    SteadyStateDirectSolver,
    SteadyStateFloquetSolver,
    SSFloquetLinearSystem,
    SSFloquetEffectiveLiouvillian

abstract type SteadyStateSolver end
abstract type SteadyStateFloquetSolver end
struct SSFloquetLinearSystem <: SteadyStateFloquetSolver end
Base.@kwdef struct SSFloquetEffectiveLiouvillian{SSST<:SteadyStateSolver} <: SteadyStateFloquetSolver
    steadystate_solver::SSST = SteadyStateDirectSolver()
end

struct SteadyStateDirectSolver <: SteadyStateSolver end
struct SteadyStateEigenSolver <: SteadyStateSolver end
Base.@kwdef struct SteadyStateLinearSolver{MT<:Union{LinearSolve.SciMLLinearSolveAlgorithm,Nothing}} <:
                   SteadyStateSolver
    alg::MT = nothing
    Pl::Union{Function,Nothing} = nothing
    Pr::Union{Function,Nothing} = nothing
end

function steadystate(
    L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject};
    solver::SteadyStateSolver = SteadyStateDirectSolver(),
    kwargs...,
) where {T}
    return _steadystate(L, solver; kwargs...)
end

function steadystate(
    H::QuantumObject{<:AbstractArray{T},OpType},
    c_ops::AbstractVector;
    solver::SteadyStateSolver = SteadyStateDirectSolver(),
    kwargs...,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    L = liouvillian(H, c_ops)

    return steadystate(L; solver = solver, kwargs...)
end

function _steadystate(
    L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::SteadyStateLinearSolver;
    kwargs...,
) where {T}
    L_tmp = L.data
    N = prod(L.dims)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _get_dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _get_dense_similar(L_tmp, N)
    cols = _get_dense_similar(L_tmp, N)
    datas = _get_dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(datas, weight)
    Tn = sparse(rows, cols, datas, N^2, N^2)
    L_tmp = L_tmp + Tn

    (haskey(kwargs, :Pl) || haskey(kwargs, :Pr)) && error("The use of preconditioners must be defined in the solver.")
    if !isnothing(solver.Pl)
        kwargs = merge((; kwargs...), (Pl = solver.Pl(L_tmp),))
    elseif isa(L_tmp, SparseMatrixCSC)
        kwargs = merge((; kwargs...), (Pl = ilu(L_tmp, τ = 0.01),))
    end
    !isnothing(solver.Pr) && (kwargs = merge((; kwargs...), (Pr = solver.Pr(L_tmp),)))

    prob = LinearProblem(L_tmp, v0)
    ρss_vec = solve(prob, solver.alg; kwargs...).u

    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator, L.dims)
end

function _steadystate(
    L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::SteadyStateEigenSolver;
    kwargs...,
) where {T}
    N = prod(L.dims)

    kwargs = merge((sigma = 1e-8, k = 1), (; kwargs...))

    ρss_vec = eigsolve(L; kwargs...).vectors[:, 1]
    ρss = reshape(ρss_vec, N, N)
    ρss /= tr(ρss)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator, L.dims)
end

function _steadystate(
    L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::SteadyStateDirectSolver,
) where {T}
    L_tmp = L.data
    N = prod(L.dims)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _get_dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _get_dense_similar(L_tmp, N)
    cols = _get_dense_similar(L_tmp, N)
    datas = _get_dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(datas, weight)
    Tn = sparse(rows, cols, datas, N^2, N^2)
    L_tmp = L_tmp + Tn

    ρss_vec = L_tmp \ v0 # This is still not supported on GPU, yet
    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator, L.dims)
end

@doc raw"""
    function steadystate_floquet(
        H_0::QuantumObject{MT,OpType1},
        H_p::QuantumObject{<:AbstractArray,OpType2},
        H_m::QuantumObject{<:AbstractArray,OpType3},
        ωd::Number,
        c_ops::AbstractVector = QuantumObject{MT,OpType1}[];
        n_max::Integer = 2,
        tol::R = 1e-8,
        solver::FSolver = SSFloquetLinearSystem,
        kwargs...,
    )

Calculates the steady state of a periodically driven system.
Here `H_0` is the Hamiltonian or the Liouvillian of the undriven system.
Considering a monochromatic drive at frequency ``\\omega_d``, we divide it into two parts,
`H_p` and `H_m`, where `H_p` oscillates
as ``e^{i \\omega t}`` and `H_m` oscillates as ``e^{-i \\omega t}``.
There are two solvers available for this function:
- `SSFloquetLinearSystem`: Solves the linear system of equations.
- `SSFloquetEffectiveLiouvillian`: Solves the effective Liouvillian.
For both cases, `n_max` is the number of Fourier components to consider, and `tol` is the tolerance for the solver.

In the case of `SSFloquetLinearSystem`, the full linear system is solved at once:

``( \mathcal{L}_0 - i n \omega_d ) \hat{\rho}_n + \mathcal{L}_1 \hat{\rho}\_{n-1} + \mathcal{L}\_{-1} \hat{\rho}\_{n+1} = 0``

This is a tridiagonal linear system of the form

``\mathbf{A} \cdot \mathbf{b} = 0``

where

``\mathbf{A} = \begin{pmatrix}
\mathcal{L}\_0 - i (-n\_{\text{max}}) \omega\_{\textrm{d}} & \mathcal{L}\_{-1} & 0 & \cdots & 0 \\
\mathcal{L}_1 & \mathcal{L}_0 - i (-n\_{\text{max}}+1) \omega\_{\textrm{d}} & \mathcal{L}\_{-1} & \cdots & 0 \\
0 & \mathcal{L}\_1 & \mathcal{L}\_0 - i (-n\_{\text{max}}+2) \omega\_{\textrm{d}} & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \mathcal{L}_0 - i n\_{\text{max}} \omega\_{\textrm{d}}
\end{pmatrix}``

and

``\mathbf{b} = \begin{pmatrix}
\hat{\rho}\_{-n_{\text{max}}} \\
\hat{\rho}\_{-n_{\text{max}}+1} \\
\vdots \\
\hat{\rho}\_{0} \\
\vdots \\
\hat{\rho}\_{n_{\text{max}}-1} \\
\hat{\rho}\_{n_{\text{max}}}
\end{pmatrix}``

This will allow to simultaneously obtain all the ``\hat{\rho}_n``.

In the case of `SSFloquetEffectiveLiouvillian`, instead, the effective Liouvillian is calculated using the matrix continued fraction method.

!!! note "different return"
    The two solvers returns different objects. The `SSFloquetLinearSystem` returns a list of `QuantumObject`, containing the density matrices for each Fourier component (`\hat{\rho}_{-n}`, with `n` from `0` to `n_max`), while the `SSFloquetEffectiveLiouvillian` returns only `\hat{\rho}_0`. 

## Arguments
- `H_0::QuantumObject`: The Hamiltonian or the Liouvillian of the undriven system.
- `H_p::QuantumObject`: The Hamiltonian or the Liouvillian of the part of the drive that oscillates as ``e^{i \\omega t}``.
- `H_m::QuantumObject`: The Hamiltonian or the Liouvillian of the part of the drive that oscillates as ``e^{-i \\omega t}``.
- `ωd::Number`: The frequency of the drive.
- `c_ops::AbstractVector = QuantumObject`: The optional collapse operators.
- `n_max::Integer = 2`: The number of Fourier components to consider.
- `tol::R = 1e-8`: The tolerance for the solver.
- `solver::FSolver = SSFloquetLinearSystem`: The solver to use.
- `kwargs...`: Additional keyword arguments to be passed to the solver.
"""
function steadystate_floquet(
    H_0::QuantumObject{MT,OpType1},
    H_p::QuantumObject{<:AbstractArray,OpType2},
    H_m::QuantumObject{<:AbstractArray,OpType3},
    ωd::Number,
    c_ops::AbstractVector = QuantumObject{MT,OpType1}[];
    n_max::Integer = 2,
    tol::R = 1e-8,
    solver::FSolver = SSFloquetLinearSystem(),
    kwargs...,
) where {
    MT<:AbstractArray,
    OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    R<:Real,
    FSolver<:SteadyStateFloquetSolver,
}
    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)
    return _steadystate_floquet(L_0, L_p, L_m, ωd, solver; n_max = n_max, tol = tol, kwargs...)
end

function _steadystate_floquet(
    L_0::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    L_p::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    L_m::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ωd::Number,
    solver::SSFloquetLinearSystem;
    n_max::Integer = 1,
    tol::R = 1e-8,
    kwargs...,
) where {T1,T2,T3,R<:Real}
    T = promote_type(T1, T2, T3)

    L_0_mat = get_data(L_0)
    L_p_mat = get_data(L_p)
    L_m_mat = get_data(L_m)

    N = size(L_0_mat, 1)
    Ns = isqrt(N)
    n_fourier = 2 * n_max + 1
    n_list = -n_max:n_max

    weight = 1
    Mn = sparse(ones(Ns), [Ns * (j - 1) + j for j in 1:Ns], fill(weight, Ns), N, N)
    L = L_0_mat + Mn

    M = spzeros(T, n_fourier * N, n_fourier * N)
    M += kron(spdiagm(1 => ones(n_fourier - 1)), L_m_mat)
    M += kron(spdiagm(-1 => ones(n_fourier - 1)), L_p_mat)
    for i in 1:n_fourier
        n = n_list[i]
        M += kron(sparse([i], [i], one(T), n_fourier, n_fourier), L - 1im * ωd * n * I)
    end

    v0 = zeros(T, n_fourier * N)
    v0[n_max*N+1] = weight

    Pl = ilu(M, τ = 0.01)
    prob = LinearProblem(M, v0)
    ρtot = solve(prob, KrylovJL_GMRES(), Pl = Pl, abstol = tol, reltol = tol).u

    offset1 = n_max * N
    offset2 = (n_max + 1) * N
    ρ0 = reshape(ρtot[offset1+1:offset2], Ns, Ns)
    ρ0_tr = tr(ρ0)
    ρ0 = ρ0 / ρ0_tr
    ρ0 = QuantumObject((ρ0 + ρ0') / 2, dims = L_0.dims)
    ρtot = ρtot / ρ0_tr

    ρ_list = [ρ0]
    for i in 0:n_max-1
        ρi_m = reshape(ρtot[offset1-(i+1)*N+1:offset1-i*N], Ns, Ns)
        ρi_m = QuantumObject(ρi_m, dims = L_0.dims)
        push!(ρ_list, ρi_m)
    end

    return ρ_list
end

function _steadystate_floquet(
    L_0::QuantumObject{<:AbstractArray,SuperOperatorQuantumObject},
    L_p::QuantumObject{<:AbstractArray,SuperOperatorQuantumObject},
    L_m::QuantumObject{<:AbstractArray,SuperOperatorQuantumObject},
    ωd::Number,
    solver::SSFloquetEffectiveLiouvillian;
    n_max::Integer = 1,
    tol::R = 1e-8,
    kwargs...,
) where {R<:Real}
    ((L_0.dims == L_p.dims) && (L_0.dims == L_m.dims)) ||
        throw(ErrorException("The operators are not of the same Hilbert dimension."))

    L_eff = liouvillian_floquet(L_0, L_p, L_m, ωd; n_max = n_max, tol = tol)

    return steadystate(L_eff; solver = solver.steadystate_solver, kwargs...)
end
