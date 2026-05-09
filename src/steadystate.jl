export steadystate, steadystate_fourier
export SteadyStateSolver,
    SteadyStateDirectSolver,
    SteadyStateEigenSolver,
    SteadyStateLinearSolver,
    SteadyStateODESolver,
    SSFloquetEffectiveLiouvillian

abstract type SteadyStateSolver end

@doc raw"""
    SteadyStateDirectSolver(
        factorization = lu
    )

A solver which solves [`steadystate`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the standard method given in `LinearAlgebra`.

# Arguments
- `factorization::Function=lu`: The factorization method to use for solving the linear system. The default is `lu`, which performs LU factorization.
"""
Base.@kwdef struct SteadyStateDirectSolver{FT <: Function} <: SteadyStateSolver
    factorization::FT = lu
end

@doc raw"""
    SteadyStateEigenSolver()

A solver which solves [`steadystate`](@ref) by finding the zero (or lowest) eigenvalue of Liouvillian [`SuperOperator`](@ref) using [`eigsolve`](@ref).
"""
struct SteadyStateEigenSolver <: SteadyStateSolver end

@doc raw"""
    SteadyStateLinearSolver(
        alg = KrylovJL_GMRES(; precs = (A, p) -> A isa SparseMatrixCSC ? (ilu(A, τ = 0.01), I) : (I, I)),
        ρ0 = nothing
    )

A solver which solves [`steadystate`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the `alg`orithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/).

# Arguments
- `alg::SciMLLinearSolveAlgorithm=KrylovJL_GMRES()`: algorithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/)
- `ρ0::Union{Nothing, QuantumObject}=nothing`: The initial guess of the `steadystate` solution. If not specified, the initial guess will be handled by the solver.

# Note
Refer to [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/) for more details about the available algorithms. For example, the preconditioners can be defined directly in the solver like: `SteadyStateLinearSolver(alg = KrylovJL_GMRES(; precs = (A, p) -> (I, Diagonal(A))))`.
"""
Base.@kwdef struct SteadyStateLinearSolver{
        MT <: Union{SciMLLinearSolveAlgorithm, Nothing},
        ST <: Union{Nothing, QuantumObject},
    } <: SteadyStateSolver
    alg::MT = KrylovJL_GMRES(; precs = (A, p) -> A isa SparseMatrixCSC ? (ilu(A, τ = 0.01), I) : (I, I))
    ρ0::ST = nothing
end

@doc raw"""
    SteadyStateODESolver(
        alg = DP5(),
        ρ0 = nothing,
        tmax = Inf,
        terminate_reltol = 1e-4,
        terminate_abstol = 1e-6
    )

An ordinary differential equation (ODE) solver for solving [`steadystate`](@ref). It solves the stationary state based on [`mesolve`](@ref) with a termination condition.

The termination condition of the stationary state ``|\rho\rangle\!\rangle`` is that either the following condition is `true`:

```math
\lVert\frac{\partial |\hat{\rho}\rangle\!\rangle}{\partial t}\rVert \leq \textrm{reltol} \times\lVert\frac{\partial |\hat{\rho}\rangle\!\rangle}{\partial t}+|\hat{\rho}\rangle\!\rangle\rVert
```

or

```math
\lVert\frac{\partial |\hat{\rho}\rangle\!\rangle}{\partial t}\rVert \leq \textrm{abstol}
```

# Arguments
- `alg::AbstractODEAlgorithm=DP5()`: The algorithm to solve the ODE.
- `ρ0::Union{Nothing,QuantumObject}=nothing`: The initial state of the system. If not specified, a random density matrix state will be generated.
- `tmax::Real`: The final time step for the steady state problem. Default to `Inf`.
- `terminate_reltol`: The relative tolerance for stationary state terminate condition. Default to `1e-4`.
- `terminate_abstol`: The absolute tolerance for stationary state terminate condition. Default to `1e-6`.
- `return_details::Union{Val, Bool}`: Whether to return the details of the ODE solution. If `Val(true)`, the [`steadystate`](#ref) function will return a 2-element tuple: the steady state and an extra `NamedTuple` containing the ODE solution details. If `Val(false)`, only the steady state will be returned. Default to `Val(false)`.

!!! warning "Tolerances for terminate condition"
    The terminate condition tolerances `terminate_reltol` and `terminate_abstol` should be larger than `reltol` and `abstol` of [`mesolve`](@ref), respectively.

For more details about the solving `alg`orithms, please refer to [`OrdinaryDiffEq.jl`](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
"""
Base.@kwdef struct SteadyStateODESolver{
        MT <: AbstractODEAlgorithm,
        ST <: Union{Nothing, QuantumObject},
        TT <: Real,
        RT <: Real,
        AT <: Real,
        DT <: Union{Val, Bool},
    } <: SteadyStateSolver
    alg::MT = DP5()
    ρ0::ST = nothing
    tmax::TT = Inf
    terminate_reltol::RT = 1.0e-4
    terminate_abstol::AT = 1.0e-6
    return_details::DT = Val(false)
end

@doc raw"""
    SSFloquetEffectiveLiouvillian(steadystate_solver = SteadyStateDirectSolver())

A solver which solves [`steadystate_fourier`](@ref) by first extracting an effective time-independent Liouvillian and then using the `steadystate_solver` to extract the steadystate..

# Parameters
- `steadystate_solver::SteadyStateSolver=SteadyStateDirectSolver()`: The solver to use for the effective Liouvillian.

!!! note
    This solver is only available for [`steadystate_fourier`](@ref).
"""
Base.@kwdef struct SSFloquetEffectiveLiouvillian{SSST <: SteadyStateSolver} <: SteadyStateSolver
    steadystate_solver::SSST = SteadyStateDirectSolver()
end

@doc raw"""
    steadystate(
        H::AbstractQuantumObject{OpType},
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        solver::SteadyStateSolver = SteadyStateDirectSolver(),
        kwargs...,
    )

Solve the stationary state based on different solvers.

# Parameters
- `H`: The Hamiltonian or the Liouvillian of the system.
- `c_ops`: The list of the collapse operators.
- `solver`: see documentation [Solving for Steady-State Solutions](@ref doc:Solving-for-Steady-State-Solutions) for different solvers.
- `kwargs`: The keyword arguments for the solver.
"""
function steadystate(
        H::AbstractQuantumObject{OpType},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        solver::SteadyStateSolver = SteadyStateDirectSolver(),
        kwargs...,
    ) where {OpType <: Union{Operator, SuperOperator}}
    solver isa SSFloquetEffectiveLiouvillian && throw(
        ArgumentError(
            "The solver `SSFloquetEffectiveLiouvillian` is only available for the `steadystate_fourier` function.",
        ),
    )

    !isendomorphic(H.dimensions) && _non_endomorphic_dims_error("Hamiltonian or Liouvillian for steadystate", H.dimensions)

    L = liouvillian(H, c_ops)

    return _steadystate(L, solver; kwargs...)
end

# this function is used for SteadyStateDirectSolver and SteadyStateLinearSolver to handle the trace preserving condition
function _handle_steady_state_trace_preserving_condition(A::AbstractMatrix, N_op::Int, N_super::Int)
    # [ solving linear problem: Ax = b ]
    # handle A (Liouvillian) to include the trace preserving condition
    weight = norm(A, 1) / length(A)
    idx_range = collect(1:N_op)
    rows = _dense_similar(A, N_op)
    cols = _dense_similar(A, N_op)
    vals = _dense_similar(A, N_op)
    fill!(rows, 1)
    copyto!(cols, N_op .* (idx_range .- 1) .+ idx_range)
    fill!(vals, weight)
    Tn = _sparse_similar(A, rows, cols, vals, N_super, N_super)
    A_new = A + Tn

    # handle b to include the trace preserving condition
    b = _dense_similar(A_new, N_super)
    fill!(b, 0)
    allowed_setindex!(b, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    return A_new, b
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateLinearSolver; kwargs...)
    (haskey(kwargs, :Pl) || haskey(kwargs, :Pr)) && throw(ArgumentError("The use of preconditioners (Pl or Pr) must be defined in the solver."))

    # handle dimensions
    state_dimensions = L.dimensions.to.op_dims
    N_op = get_size(state_dimensions)[1]
    N_super = N_op^2

    # handle initial guess of u (steady state) and return dimensions
    # u0 can be useful for parameter sweeps when the steady state changes smoothly with the parameters
    u0 = if isnothing(solver.ρ0)
        nothing
    else
        _, u0_data, _, _ = _handle_init_state_and_sol_type_dims(L, solver.ρ0)
        u0_data
    end

    # linear problem: Au = b
    A, b = _handle_steady_state_trace_preserving_condition(L.data, N_op, N_super)

    prob = LinearProblem{true}(A, b, u0 = u0)
    ρss_vec = solve(prob, solver.alg; kwargs...).u

    ρss = reshape(ρss_vec, N_op, N_op)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator(), state_dimensions)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateEigenSolver; kwargs...)
    state_dimensions = L.dimensions.to.op_dims
    N = get_size(state_dimensions)[1]

    kwargs = merge((sigma = 1.0e-8, eigvals = 1), (; kwargs...))

    ρss_vec = eigsolve(L; kwargs...).vectors[:, 1]
    ρss = reshape(ρss_vec, N, N)
    ρss /= tr(ρss)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator(), state_dimensions)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateDirectSolver)
    # handle dimensions
    state_dimensions = L.dimensions.to.op_dims
    N_op = get_size(state_dimensions)[1]
    N_super = N_op^2

    # linear problem: Au = b
    A, b = _handle_steady_state_trace_preserving_condition(L.data, N_op, N_super)

    F = solver.factorization(A)
    ρss_vec = F \ b # This is still not supported on GPU, yet
    ρss = reshape(ρss_vec, N_op, N_op)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator(), state_dimensions)
end

function _steadystate(L::AbstractQuantumObject{SuperOperator}, solver::SteadyStateODESolver; kwargs...)
    ρ0 = isnothing(solver.ρ0) ? rand_dm(eltype(L), L.dimensions.to.op_dims) : solver.ρ0 # the validity of L and ρ0 will be checked in mesolve
    ftype = _float_type(ρ0)
    tlist = [ftype(0), ftype(solver.tmax)]

    # overwrite some kwargs and throw warning message to tell the users that we are ignoring these settings
    haskey(kwargs, :progress_bar) && @warn "Ignore keyword argument 'progress_bar' for SteadyStateODESolver"
    haskey(kwargs, :save_everystep) && @warn "Ignore keyword argument 'save_everystep' for SteadyStateODESolver"
    haskey(kwargs, :saveat) && @warn "Ignore keyword argument 'saveat' for SteadyStateODESolver"
    kwargs2 = merge(
        NamedTuple(kwargs), # we convert to NamedTuple just in case if kwargs is empty
        (progress_bar = Val(false), save_everystep = false, saveat = ftype[]),
    )

    # add terminate condition (callback)
    cb = TerminateSteadyState(
        solver.terminate_abstol,
        solver.terminate_reltol,
        SteadyStateODECondition(similar(mat2vec(ket2dm(ρ0)).data)),
    )
    kwargs3 = _merge_kwargs_with_callback(kwargs2, cb)

    sol = mesolve(L, ρ0, tlist; kwargs3...)
    ρss = sol.states[end]

    if getVal(solver.return_details)
        details = (
            t_final = sol.times_states[end], # sol.times_states relates to sol.t in OrdinaryDiffEq.jl
        )
        return ρss, details
    else
        return ρss
    end
end

_steadystate(
    L::QuantumObjectEvolution{SuperOperator},
    solver::T;
    kwargs...,
) where {T <: Union{SteadyStateDirectSolver, SteadyStateEigenSolver, SteadyStateLinearSolver}} =
    throw(ArgumentError("$(get_typename_wrapper(solver)) does not support QobjEvo."))

struct SteadyStateODECondition{CT <: AbstractArray}
    cache::CT
end

function (f::SteadyStateODECondition)(integrator, abstol, reltol, min_t)
    # this condition is same as DiffEqBase.NormTerminationMode

    f.cache .= (integrator.u .- integrator.uprev) ./ integrator.dt
    norm_du_dt = norm(f.cache)
    f.cache .+= integrator.u
    if norm_du_dt <= reltol * norm(f.cache) || norm_du_dt <= abstol
        return true
    else
        return false
    end
end

@doc raw"""
    steadystate_fourier(
        H_0::QuantumObject,
        H_p::QuantumObject,
        H_m::QuantumObject,
        ωd::Number,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        n_max::Integer = 2,
        tol::R = 1e-8,
        solver::FSolver = SteadyStateLinearSolver(),
        kwargs...,
    )

Calculates the steady state of a periodically driven system.
Here `H_0` is the Hamiltonian or the Liouvillian of the undriven system.
Considering a monochromatic drive at frequency ``\omega_d``, we divide it into two parts,
`H_p` and `H_m`, where `H_p` oscillates
as ``e^{i \omega t}`` and `H_m` oscillates as ``e^{-i \omega t}``.
There are two solvers available for this function:
- `SteadyStateLinearSolver`: Solves the linear system of equations.
- `SSFloquetEffectiveLiouvillian`: Solves the effective Liouvillian.
For both cases, `n_max` is the number of Fourier components to consider, and `tol` is the tolerance for the solver.

In the case of `SteadyStateLinearSolver`, the full linear system is solved at once:

```math
( \mathcal{L}_0 - i n \omega_d ) \hat{\rho}_n + \mathcal{L}_1 \hat{\rho}_{n-1} + \mathcal{L}_{-1} \hat{\rho}_{n+1} = 0
```

This is a tridiagonal linear system of the form

```math
\mathbf{A} \cdot \mathbf{b} = 0
```

where

```math
\mathbf{A} = \begin{pmatrix}
\mathcal{L}_0 - i (-n_{\textrm{max}}) \omega_{\textrm{d}} & \mathcal{L}_{-1} & 0 & \cdots & 0 \\
\mathcal{L}_1 & \mathcal{L}_0 - i (-n_{\textrm{max}}+1) \omega_{\textrm{d}} & \mathcal{L}_{-1} & \cdots & 0 \\
0 & \mathcal{L}_1 & \mathcal{L}_0 - i (-n_{\textrm{max}}+2) \omega_{\textrm{d}} & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \mathcal{L}_0 - i n_{\textrm{max}} \omega_{\textrm{d}}
\end{pmatrix}
```

and

```math
\mathbf{b} = \begin{pmatrix}
\hat{\rho}_{-n_{\textrm{max}}} \\
\hat{\rho}_{-n_{\textrm{max}}+1} \\
\vdots \\
\hat{\rho}_{0} \\
\vdots \\
\hat{\rho}_{n_{\textrm{max}}-1} \\
\hat{\rho}_{n_{\textrm{max}}}
\end{pmatrix}
```

This will allow to simultaneously obtain all the ``\hat{\rho}_n``.

In the case of `SSFloquetEffectiveLiouvillian`, instead, the effective Liouvillian is calculated using the matrix continued fraction method.

!!! note "different return"
    The two solvers returns different objects. The `SteadyStateLinearSolver` returns a list of [`QuantumObject`](@ref), containing the density matrices for each Fourier component (``\hat{\rho}_{-n}``, with ``n`` from ``0`` to ``n_\textrm{max}``), while the `SSFloquetEffectiveLiouvillian` returns only ``\hat{\rho}_0``. 

## Arguments
- `H_0::QuantumObject`: The Hamiltonian or the Liouvillian of the undriven system.
- `H_p::QuantumObject`: The Hamiltonian or the Liouvillian of the part of the drive that oscillates as ``e^{i \omega t}``.
- `H_m::QuantumObject`: The Hamiltonian or the Liouvillian of the part of the drive that oscillates as ``e^{-i \omega t}``.
- `ωd::Number`: The frequency of the drive.
- `c_ops::Union{Nothing,AbstractVector} = nothing`: The optional collapse operators.
- `n_max::Integer = 2`: The number of Fourier components to consider.
- `tol::R = 1e-8`: The tolerance for the solver.
- `solver::FSolver = SteadyStateLinearSolver`: The solver to use.
- `kwargs...`: Additional keyword arguments to be passed to the solver.
"""
function steadystate_fourier(
        H_0::QuantumObject{OpType1},
        H_p::QuantumObject{OpType2},
        H_m::QuantumObject{OpType3},
        ωd::Number,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        n_max::Integer = 2,
        tol::R = 1.0e-8,
        solver::FSolver = SteadyStateLinearSolver(),
        kwargs...,
    ) where {
        OpType1 <: Union{Operator, SuperOperator},
        OpType2 <: Union{Operator, SuperOperator},
        OpType3 <: Union{Operator, SuperOperator},
        R <: Real,
        FSolver <: SteadyStateSolver,
    }
    !isendomorphic(H_0.dimensions) && _non_endomorphic_dims_error("Hamiltonian or Liouvillian for steadystate_fourier", H_0.dimensions)
    check_dimensions(H_0, H_p, H_m)

    L_0 = liouvillian(H_0, c_ops)
    (L_0 isa QuantumObject) || throw(ArgumentError("steadystate_fourier only supports (time-independent) QuantumObject in c_ops"))

    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)
    return _steadystate_fourier(L_0, L_p, L_m, ωd, solver; n_max = n_max, tol = tol, kwargs...)
end

function _steadystate_fourier(
        L_0::QuantumObject{SuperOperator},
        L_p::QuantumObject{SuperOperator},
        L_m::QuantumObject{SuperOperator},
        ωd::Number,
        solver::SteadyStateLinearSolver;
        n_max::Integer = 1,
        tol::R = 1.0e-8,
        kwargs...,
    ) where {R <: Real}
    T1 = eltype(L_0)
    T2 = eltype(L_p)
    T3 = eltype(L_m)
    T = promote_type(T1, T2, T3)

    L_0_mat = get_data(L_0)
    L_p_mat = get_data(L_p)
    L_m_mat = get_data(L_m)

    state_dimensions = L_0.dimensions.to.op_dims
    N = size(L_0_mat, 1)
    Ns = get_size(state_dimensions)[1]
    n_fourier = 2 * n_max + 1
    n_list = (-n_max):n_max

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
    v0[n_max * N + 1] = weight

    (haskey(kwargs, :Pl) || haskey(kwargs, :Pr)) && error("The use of preconditioners must be defined in the solver.")
    !haskey(kwargs, :abstol) && (kwargs = merge((; kwargs...), (abstol = tol,)))
    !haskey(kwargs, :reltol) && (kwargs = merge((; kwargs...), (reltol = tol,)))

    prob = LinearProblem{true}(M, v0)
    ρtot = solve(prob, solver.alg; kwargs...).u

    offset1 = n_max * N
    offset2 = (n_max + 1) * N
    ρ0 = reshape(ρtot[(offset1 + 1):offset2], Ns, Ns)
    ρ0_tr = tr(ρ0)
    ρ0 = ρ0 / ρ0_tr
    ρ0 = QuantumObject((ρ0 + ρ0') / 2, type = Operator(), dims = state_dimensions)
    ρtot = ρtot / ρ0_tr

    ρ_list = [ρ0]
    for i in 0:(n_max - 1)
        ρi_m = reshape(ρtot[(offset1 - (i + 1) * N + 1):(offset1 - i * N)], Ns, Ns)
        ρi_m = QuantumObject(ρi_m, type = Operator(), dims = state_dimensions)
        push!(ρ_list, ρi_m)
    end

    return ρ_list
end

function _steadystate_fourier(
        L_0::QuantumObject{SuperOperator},
        L_p::QuantumObject{SuperOperator},
        L_m::QuantumObject{SuperOperator},
        ωd::Number,
        solver::SSFloquetEffectiveLiouvillian;
        n_max::Integer = 1,
        tol::R = 1.0e-8,
        kwargs...,
    ) where {R <: Real}
    L_eff = liouvillian_floquet(L_0, L_p, L_m, ωd; n_max = n_max, tol = tol)

    return steadystate(L_eff; solver = solver.steadystate_solver, kwargs...)
end
