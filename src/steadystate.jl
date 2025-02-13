export steadystate, steadystate_fourier, steadystate_floquet
export SteadyStateSolver,
    SteadyStateDirectSolver,
    SteadyStateEigenSolver,
    SteadyStateLinearSolver,
    SteadyStateODESolver,
    SSFloquetEffectiveLiouvillian

abstract type SteadyStateSolver end

@doc raw"""
    SteadyStateDirectSolver()

A solver which solves [`steadystate`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the standard method given in `LinearAlgebra`.
"""
struct SteadyStateDirectSolver <: SteadyStateSolver end

@doc raw"""
    SteadyStateEigenSolver()

A solver which solves [`steadystate`](@ref) by finding the zero (or lowest) eigenvalue of Liouvillian [`SuperOperator`](@ref) using [`eigsolve`](@ref).
"""
struct SteadyStateEigenSolver <: SteadyStateSolver end

@doc raw"""
    SteadyStateLinearSolver(alg = KrylovJL_GMRES(), Pl = nothing, Pr = nothing)

A solver which solves [`steadystate`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the `alg`orithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/).

# Parameters
- `alg::SciMLLinearSolveAlgorithm=KrylovJL_GMRES()`: algorithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/)
- `Pl::Union{Function,Nothing}=nothing`: left preconditioner, see documentation [Solving for Steady-State Solutions](@ref doc:Solving-for-Steady-State-Solutions) for more details.
- `Pr::Union{Function,Nothing}=nothing`: right preconditioner, see documentation [Solving for Steady-State Solutions](@ref doc:Solving-for-Steady-State-Solutions) for more details.
"""
Base.@kwdef struct SteadyStateLinearSolver{
    MT<:Union{SciMLLinearSolveAlgorithm,Nothing},
    PlT<:Union{Function,Nothing},
    PrT<:Union{Function,Nothing},
} <: SteadyStateSolver
    alg::MT = KrylovJL_GMRES()
    Pl::PlT = nothing
    Pr::PrT = nothing
end

@doc raw"""
    SteadyStateODESolver(alg = Tsit5())

An ordinary differential equation (ODE) solver for solving [`steadystate`](@ref).

It includes a field (attribute) `SteadyStateODESolver.alg` that specifies the solving algorithm. Default to `Tsit5()`.

For more details about the solvers, please refer to [`OrdinaryDiffEq.jl`](https://docs.sciml.ai/OrdinaryDiffEq/stable/)
"""
Base.@kwdef struct SteadyStateODESolver{MT<:OrdinaryDiffEqAlgorithm} <: SteadyStateSolver
    alg::MT = Tsit5()
end

@doc raw"""
    SSFloquetEffectiveLiouvillian(steadystate_solver = SteadyStateDirectSolver())

A solver which solves [`steadystate_fourier`](@ref) by first extracting an effective time-independent Liouvillian and then using the `steadystate_solver` to extract the steadystate..

# Parameters
- `steadystate_solver::SteadyStateSolver=SteadyStateDirectSolver()`: The solver to use for the effective Liouvillian.

!!! note
    This solver is only available for [`steadystate_fourier`](@ref).
"""
Base.@kwdef struct SSFloquetEffectiveLiouvillian{SSST<:SteadyStateSolver} <: SteadyStateSolver
    steadystate_solver::SSST = SteadyStateDirectSolver()
end

@doc raw"""
    steadystate(
        H::QuantumObject{OpType},
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
    H::QuantumObject{OpType},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    solver::SteadyStateSolver = SteadyStateDirectSolver(),
    kwargs...,
) where {OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    solver isa SSFloquetEffectiveLiouvillian && throw(
        ArgumentError(
            "The solver `SSFloquetEffectiveLiouvillian` is only available for the `steadystate_fourier` function.",
        ),
    )

    L = liouvillian(H, c_ops)

    return _steadystate(L, solver; kwargs...)
end

function _steadystate(L::QuantumObject{SuperOperatorQuantumObject}, solver::SteadyStateLinearSolver; kwargs...)
    L_tmp = L.data
    N = prod(L.dimensions)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _get_dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _get_dense_similar(L_tmp, N)
    cols = _get_dense_similar(L_tmp, N)
    vals = _get_dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(vals, weight)
    Tn = sparse(rows, cols, vals, N^2, N^2)
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
    return QuantumObject(ρss, Operator, L.dimensions)
end

function _steadystate(L::QuantumObject{SuperOperatorQuantumObject}, solver::SteadyStateEigenSolver; kwargs...)
    N = prod(L.dimensions)

    kwargs = merge((sigma = 1e-8, k = 1), (; kwargs...))

    ρss_vec = eigsolve(L; kwargs...).vectors[:, 1]
    ρss = reshape(ρss_vec, N, N)
    ρss /= tr(ρss)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator, L.dimensions)
end

function _steadystate(L::QuantumObject{SuperOperatorQuantumObject}, solver::SteadyStateDirectSolver)
    L_tmp = L.data
    N = prod(L.dimensions)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _get_dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _get_dense_similar(L_tmp, N)
    cols = _get_dense_similar(L_tmp, N)
    vals = _get_dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(vals, weight)
    Tn = sparse(rows, cols, vals, N^2, N^2)
    L_tmp = L_tmp + Tn

    ρss_vec = L_tmp \ v0 # This is still not supported on GPU, yet
    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator, L.dimensions)
end

_steadystate(L::QuantumObject{SuperOperatorQuantumObject}, solver::SteadyStateODESolver; kwargs...) = throw(
    ArgumentError(
        "The initial state ψ0 is required for SteadyStateODESolver, use the following call instead: `steadystate(H, ψ0, tmax, c_ops)`.",
    ),
)

@doc raw"""
    steadystate(
        H::QuantumObject{HOpType},
        ψ0::QuantumObject{StateOpType},
        tmax::Real = Inf,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        solver::SteadyStateODESolver = SteadyStateODESolver(),
        reltol::Real = 1.0e-8,
        abstol::Real = 1.0e-10,
        kwargs...,
    )

Solve the stationary state based on time evolution (ordinary differential equations; `OrdinaryDiffEq.jl`) with a given initial state.

The termination condition of the stationary state ``|\rho\rangle\rangle`` is that either the following condition is `true`:

```math
\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}\rVert \leq \textrm{reltol} \times\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}+|\hat{\rho}\rangle\rangle\rVert
```

or

```math
\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}\rVert \leq \textrm{abstol}
```

# Parameters
- `H`: The Hamiltonian or the Liouvillian of the system.
- `ψ0`: The initial state of the system.
- `tmax=Inf`: The final time step for the steady state problem.
- `c_ops=nothing`: The list of the collapse operators.
- `solver`: see [`SteadyStateODESolver`](@ref) for more details.
- `reltol=1.0e-8`: Relative tolerance in steady state terminate condition and solver adaptive timestepping.
- `abstol=1.0e-10`: Absolute tolerance in steady state terminate condition and solver adaptive timestepping.
- `kwargs`: The keyword arguments for the ODEProblem.
"""
function steadystate(
    H::QuantumObject{HOpType},
    ψ0::QuantumObject{StateOpType},
    tmax::Real = Inf,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    solver::SteadyStateODESolver = SteadyStateODESolver(),
    reltol::Real = 1.0e-8,
    abstol::Real = 1.0e-10,
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    ftype = _FType(ψ0)
    cb = TerminateSteadyState(abstol, reltol, _steadystate_ode_condition)
    sol = mesolve(
        H,
        ψ0,
        [ftype(0), ftype(tmax)],
        c_ops,
        progress_bar = Val(false),
        save_everystep = false,
        saveat = ftype[],
        callback = cb,
    )

    ρss = sol.states[end]
    return ρss
end

function _steadystate_ode_condition(integrator, abstol, reltol, min_t)
    # this condition is same as DiffEqBase.NormTerminationMode

    du_dt = (integrator.u - integrator.uprev) / integrator.dt
    norm_du_dt = norm(du_dt)
    if (norm_du_dt <= reltol * norm(du_dt + integrator.u)) || (norm_du_dt <= abstol)
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

!!! note
    `steadystate_floquet` is a synonym of `steadystate_fourier`.

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
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    n_max::Integer = 2,
    tol::R = 1e-8,
    solver::FSolver = SteadyStateLinearSolver(),
    kwargs...,
) where {
    OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    R<:Real,
    FSolver<:SteadyStateSolver,
}
    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)
    return _steadystate_fourier(L_0, L_p, L_m, ωd, solver; n_max = n_max, tol = tol, kwargs...)
end

function _steadystate_fourier(
    L_0::QuantumObject{SuperOperatorQuantumObject},
    L_p::QuantumObject{SuperOperatorQuantumObject},
    L_m::QuantumObject{SuperOperatorQuantumObject},
    ωd::Number,
    solver::SteadyStateLinearSolver;
    n_max::Integer = 1,
    tol::R = 1e-8,
    kwargs...,
) where {R<:Real}
    T1 = eltype(L_0)
    T2 = eltype(L_p)
    T3 = eltype(L_m)
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

    (haskey(kwargs, :Pl) || haskey(kwargs, :Pr)) && error("The use of preconditioners must be defined in the solver.")
    if !isnothing(solver.Pl)
        kwargs = merge((; kwargs...), (Pl = solver.Pl(M),))
    elseif isa(M, SparseMatrixCSC)
        kwargs = merge((; kwargs...), (Pl = ilu(M, τ = 0.01),))
    end
    !isnothing(solver.Pr) && (kwargs = merge((; kwargs...), (Pr = solver.Pr(M),)))
    !haskey(kwargs, :abstol) && (kwargs = merge((; kwargs...), (abstol = tol,)))
    !haskey(kwargs, :reltol) && (kwargs = merge((; kwargs...), (reltol = tol,)))

    prob = LinearProblem(M, v0)
    ρtot = solve(prob, solver.alg; kwargs...).u

    offset1 = n_max * N
    offset2 = (n_max + 1) * N
    ρ0 = reshape(ρtot[offset1+1:offset2], Ns, Ns)
    ρ0_tr = tr(ρ0)
    ρ0 = ρ0 / ρ0_tr
    ρ0 = QuantumObject((ρ0 + ρ0') / 2, type = Operator, dims = L_0.dimensions)
    ρtot = ρtot / ρ0_tr

    ρ_list = [ρ0]
    for i in 0:n_max-1
        ρi_m = reshape(ρtot[offset1-(i+1)*N+1:offset1-i*N], Ns, Ns)
        ρi_m = QuantumObject(ρi_m, type = Operator, dims = L_0.dimensions)
        push!(ρ_list, ρi_m)
    end

    return ρ_list
end

function _steadystate_fourier(
    L_0::QuantumObject{SuperOperatorQuantumObject},
    L_p::QuantumObject{SuperOperatorQuantumObject},
    L_m::QuantumObject{SuperOperatorQuantumObject},
    ωd::Number,
    solver::SSFloquetEffectiveLiouvillian;
    n_max::Integer = 1,
    tol::R = 1e-8,
    kwargs...,
) where {R<:Real}
    check_dimensions(L_0, L_p, L_m)

    L_eff = liouvillian_floquet(L_0, L_p, L_m, ωd; n_max = n_max, tol = tol)

    return steadystate(L_eff; solver = solver.steadystate_solver, kwargs...)
end

# TODO: Synonym to align with QuTiP. Track https://github.com/qutip/qutip/issues/2632 when this can be removed.
const steadystate_floquet = steadystate_fourier
