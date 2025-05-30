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

# Arguments
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
    SteadyStateODESolver(
        alg = Tsit5(),
        ψ0 = nothing,
        tmax = Inf,
        )

An ordinary differential equation (ODE) solver for solving [`steadystate`](@ref).

Solve the stationary state based on time evolution (ordinary differential equations; `OrdinaryDiffEq.jl`) with a given initial state.

The termination condition of the stationary state ``|\rho\rangle\rangle`` is that either the following condition is `true`:

```math
\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}\rVert \leq \textrm{reltol} \times\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}+|\hat{\rho}\rangle\rangle\rVert
```

or

```math
\lVert\frac{\partial |\hat{\rho}\rangle\rangle}{\partial t}\rVert \leq \textrm{abstol}
```

# Arguments
- `alg::OrdinaryDiffEqAlgorithm=Tsit5()`: The algorithm to solve the ODE.
- `ψ0::Union{Nothing,QuantumObject}=nothing`: The initial state of the system. If not specified, a random pure state will be generated.
- `tmax::Real=Inf`: The final time step for the steady state problem.

For more details about the solvers, please refer to [`OrdinaryDiffEq.jl`](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
"""
Base.@kwdef struct SteadyStateODESolver{MT<:OrdinaryDiffEqAlgorithm,ST<:Union{Nothing,QuantumObject},T<:Real} <:
                   SteadyStateSolver
    alg::MT = Tsit5()
    ψ0::ST = nothing
    tmax::T = Inf
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
) where {OpType<:Union{Operator,SuperOperator}}
    solver isa SSFloquetEffectiveLiouvillian && throw(
        ArgumentError(
            "The solver `SSFloquetEffectiveLiouvillian` is only available for the `steadystate_fourier` function.",
        ),
    )

    L = liouvillian(H, c_ops)

    return _steadystate(L, solver; kwargs...)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateLinearSolver; kwargs...)
    L_tmp = L.data
    N = prod(L.dimensions)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _dense_similar(L_tmp, N)
    cols = _dense_similar(L_tmp, N)
    vals = _dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(vals, weight)
    Tn = _sparse_similar(L_tmp, rows, cols, vals, N^2, N^2)
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
    return QuantumObject(ρss, Operator(), L.dimensions)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateEigenSolver; kwargs...)
    N = prod(L.dimensions)

    kwargs = merge((sigma = 1e-8, eigvals = 1), (; kwargs...))

    ρss_vec = eigsolve(L; kwargs...).vectors[:, 1]
    ρss = reshape(ρss_vec, N, N)
    ρss /= tr(ρss)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator(), L.dimensions)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateDirectSolver)
    L_tmp = L.data
    N = prod(L.dimensions)
    weight = norm(L_tmp, 1) / length(L_tmp)

    v0 = _dense_similar(L_tmp, N^2)
    fill!(v0, 0)
    allowed_setindex!(v0, weight, 1) # Because scalar indexing is not allowed on GPU arrays

    idx_range = collect(1:N)
    rows = _dense_similar(L_tmp, N)
    cols = _dense_similar(L_tmp, N)
    vals = _dense_similar(L_tmp, N)
    fill!(rows, 1)
    copyto!(cols, N .* (idx_range .- 1) .+ idx_range)
    fill!(vals, weight)
    Tn = sparse(rows, cols, vals, N^2, N^2)
    L_tmp = L_tmp + Tn

    ρss_vec = L_tmp \ v0 # This is still not supported on GPU, yet
    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    return QuantumObject(ρss, Operator(), L.dimensions)
end

function _steadystate(L::QuantumObject{SuperOperator}, solver::SteadyStateODESolver; kwargs...)
    tmax = solver.tmax

    ψ0 = isnothing(solver.ψ0) ? rand_ket(L.dimensions) : solver.ψ0
    abstol = haskey(kwargs, :abstol) ? kwargs[:abstol] : DEFAULT_ODE_SOLVER_OPTIONS.abstol
    reltol = haskey(kwargs, :reltol) ? kwargs[:reltol] : DEFAULT_ODE_SOLVER_OPTIONS.reltol

    ftype = _FType(ψ0)
    _terminate_func = SteadyStateODECondition(similar(mat2vec(ket2dm(ψ0)).data))
    cb = TerminateSteadyState(abstol, reltol, _terminate_func)
    sol = mesolve(
        L,
        ψ0,
        [ftype(0), ftype(tmax)],
        progress_bar = Val(false),
        save_everystep = false,
        saveat = ftype[],
        callback = cb,
    )

    ρss = sol.states[end]
    return ρss
end

struct SteadyStateODECondition{CT<:AbstractArray}
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
    OpType1<:Union{Operator,SuperOperator},
    OpType2<:Union{Operator,SuperOperator},
    OpType3<:Union{Operator,SuperOperator},
    R<:Real,
    FSolver<:SteadyStateSolver,
}
    L_0 = liouvillian(H_0, c_ops)
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
    tol::R = 1e-8,
    kwargs...,
) where {R<:Real}
    L0_mat = get_data(L_0)
    Lp_mat = get_data(L_p)
    Lm_mat = get_data(L_m)
    T = promote_type(eltype(L0_mat), eltype(Lp_mat), eltype(Lm_mat))
    N = size(L0_mat, 1)
    Ns = isqrt(N)
    n_fourier = 2 * n_max + 1
    n_list = (-n_max):n_max
    weight = one(T)
    diag_indices = Ns * (0:(Ns-1)) .+ (1:Ns)
    Mn = _sparse_similar(L0_mat, ones(Int, Ns), diag_indices, fill(weight, Ns), N, N)
    L = L0_mat + Mn
    Kp = _sparse_similar(L0_mat, spdiagm(-1 => ones(T, n_fourier - 1)))
    Km = _sparse_similar(L0_mat, spdiagm(1 => ones(T, n_fourier - 1)))
    M = kron(Kp, Lm_mat) + kron(Km, Lp_mat)
    n_vals = -1im * ωd * T.(n_list)
    I_N = sparse(I, N, N)
    I_F = spdiagm(0 => ones(T, n_fourier))
    D_F = spdiagm(0 => n_vals)
    block_diag = _spre(L, I_F) + _spost(D_F, I_N)
    M += block_diag
    v0 = zeros(T, n_fourier * N)
    allowed_setindex!(v0, weight, n_max * N + 1)
    (haskey(kwargs, :Pl) || haskey(kwargs, :Pr)) &&
        error("The use of preconditioners must be defined in the solver.")
    if !isnothing(solver.Pl)
        kwargs = (; kwargs..., Pl = solver.Pl(M))
    elseif isa(M, SparseMatrixCSC)
        kwargs = (; kwargs..., Pl = ilu(M, τ = 0.01))
    end
    if !isnothing(solver.Pr)
        kwargs = (; kwargs..., Pr = solver.Pr(M))
    end
    kwargs = merge((; abstol = tol, reltol = tol), kwargs)
    prob = LinearProblem(M, v0)
    ρtot = solve(prob, solver.alg; kwargs...).u
    offset1 = n_max * N
    offset2 = (n_max + 1) * N
    ρ0 = Matrix(reshape(view(ρtot, (offset1 + 1):offset2), Ns, Ns))
    ρ0 ./= tr(ρ0)
    ρ0 = QuantumObject((ρ0 + ρ0') / 2, type = Operator(), dims = L_0.dimensions)
    idx_ranges = [(offset1 - (i + 1) * N + 1):(offset1 - i * N) for i in 0:(n_max - 1)]
    ρ_components = map(idx_range ->
        QuantumObject(Matrix(reshape(view(ρtot, idx_range), Ns, Ns)),
                      type = Operator(),
                      dims = L_0.dimensions),
        idx_ranges)
    return vcat([ρ0], ρ_components)
end

function _steadystate_fourier(
    L_0::QuantumObject{SuperOperator},
    L_p::QuantumObject{SuperOperator},
    L_m::QuantumObject{SuperOperator},
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
