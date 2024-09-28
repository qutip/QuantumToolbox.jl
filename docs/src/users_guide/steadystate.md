# [Solving for Steady-State Solutions](@id doc:Solving-for-Steady-State-Solutions)

## Introduction

For time-independent open quantum systems with decay rates larger than the corresponding excitation rates, the system will tend toward a steady state as ``t\rightarrow\infty`` that satisfies the equation

```math
\frac{d\hat{\rho}_{\textrm{ss}}}{dt} = \mathcal{L}\hat{\rho}_{\textrm{ss}}=0.
```

Although the requirement for time-independence seems quite restrictive, one can often employ a transformation to the interaction picture that yields a time-independent Hamiltonian. For many these systems, solving for the asymptotic density matrix ``\hat{\rho}_{\textrm{ss}}`` can be achieved using direct or iterative solution methods faster than using master equation or Monte-Carlo simulations. Although the steady state equation has a simple mathematical form, the properties of the Liouvillian operator are such that the solutions to this equation are anything but straightforward to find.

## Steady State solvers in `QuantumToolbox.jl`
In `QuantumToolbox.jl`, the steady-state solution for a system Hamiltonian or Liouvillian is given by [`steadystate`](@ref). This function implements a number of different methods for finding the steady state, each with their own pros and cons, where the method used can be chosen using the `solver` keyword argument.

| **Solver** | **Description** |
|:-----------|:----------------|
| [`SteadyStateDirectSolver()`](@ref SteadyStateDirectSolver) | Directly solve ``Ax=b`` using the standard method given in `Julia` `LinearAlgebra` |
| [`SteadyStateLinearSolver()`](@ref SteadyStateLinearSolver) | Directly solve ``Ax=b`` using the algorithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/) |
| [`SteadyStateEigenSolver()`](@ref SteadyStateEigenSolver) | Find the zero (or lowest) eigenvalue of ``\mathcal{L}`` using [`eigsolve`](@ref) |
| [`SteadyStateODESolver()`](@ref SteadyStateODESolver) | Solving time evolution with algorithms given in [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) |

## Using Steady State solvers

The function [`steadystate`](@ref) can take either a Hamiltonian and a list of collapse operators `c_ops` as input, generating internally the corresponding Liouvillian ``\mathcal{L}`` in Lindblad form, or alternatively, a Liouvillian ``\mathcal{L}`` passed by the user.

```julia
ρ_ss = steadystate(H)        # Hamiltonian
ρ_ss = steadystate(H, c_ops) # Hamiltonian and collapse operators
ρ_ss = steadystate(L)        # Liouvillian
```
where `H` is a quantum object representing the system Hamiltonian ([`Operator`](@ref)) or Liouvillian ([`SuperOperator`](@ref)), and `c_ops` (defaults to `nothing`) can be a list of [`QuantumObject`](@ref) for the system collapse operators. The output, labelled as `ρ_ss`, is the steady-state solution for the systems. If no other keywords are passed to the function, the default solver [`SteadyStateDirectSolver()`](@ref SteadyStateDirectSolver) is used.

To change a solver, use the keyword argument `solver`, for example:

```julia
ρ_ss = steadystate(H, c_ops; solver = SteadyStateLinearSolver())
```

!!! note "Initial state for `SteadyStateODESolver()`"
    It is necessary to provide the initial state `ψ0` for ODE solving method, namely
    `steadystate(H, ψ0, tspan, c_ops, solver = SteadyStateODESolver())`, where `tspan::Real` represents the final time step, defaults to `Inf` (infinity).

Although it is not obvious, the [`SteadyStateDirectSolver()`](@ref SteadyStateDirectSolver) and [`SteadyStateEigenSolver()`](@ref SteadyStateEigenSolver) methods all use an LU decomposition internally and thus can have a large memory overhead. In contrast, for [`SteadyStateLinearSolver()`](@ref SteadyStateLinearSolver), iterative algorithms provided by [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/), such as `KrylovJL_GMRES()` and `KrylovJL_BICGSTAB()`, do not factor the matrix and thus take less memory than the LU methods and allow, in principle, for extremely large system sizes. The downside is that these methods can take much longer than the direct method as the condition number of the Liouvillian matrix is large, indicating that these iterative methods require a large number of iterations for convergence. 

To overcome this, one can provide preconditioner that solves for an approximate inverse for the (modified) Liouvillian, thus better conditioning the problem, leading to faster convergence. The left and right preconditioner can be specified as the keyword argument `Pl` and `Pr`, respectively:
```julia
solver = SteadyStateLinearSolver(alg = KrylovJL_GMRES(), Pl = Pl, Pr = Pr)
```
The use of a preconditioner can actually make these iterative methods faster than the other solution methods. The problem with precondioning is that it is only well defined for Hermitian matrices. Since the Liouvillian is non-Hermitian, the ability to find a good preconditioner is not guaranteed. Moreover, if a preconditioner is found, it is not guaranteed to have a good condition number. 

Furthermore, `QuantiumToolbox` can take advantage of the Intel (Pardiso) LU solver in the Intel Math Kernel library (Intel-MKL) by using [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/) and either [`Pardiso.jl`](https://github.com/JuliaSparse/Pardiso.jl) or [`MKL_jll.jl`](https://github.com/JuliaBinaryWrappers/MKL_jll.jl):

```julia
using QuantumToolbox
using LinearSolve # must be loaded

using Pardiso
solver = SteadyStateLinearSolver(alg = MKLPardisoFactorize())

using MKL_jll
solver = SteadyStateLinearSolver(alg = MKLLUFactorization())
```

See [`LinearSolve.jl` Solvers](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/) for more details.

## Example: Harmonic Oscillator in Thermal Bath

```@example steady_state_example
using QuantumToolbox
using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())


# Define parameters
N = 20  # number of basis states to consider
a = destroy(N)
H = a' * a
ψ0 = basis(N, 10)  # initial state
κ = 0.1  # coupling to oscillator
n_th = 2  # temperature with average of 2 excitations

# collapse operators 
# c_op_list = [        emission        ;     absorption      ]
c_op_list = [ sqrt(κ * (1 + n_th)) * a ; sqrt(κ * n_th) * a' ]

# find steady-state solution
ρ_ss = steadystate(H, c_op_list)

# find expectation value for particle number in steady state
e_ops = [a' * a]
exp_ss = real(expect(e_ops[1], ρ_ss))

tlist = LinRange(0, 50, 100)

# monte-carlo
sol_mc = mcsolve(H, ψ0, tlist, c_op_list, e_ops=e_ops, n_traj=100, progress_bar=false)
exp_mc = real(sol_mc.expect[1, :])

# master eq.
sol_me = mesolve(H, ψ0, tlist, c_op_list, e_ops=e_ops, progress_bar=false)
exp_me = real(sol_me.expect[1, :])

# plot the results
fig = Figure(size = (800, 400), fontsize = 15)
ax = Axis(fig[1, 1], 
    title = L"Decay of Fock state $|10\rangle$ in a thermal environment with $\langle n\rangle=2$",
    xlabel = "Time", 
    ylabel = "Number of excitations", 
    titlesize = 24,
    xlabelsize = 20, 
    ylabelsize = 20
)
lines!(ax, tlist, exp_mc, label = "Monte-Carlo", linewidth = 2, color = :blue)
lines!(ax, tlist, exp_me, label = "Master Equation", linewidth = 2, color = :orange, linestyle = :dash)
lines!(ax, tlist, exp_ss .* ones(length(tlist)), label = "Steady State", linewidth = 2, color = :red)
axislegend(ax, position = :rt)

fig
```

## Calculate steady state for periodically driven systems

See the docstring of [`steadystate_floquet`](@ref) for more details.
