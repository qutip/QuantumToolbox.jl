# [Time Evolution Solutions](@id doc-TE:Time-Evolution-Solutions)

```@setup TE-solution
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## [Solution](@id doc-TE:Solution)
`QuantumToolbox` utilizes the powerful [`DifferentialEquation.jl`](https://docs.sciml.ai/DiffEqDocs/stable/) to simulate different kinds of quantum system dynamics. Thus, we will first look at the data structure used for returning the solution (`sol`) from [`DifferentialEquation.jl`](https://docs.sciml.ai/DiffEqDocs/stable/). The solution stores all the crucial data needed for analyzing and plotting the results of a simulation. A generic structure [`TimeEvolutionSol`](@ref) contains the following properties for storing simulation data:

| **Fields (Attributes)** | **Description** |
|:------------------------|:----------------|
| `sol.times` | The list of time points at which the expectation values are calculated during the evolution. |
| `sol.times_states` | The list of time points at which the states are stored during the evolution. |
| `sol.states` | The list of result states corresponding to each time point in `sol.times_states`. |
| `sol.expect` | The expectation values corresponding to each time point in `sol.times`. |
| `sol.alg` | The algorithm which is used during the solving process. |
| `sol.abstol` | The absolute tolerance which is used during the solving process. |
| `sol.reltol` | The relative tolerance which is used during the solving process. |
| `sol.retcode` (or `sol.converged`) | The returned status from the solver. |

## [Accessing data in solutions](@id doc-TE:Accessing-data-in-solutions)

To understand how to access the data in solution, we will use an example as a guide, although we do not worry about the simulation details at this stage. The Schrödinger equation solver ([`sesolve`](@ref)) used in this example returns [`TimeEvolutionSol`](@ref):

```@example TE-solution
H = 0.5 * sigmay()
ψ0 = basis(2, 0)
e_ops = [
    proj(basis(2, 0)),
    proj(basis(2, 1)),
    basis(2, 0) * basis(2, 1)'
]
tlist = LinRange(0, 10, 100)
sol = sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
nothing # hide
```

To see what is contained inside the solution, we can use the `print` function:

```@example TE-solution
print(sol)
```

It tells us the number of expectation values are computed and the number of states are stored. Now we have all the information needed to analyze the simulation results. To access the data for the three expectation values, one can do:

```@example TE-solution
expt1 = real(sol.expect[1,:])
expt2 = real(sol.expect[2,:])
expt3 = real(sol.expect[3,:])
nothing # hide
```

Recall that `Julia` uses `Fortran`-style indexing that begins with one (i.e., `[1,:]` represents the 1-st observable, where `:` represents all values corresponding to `tlist`).

Together with the list of time points at which these expectation values are calculated:

```@example TE-solution
times = sol.times
nothing # hide
```

we can plot the resulting expectation values:

```@example TE-solution
# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = L"t")
lines!(ax, times, expt1, label = L"\langle 0 | \rho(t) | 0 \rangle")
lines!(ax, times, expt2, label = L"\langle 1 | \rho(t) | 1 \rangle")
lines!(ax, times, expt3, label = L"\langle 0 | \rho(t) | 1 \rangle")

ylims!(ax, (-0.5, 1.0))
axislegend(ax, position = :lb)

fig
```

State vectors, or density matrices, are accessed in a similar manner:

```@example TE-solution
sol.states
```

Together with the list of time points at which these states are stored:

```@example TE-solution
times = sol.times_states
nothing # hide
```

Here, the solution contains only one (final) state. Because the `states` will be saved depend on the keyword argument `saveat` in `kwargs`. If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). One can also specify `e_ops` and `saveat` separately.

Some other solvers can have other output.

## [Multiple trajectories solution](@id doc-TE:Multiple-trajectories-solution)

The solutions are different for solvers which compute multiple trajectories, such as the [`TimeEvolutionMCSol`](@ref) (Monte Carlo) or the [`TimeEvolutionStochasticSol`](@ref) (stochastic methods). The storage of expectation values and states depends on the keyword argument `keep_runs_results`, which determines whether the results of all trajectories are stored in the solution.

When the keyword argument `keep_runs_results` is passed as `Val(false)` to a multi-trajectory solver, the `states` and `expect` fields store only the average results (averaged over all trajectories). The results can be accessed by the following index-order:

```julia
sol.states[time_idx]
sol.expect[e_op,time_idx]
```

For example:

```@example TE-solution
tlist = LinRange(0, 1, 11)
c_ops = (destroy(2),)
e_ops = (num(2),)

sol_mc1 = mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, ntraj=25, keep_runs_results=Val(false), progress_bar=Val(false))

size(sol_mc1.expect)
```

If the keyword argument `keep_runs_results = Val(true)`, the results for each trajectory and each time are stored, and the index-order of the elements in fields `states` and `expect` are:


```julia
sol.states[trajectory,time_idx]
sol.expect[e_op,trajectory,time_idx]
```

For example:

```@example TE-solution
sol_mc2 = mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, ntraj=25, keep_runs_results=Val(true), progress_bar=Val(false))

size(sol_mc2.expect)
```

We also provide the following functions for statistical analysis of multi-trajectory `sol`utions:

| **Functions** | **Description** |
|:------------|:----------------|
| [`average_states(sol)`](@ref average_states) | Return the trajectory-averaged result states (as density [`Operator`](@ref)) at each time point. |
| [`average_expect(sol)`](@ref average_expect) | Return the trajectory-averaged expectation values at each time point. |
| [`std_expect(sol)`](@ref std_expect) | Return the trajectory-wise standard deviation of the expectation values at each time point. |
