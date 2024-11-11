# [Schrödinger Equation Solver](@id doc-TE:Schrödinger-Equation-Solver)

## [Unitary evolution](@id doc-TE:Unitary-evolution)

The dynamics of a closed (pure) quantum system is governed by the Schrödinger equation

```math
i\hbar\frac{\partial}{\partial t}\Psi(\vec{x}, t) = \hat{H}\Psi(\vec{x}, t),
```

where ``\Psi(\vec{x}, t)`` is the wave function, ``\hat{H}`` is the Hamiltonian, and ``\hbar`` is reduced Planck constant. In general, the Schrödinger equation is a partial differential equation (PDE) where both 
``\Psi`` and ``\hat{H}`` are functions of space ``\vec{x}`` and time ``t``. For computational purposes it is useful to expand the PDE in a set of basis functions that span the Hilbert space of the Hamiltonian, and to write the equation in matrix and vector form, namely

```math
i\hbar\frac{d}{dt}|\psi(t)\rangle = \hat{H}|\psi(t)\rangle,
```

where ``|\psi(t)\rangle`` is the state vector, and the Hamiltonian ``\hat{H}`` is now under matrix representation. This matrix equation can, in principle, be solved by diagonalizing the Hamiltonian matrix ``\hat{H}``. In practice, however, it is difficult to perform this diagonalization unless the size of the Hilbert space (dimension of the matrix ``\hat{H}``) is small. Analytically, it is a formidable task to calculate the dynamics for systems with more than two states. If, in addition, we consider dissipation due to the inevitable interaction with a surrounding environment, the computational complexity grows even larger, and we have to resort to numerical calculations in all realistic situations. This illustrates the importance of numerical calculations in describing the dynamics of open quantum systems, and the need for efficient and accessible tools for this task.

The Schrödinger equation, which governs the time-evolution of closed quantum systems, is defined by its Hamiltonian and state vector. In the previous sections, [Manipulating States and Operators](@ref doc:Manipulating-States-and-Operators) and [Tensor Products and Partial Traces](@ref doc:Tensor-products-and-Partial-Traces), we showed how Hamiltonians and state vectors are constructed in `QuantumToolbox.jl`. Given a Hamiltonian, we can calculate the unitary (non-dissipative) time-evolution of an arbitrary initial state vector ``|\psi(0)\rangle`` using the `QuantumToolbox` time evolution problem [`sesolveProblem`](@ref) or directly call the function [`sesolve`](@ref). It evolves the state vector ``|\psi(t)\rangle`` and evaluates the expectation values for a set of operators `e_ops` at each given time points, using an ordinary differential equation solver provided by the powerful julia package [`DifferentialEquation.jl`](https://docs.sciml.ai/DiffEqDocs/stable/).

## [Example: Spin dynamics](@id doc-TE:Example:Spin-dynamics)

```@setup sesolve
using QuantumToolbox
```

For example, the time evolution of a quantum spin-``\frac{1}{2}`` system (initialized in spin-``\uparrow``) with tunneling rate ``0.1`` is calculated, and the expectation values of the Pauli-Z operator ``\hat{\sigma}_z`` is also evaluated, with the following code

```@example sesolve
H = 2 * π * 0.1 * sigmax()
ψ0 = basis(2, 0) # spin-up
tlist = LinRange(0.0, 10.0, 20)

prob = sesolveProblem(H, ψ0, tlist, e_ops = [sigmaz()])
sol = sesolve(prob)
```

!!! note "Note"
    Here, we generate the time evolution problem by [`sesolveProblem`](@ref) first and then put it into the function [`sesolve`](@ref). One can also directly call [`sesolve`](@ref), which we also demonstrates in below.

The function returns [`TimeEvolutionSol`](@ref), as described in the previous section [Time Evolution Solutions](@ref doc-TE:Time-Evolution-Solutions). The attribute `expect` in `sol`ution is a list of expectation values for the operator(s) that are passed to the `e_ops` keyword argument. 

```@example sesolve
sol.expect
```

Passing multiple operators to `e_ops` as a `Vector` results in the expectation values for each operators at each time points.

```@example sesolve
tlist = LinRange(0.0, 10.0, 100)
sol = sesolve(H, ψ0, tlist, e_ops = [sigmaz(), sigmay()])
```

!!! note "Note"
    Here, we call [`sesolve`](@ref) directly instead of pre-defining [`sesolveProblem`](@ref) first (as shown previously).

```@example sesolve
times = sol.times
print(size(times))
```

```@example sesolve
expt = sol.expect
print(size(expt))
```

We can therefore plot the expectation values:

```@example sesolve
using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())

expt_z = real(expt[1,:])
expt_y = real(expt[2,:])

fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Expectation values")
lines!(ax, times, expt_z, label = L"\langle\hat{\sigma}_z\rangle", linestyle = :solid)
lines!(ax, times, expt_y, label = L"\langle\hat{\sigma}_y\rangle", linestyle = :dash)

axislegend(ax, position = :rb)

fig
```

If the keyword argument `e_ops` is not specified (or given as an empty `Vector`), the [`sesolve`](@ref) and functions return a [`TimeEvolutionSol`](@ref) that contains a list of state vectors which corresponds to the time points specified in `tlist`:

```@example sesolve
tlist = [0, 10]
sol = sesolve(H, ψ0, tlist) # or specify: e_ops = []

sol.states
```

This is because the `states` will be saved depend on the keyword argument `saveat` in `kwargs`. If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). 

One can also specify `e_ops` and `saveat` separately:

```@example sesolve
tlist = [0, 5, 10]
sol = sesolve(H, ψ0, tlist, e_ops = [sigmay()], saveat = tlist)
```

```@example sesolve
sol.expect
```

```@example sesolve
sol.states
```
