# [Stochastic Solver](@id doc-TE:Stochastic-Solver)

When a quantum system is subjected to continuous measurement, through homodyne detection for example, it is possible to simulate the conditional quantum state using stochastic Schrödinger and master equations. The solution of these stochastic equations are quantum trajectories, which represent the conditioned evolution of the system given a specific measurement record.

## [Stochastic Schrödinger equation](@id doc-TE:Stochastic-Schrödinger-equation)

The stochastic Schrödinger time evolution of a quantum system is defined by the following stochastic differential equation [Wiseman2009Quantum; section 4.4](@cite):

```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 

```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{S}_n - \frac{1}{2} \hat{S}_n^\dagger \hat{S}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{S}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \psi(t) | \hat{S}_n + \hat{S}_n^\dagger | \psi(t) \rangle.
```

Above, ``\hat{H}`` is the Hamiltonian, ``\hat{S}_n`` are the stochastic collapse operators, and  ``dW_n(t)`` is the real Wiener increment (associated to ``\hat{S}_n``) which has the expectation values of ``E[dW_n]=0`` and ``E[dW_n^2]=dt``.

The solver [`ssesolve`](@ref) will construct the operators ``\hat{K}`` and ``\hat{M}_n``. Once the user passes the Hamiltonian (``\hat{H}``) and the stochastic collapse operators list (`sc_ops`; ``\{\hat{S}_n\}_n``). As with the [`mcsolve`](@ref), the number of trajectories and the random number generator for the noise realization can be fixed using the arguments: `ntraj` and `rng`, respectively.

## [Stochastic master equation](@id doc-TE:Stochastic-master-equation)

When the initial state of the system is a density matrix ``\rho(0)``, or when additional loss channels are included, the stochastic master equation solver [`smesolve`](@ref) must be used. The stochastic master equation is given by [Wiseman2009Quantum; section 4.4](@cite):

```math
d \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_i \mathcal{D}[\hat{C}_i] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,
```

The above implementation takes into account 2 types of collapse operators. ``\hat{C}_i`` (`c_ops`) represent the collapse operators related to pure dissipation, while ``\hat{S}_n`` (`sc_ops`) are the stochastic collapse operators. The first three terms on the right-hand side of the above equation is the deterministic part of the evolution which takes into account all operators ``\hat{C}_i`` and ``\hat{S}_n``. The last term (``\mathcal{H}[\hat{S}_n] \rho(t)``) is the stochastic part given solely by the operators ``\hat{S}_n``.


## [Example: Homodyne detection](@id doc-TE:Example:Homodyne-detection)

First, we solve the dynamics for an optical cavity at absolute zero (``0K``) whose output is monitored using homodyne detection. The cavity decay rate is given by ``\kappa`` and the ``\Delta`` is the cavity detuning with respect to the driving field. The homodyne current ``J_x`` is calculated using

```math
J_x = \langle \hat{x} \rangle + dW/dt,
```

where ``\hat{x}`` is the operator build from the `sc_ops` as

```math
\hat{x} = \hat{S} + \hat{S}^\dagger
```

```@setup stochastic-solve
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

```@example stochastic-solve
# parameters
N = 20         # Fock space dimension
Δ = 5 * 2 * π  # cavity detuning
κ = 2          # cavity decay rate
α = 4          # intensity of initial state
ntraj = 500    # number of trajectories

tlist = 0:0.0025:1

# operators
a = destroy(N)
x = a + a'
H = Δ * a' * a

# initial state
ψ0 = coherent(N, √α)

# temperature with average of 0 excitations (absolute zero)
n_th = 0
# c_ops  = [√(κ * n_th) * a'] -> nothing
sc_ops = [√(κ * (n_th + 1)) * a]
```

In this case, there is no additional dissipation (`c_ops = nothing`), and thus, we can use the [`ssesolve`](@ref):

```@example stochastic-solve
sse_sol = ssesolve(
    H,
    ψ0,
    tlist,
    sc_ops,
    e_ops = [x],
    ntraj = ntraj,
)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "Time")
#lines!(ax, tlist, real(sse_sol.xxxxxx), label = L"J_x", color = :red, linestyle = :solid) TODO: add this in the future
lines!(ax, tlist, real(sse_sol.expect[1,:]),  label = L"\langle x \rangle", color = :black, linestyle = :solid)

axislegend(ax, position = :rt)

fig
```

Next, we consider the same model but at a finite temperature to demonstrate [`smesolve`](@ref):

```@example stochastic-solve
# temperature with average of 1 excitations
n_th = 1
c_ops  = [√(κ * n_th) * a']
sc_ops = [√(κ * (n_th + 1)) * a]

sme_sol = smesolve(
    H,
    ψ0,
    tlist,
    c_ops,
    sc_ops,
    e_ops = [x],
    ntraj = ntraj,
)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "Time")
#lines!(ax, tlist, real(sme_sol.xxxxxx), label = L"J_x", color = :red, linestyle = :solid) TODO: add this in the future
lines!(ax, tlist, real(sme_sol.expect[1,:]),  label = L"\langle x \rangle", color = :black, linestyle = :solid)

axislegend(ax, position = :rt)

fig
```
