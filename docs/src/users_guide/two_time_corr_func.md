# [Two-time Correlation Functions](@id doc:Two-time-Correlation-Functions)

## Introduction

With the `QuantumToolbox.jl` time-evolution function [`mesolve`](@ref), a state vector ([`Ket`](@ref)) or density matrix ([`Operator`](@ref)) can be evolved from an initial state at ``t_0`` to an arbitrary time ``t``, namely

```math
\hat{\rho}(t) = \mathcal{G}(t, t_0)\{\hat{\rho}(t_0)\},
```
where ``\mathcal{G}(t, t_0)\{\cdot\}`` is the propagator defined by the equation of motion. The resulting density matrix can then be used to evaluate the expectation values of arbitrary combinations of same-time operators.

To calculate two-time correlation functions on the form ``\left\langle \hat{A}(t+\tau) \hat{B}(t) \right\rangle``, we can use the quantum regression theorem [see, e.g., [Gardiner-Zoller2004](@citet)] to write

```math
\left\langle \hat{A}(t+\tau) \hat{B}(t) \right\rangle = \textrm{Tr} \left[\hat{A} \mathcal{G}(t+\tau, t)\{\hat{B}\hat{\rho}(t)\} \right] = \textrm{Tr} \left[\hat{A} \mathcal{G}(t+\tau, t)\{\hat{B} \mathcal{G}(t, 0)\{\hat{\rho}(0)\}\} \right],
```

We therefore first calculate ``\hat{\rho}(t) = \mathcal{G}(t, 0)\{\hat{\rho}(0)\}`` using [`mesolve`](@ref) with ``\hat{\rho}(0)`` as initial state, and then again use [`mesolve`](@ref) to calculate ``\mathcal{G}(t+\tau, t)\{\hat{B}\hat{\rho}(t)\}`` using ``\hat{B}\hat{\rho}(t)`` as initial state.

Note that if the initial state is the steady state, then ``\hat{\rho}(t) = \mathcal{G}(t, 0)\{\hat{\rho}_{\textrm{ss}}\} = \hat{\rho}_{\textrm{ss}}`` and

```math
\left\langle \hat{A}(t+\tau) \hat{B}(t) \right\rangle = \textrm{Tr} \left[\hat{A} \mathcal{G}(t+\tau, t)\{\hat{B}\hat{\rho}_{\textrm{ss}}\} \right] = \textrm{Tr} \left[\hat{A} \mathcal{G}(\tau, 0)\{\hat{B} \hat{\rho}_{\textrm{ss}}\} \right] = \left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle,
```
which is independent of ``t``, so that we only have one time coordinate ``\tau``.

`QuantumToolbox.jl` provides a family of functions that assists in the process of calculating two-time correlation functions. The available functions and their usage is shown in the table below.

| **Function call** | **Correlation function** |
|:------------------|:-------------------------|
| [`correlation_2op_2t`](@ref) | ``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle`` or ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle`` |
| [`correlation_2op_1t`](@ref) | ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` or ``\left\langle \hat{A}(0) \hat{B}(\tau) \right\rangle`` |
| [`correlation_3op_1t`](@ref) | ``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` |
| [`correlation_3op_2t`](@ref) | ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle`` |

The most common used case is to calculate the two time correlation function ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle``, which can be done by [`correlation_2op_1t`](@ref).

```@setup correlation_and_spectrum
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## Steadystate correlation function

The following code demonstrates how to calculate the ``\langle \hat{x}(t) \hat{x}(0)\rangle`` correlation for a leaky cavity with three different relaxation rates ``\gamma``.

```@example correlation_and_spectrum
tlist = LinRange(0, 10, 200)
a = destroy(10)
x = a' + a
H = a' * a

# if the initial state is specified as `nothing`, the steady state will be calculated and used as the initial state.
corr1 = correlation_2op_1t(H, nothing, tlist, [sqrt(0.5) * a], x, x)
corr2 = correlation_2op_1t(H, nothing, tlist, [sqrt(1.0) * a], x, x)
corr3 = correlation_2op_1t(H, nothing, tlist, [sqrt(2.0) * a], x, x)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = L"Time $t$", ylabel = L"\langle \hat{x}(t) \hat{x}(0) \rangle")
lines!(ax, tlist, real(corr1), label = L"\gamma = 0.5", linestyle = :solid)
lines!(ax, tlist, real(corr2), label = L"\gamma = 1.0", linestyle = :dash)
lines!(ax, tlist, real(corr3), label = L"\gamma = 2.0", linestyle = :dashdot)

axislegend(ax, position = :rt)

fig
```

## Emission spectrum

Given a correlation function ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle``, we can define the corresponding power spectrum as

```math
S(\omega) = \int_{-\infty}^\infty \left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle e^{-i \omega \tau} d \tau
```

In `QuantumToolbox.jl`, we can calculate ``S(\omega)`` using either [`spectrum`](@ref), which provides several solvers to perform the Fourier transform semi-analytically, or we can use the function [`spectrum_correlation_fft`](@ref) to numerically calculate the fast Fourier transform (FFT) of a given correlation data.

The following example demonstrates how these methods can be used to obtain the emission (``\hat{A} = \hat{a}^\dagger`` and ``\hat{B} = \hat{a}``) power spectrum.

```@example correlation_and_spectrum
N = 4             # number of cavity fock states
ωc = 1.0 * 2 * π  # cavity frequency
ωa = 1.0 * 2 * π  # atom frequency
g  = 0.1 * 2 * π  # coupling strength
κ  = 0.75         # cavity dissipation rate
γ  = 0.25         # atom dissipation rate

# Jaynes-Cummings Hamiltonian
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
H = ωc * a' * a + ωa * sm' * sm + g * (a' * sm + a * sm')

# collapse operators
n_th = 0.25
c_ops = [
    sqrt(κ * (1 + n_th)) * a,
    sqrt(κ *      n_th)  * a',
    sqrt(γ)              * sm,
];

# calculate the correlation function using mesolve, and then FFT to obtain the spectrum. 
# Here we need to make sure to evaluate the correlation function for a sufficient long time and 
# sufficiently high sampling rate so that FFT captures all the features in the resulting spectrum.
tlist = LinRange(0, 100, 5000)
corr = correlation_2op_1t(H, nothing, tlist, c_ops, a', a; progress_bar = Val(false))
ωlist1, spec1 = spectrum_correlation_fft(tlist, corr)

# calculate the power spectrum using spectrum
# using Exponential Series (default) method
ωlist2 = LinRange(0.25, 1.75, 200) * 2 * π
spec2 = spectrum(H, ωlist2, c_ops, a', a; solver = ExponentialSeries())

# calculate the power spectrum using spectrum
# using Pseudo-Inverse method
spec3 = spectrum(H, ωlist2, c_ops, a', a; solver = PseudoInverse())

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], title = "Vacuum Rabi splitting", xlabel = "Frequency", ylabel = "Emission power spectrum")
lines!(ax, ωlist1 / (2 * π), spec1, label = "mesolve + FFT", linestyle = :solid)
lines!(ax, ωlist2 / (2 * π), spec2, label = "Exponential Series", linestyle = :dash)
lines!(ax, ωlist2 / (2 * π), spec3, label = "Pseudo-Inverse", linestyle = :dashdot)

xlims!(ax, ωlist2[1] / (2 * π), ωlist2[end] / (2 * π))
axislegend(ax, position = :rt)

fig
```

## Non-steadystate correlation function

More generally, we can also calculate correlation functions of the kind ``\left\langle \hat{A}(t_1 + t_2) \hat{B}(t_1) \right\rangle``, i.e., the correlation function of a system that is not in its steady state. In `QuantumToolbox.jl`, we can evaluate such correlation functions using the function [`correlation_2op_2t`](@ref). The default behavior of this function is to return a matrix with the correlations as a function of the two time coordinates (``t_1`` and ``t_2``).

```@example correlation_and_spectrum
t1_list = LinRange(0, 10.0, 200)
t2_list = LinRange(0, 10.0, 200)

N = 10
a = destroy(N)
x = a' + a
H = a' * a

c_ops = [sqrt(0.25) * a]

α = 2.5
ρ0 = coherent_dm(N, α)

corr = correlation_2op_2t(H, ρ0, t1_list, t2_list, c_ops, x, x; progress_bar = Val(false))

# plot by CairoMakie.jl
fig = Figure(size = (500, 400))

ax = Axis(fig[1, 1], title = L"\langle \hat{x}(t_1 + t_2) \hat{x}(t_1) \rangle", xlabel = L"Time $t_1$", ylabel = L"Time $t_2$")

heatmap!(ax, t1_list, t2_list, real(corr))

fig
```

### Example: first-order optical coherence function

This example demonstrates how to calculate a correlation function on the form ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` for a non-steady initial state. Consider an oscillator that is interacting with a thermal environment. If the oscillator initially is in a coherent state, it will gradually decay to a thermal (incoherent) state. The amount of coherence can be quantified using the first-order optical coherence function

```math
g^{(1)}(\tau) = \frac{\left\langle \hat{a}^\dagger(\tau) \hat{a}(0) \right\rangle}{\sqrt{\left\langle \hat{a}^\dagger(\tau) \hat{a}(\tau) \right\rangle \left\langle \hat{a}^\dagger(0) \hat{a}(0)\right\rangle}}.
```
For a coherent state ``\vert g^{(1)}(\tau) \vert = 1``, and for a completely incoherent (thermal) state ``g^{(1)}(\tau) = 0``. The following code calculates and plots ``g^{(1)}(\tau)`` as a function of ``\tau``:

```@example correlation_and_spectrum
τlist = LinRange(0, 10, 200)

# Hamiltonian
N = 15
a = destroy(N)
H = 2 * π * a' * a

# collapse operator
G1 = 0.75
n_th = 2.00  # bath temperature in terms of excitation number
c_ops = [
    sqrt(G1 * (1 + n_th)) * a,
    sqrt(G1 *      n_th)  * a'
]

# start with a coherent state of α = 2.0
ρ0 = coherent_dm(N, 2.0)

# first calculate the occupation number as a function of time
n = mesolve(H, ρ0, τlist, c_ops, e_ops = [a' * a], progress_bar = Val(false)).expect[1,:]
n0 = n[1] # occupation number at τ = 0

# calculate the correlation function G1 and normalize with n to obtain g1
g1 = correlation_2op_1t(H, ρ0, τlist, c_ops, a', a, progress_bar = Val(false))
g1 = g1 ./ sqrt.(n .* n0)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], title = "Decay of a coherent state to an incoherent (thermal) state", xlabel = L"Time $\tau$")
lines!(ax, τlist, real(g1), label = L"g^{(1)}(\tau)", linestyle = :solid)
lines!(ax, τlist, real(n), label = L"n(\tau)", linestyle = :dash)

axislegend(ax, position = :rt)

fig
```

### Example: second-order optical coherence function

The second-order optical coherence function, with time-delay ``\tau``, is defined as

```math
g^{(2)}(\tau) = \frac{\left\langle \hat{a}^\dagger(0) \hat{a}^\dagger(\tau) \hat{a}(\tau) \hat{a}(0) \right\rangle}{\left\langle \hat{a}^\dagger(0) \hat{a}(0) \right\rangle^2}.
```

For a coherent state ``g^{(2)}(\tau) = 1``, for a thermal state ``g^{(2)}(\tau = 0) = 2`` and it decreases as a function of time (bunched photons, they tend to appear together), and for a Fock state with ``n``-photons ``g^{(2)}(\tau = 0) = n(n-1)/n^2 < 1`` and it increases with time (anti-bunched photons, more likely to arrive separated in time).

To calculate this type of correlation function with `QuantumToolbox.jl`, we can use [`correlation_3op_1t`](@ref), which computes a correlation function on the form ``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` (three operators and one delay-time vector). We first have to combine the central two operators into one single one as they are evaluated at the same time, e.g. here we do ``\hat{B}(\tau) = \hat{a}^\dagger(\tau) \hat{a}(\tau) = (\hat{a}^\dagger\hat{a})(\tau)``.

The following code calculates and plots ``g^{(2)}(\tau)`` as a function of ``\tau`` for a coherent, thermal and Fock state:

```@example correlation_and_spectrum
τlist = LinRange(0, 25, 200)

# Hamiltonian
N = 25
a = destroy(N)
H = 2 * π * a' * a

κ = 0.25
n_th = 2.0  # bath temperature in terms of excitation number
c_ops = [
    sqrt(κ * (1 + n_th)) * a,
    sqrt(κ *      n_th)  * a'
]

cases = [
    Dict("state" => coherent_dm(N, sqrt(2)), "label" => "coherent state", "lstyle" => :solid),
    Dict("state" => thermal_dm(N, 2), "label" => "thermal state", "lstyle" => :dash),
    Dict("state" => fock_dm(N, 2), "label" => "Fock state", "lstyle" => :dashdot),
]

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = L"Time $\tau$", ylabel = L"g^{(2)}(\tau)")

for case in cases
    ρ0 = case["state"]

    # calculate the occupation number at τ = 0
    n0 = expect(a' * a, ρ0)

    # calculate the correlation function g2
    g2 = correlation_3op_1t(H, ρ0, τlist, c_ops, a', a' * a, a, progress_bar = Val(false))
    g2 = g2 ./ n0^2

    lines!(ax, τlist, real(g2), label = case["label"], linestyle = case["lstyle"])
end

axislegend(ax, position = :rt)

fig
```
