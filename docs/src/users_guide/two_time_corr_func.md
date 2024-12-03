# [Two-time Correlation Functions](@id doc:Two-time-Correlation-Functions)

## Introduction

With the `QuantumToolbox.jl` time-evolution function [`mesolve`](@ref), a state vector ([`Ket`](@ref)) or density matrix ([`Operator`](@ref)) can be evolved from an initial state at ``t_0`` to an arbitrary time ``t``, namely

```math
\hat{\rho}(t) = \mathcal{G}(t, t_0)\{\hat{\rho}(t_0)\},
```
where ``\mathcal{G}(t, t_0)\{\cdot\}`` is the propagator defined by the equation of motion. The resulting density matrix can then be used to evaluate the expectation values of arbitrary combinations of same-time operators.

To calculate two-time correlation functions on the form ``\left\langle \hat{A}(t+\tau) \hat{B}(t) \right\rangle``, we can use the quantum regression theorem (see, e.g., [Gardiner-Zoller2004](@cite)) to write

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

Given a correlation function ``\langle \hat{A}(\tau) \hat{B}(0) \rangle``, we can define the corresponding power spectrum as

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

The following part of this page is still under construction, please visit [API](@ref doc-API) first.

### Example: first-order optical coherence function

### Example: second-order optical coherence function
