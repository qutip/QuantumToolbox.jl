# QuantumToolbox.jl Logo Tutorial

## Introduction

In this tutorial, we will demonstrate how to create the logo for the **QuantumToolbox.jl** package. The logo represents the Wigner function of the triangular cat state, which is a linear superposition of three coherent states. The resulting Wigner function has a triangular shape that resembles the Julia logo. We will also define a custom colormap that varies based on the value of the Wigner function and the spatial coordinates, such that the three blobs corresponding to the coherent states have different colors (matching the colors of the Julia logo).

### Triangular Cat State

A cat state, often referred to as a Schrödinger cat state, is a quantum state that is a superposition of two coherent states with opposite phases:

```math
| \psi_{\text{cat}} \rangle = \frac{1}{\sqrt{2}} \left( | \alpha \rangle + | -\alpha \rangle \right)
```

where ``| \alpha \rangle`` is a coherent state with amplitude ``\alpha``.

The triangular cat state is a generalization of the standard cat state. It is a superposition of three coherent states with phases ``\theta_0, \theta_1, \theta_2``separated by ``120^\circ``(or ``2\pi/3``radians):

```math
| \psi_{\text{tri-cat}} \rangle = \frac{1}{\sqrt{3}} \left( | \alpha_0 \rangle + | \alpha_1 \rangle + | \alpha_2 \rangle \right)
```

where ``\alpha_j = \rho e^{i\theta_j}``with ``\theta_j = \frac{\pi}{2} + \frac{2\pi j}{3}``and ``j = 0, 1, 2``.

### Wigner Function

The Wigner function ``W(x, p)``is a quasi-probability distribution used in quantum mechanics to represent quantum states in phase space. It is defined as:

```math
W(x, p) = \frac{1}{\pi \hbar} \int_{-\infty}^{\infty} \psi^*(x + y) \psi(x - y) e^{2ipy / \hbar} \, dy
```

where ``\psi(x)``is the wave function of the quantum state, ``x``is the position, ``p``is the momentum, and ``\hbar``is the reduced Planck constant. Unlike classical probability distributions, the Wigner function can take negative values, which indicates non-classical behavior.

## Generating the Logo

First, let's load the required packages:

```@example logo
using QuantumToolbox
using CairoMakie
CairoMakie.activate!(type = "svg", pt_per_unit = 1)
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

### Parameters

Here we define the parameters for the triangular cat state:

```@example logo
N = 30  # Cutoff of the Hilbert space for the harmonic oscillator
ρ = 2.5  # Amplitude of the coherent state
θ1 = π / 2
θ2 = π / 2 + 2π / 3
θ3 = π / 2 + 4π / 3
α1 = ρ * exp(im * θ1)
α2 = ρ * exp(im * θ2)
α3 = ρ * exp(im * θ3)
```

### Constructing the State

Next, we construct the triangular cat state as a normalized superposition of three coherent states:

```@example logo
ψ = coherent(N, α1) + coherent(N, α2) + coherent(N, α3)
normalize!(ψ)
```

### Defining the Grid and calculating the Wigner function

We define the grid for the Wigner function and calculate it using the [`wigner`](@ref) function. We shift the grid in the imaginary direction to ensure that the Wigner function is centered around the origin of the figure. The [`wigner`](@ref) function also supports the `g` scaling factor, which we put here equal to ``2``.

```@example logo
xvec = range(-ρ, ρ, 500) .* 1.5
yvec = xvec .+ (abs(imag(α1)) - abs(imag(α2))) / 2

wig = wigner(ψ, xvec, yvec, g = 2)
```

### Plotting the Wigner function

Finally, we plot the Wigner function using the `heatmap` function from the `CairoMakie` package.

```@example logo
fig = Figure(size = (500, 500), figure_padding = 0)
ax = Axis(fig[1, 1])
heatmap!(ax, xvec, yvec, wig', colormap = :RdBu, interpolate = true, rasterize = 2)
hidespines!(ax)
hidexdecorations!(ax)
hideydecorations!(ax)
fig
```

#### Introducing some decoherence

The figure obtained above coulb be already a potential logo for the package. However, we see that the fringe patterns are more intense than the three coherent gaussian amplitudes. We can introduce some decoherence to reduce this effect. Thus, we evolve the system under the evolution of a damped quantum harmonic oscillator, which is described by the Lindblad master equation:

```math
\frac{d \hat{\rho}}{dt} = -i [\hat{H}, \hat{\rho}] + \gamma \left( 2 \hat{a} \hat{\rho} \hat{a}^\dagger - \hat{a}^\dagger \hat{a} \hat{\rho} - \hat{\rho} \hat{a}^\dagger \hat{a} \right)
```

where ``\hat{\rho}`` is the density matrix, ``\hat{H} = \omega \hat{a}^\dagger \hat{a}``is the Hamiltonian of the harmonic oscillator (``\hbar = 1``), ``\hat{a}``and ``\hat{a}^\dagger``are the annihilation and creation operators, and ``\gamma``is the damping rate. Thus, we initialize the system in the triangular cat state and evolve it under the Lindblad master equation, using the [`mesolve`](@ref) function.

```@example logo
γ = 0.012

a = destroy(N)
H = a' * a
c_ops = [sqrt(γ) * a]

tlist = range(0, 2π, 100)

sol = mesolve(H, ψ, tlist, c_ops, progress_bar = false)
nothing # hide
```

And the Wigner function becomes more uniform:

```@example logo
wig = wigner(sol.states[end], xvec, yvec, g = 2)

fig = Figure(size = (500, 500), figure_padding = 0)
ax = Axis(fig[1, 1])

img_wig = heatmap!(ax, xvec, yvec, wig', colormap = :RdBu, interpolate = true, rasterize = 2)
hidespines!(ax)
hidexdecorations!(ax)
hideydecorations!(ax)

fig
```

At this stage, we have finished to use the `QuantumToolbox` package. From now on, we will use the `CairoMakie` package to define custom colormaps and plot the Wigner function in a Julia logo style.

### Custom Colormap

We define a custom colormap that changes depending on the Wigner function and spatial coordinates. Indeed, we want the three different colormaps, in the regions corresponding to the three coherent states, to match the colors of the Julia logo. We also want the colormap change to be smooth, so we use a Gaussian function to blend the colors. We introduce also a Wigner function dependent transparency to make the logo more appealing.

```@example logo
function set_color_julia(x, y, wig::T, α1, α2, α3, cmap1, cmap2, cmap3, δ) where {T}
    amp1 = gaussian(x, real(α1), δ) * gaussian(y, imag(α1), δ)
    amp2 = gaussian(x, real(α2), δ) * gaussian(y, imag(α2), δ)
    amp3 = gaussian(x, real(α3), δ) * gaussian(y, imag(α3), δ)

    c1 = get(cmap1, wig)
    c2 = get(cmap2, wig)
    c3 = get(cmap3, wig)

    c_tot = (amp1 * c1 + amp2 * c2 + amp3 * c3) / (amp1 + amp2 + amp3)

    wig_abs = abs(2 * (wig - 1 / 2))
    # We introduce some non-linearity to increase the contrast
    alpha = 2 * (1 / (1 + exp(-5 * wig_abs)) - 1 / 2)

    return RGBAf(c_tot.r, c_tot.g, c_tot.b, alpha)
end

X, Y = meshgrid(xvec, yvec)
δ = 1.25 # Smoothing parameter for the Gaussian functions
```

#### Colormaps from the Julia colors

We define the colormaps for the three coherent states using the colors of the Julia logo. We use the `cgrad` function from the `CairoMakie` package to create the colormaps.

```@example logo
julia_red = RGBAf(0.796, 0.235, 0.2, 1.0)
julia_green = RGBAf(0.22, 0.596, 0.149, 1.0)
julia_blue = RGBAf(0.251, 0.388, 0.847, 1.0)
julia_purple = RGBAf(0.584, 0.345, 0.698, 1.0)
n_repeats = 2

cmap1 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_green, n_repeats)))
cmap2 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_red, n_repeats)))
cmap3 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_purple, n_repeats)))
```

### Normalizing the Wigner function and applying the custom colormap

The colormaps require the input to be in the range ``[0, 1]``. We normalize the Wigner function such that the maximum value is ``1``and the zeros are set to ``0.5``.

```@example logo
vmax = maximum(wig)
wig_normalized = wig ./ (vmax * 2) .+ 1 / 2
nothing # hide
```

And we now apply this custom colormap to make an image (a `Matrix{RGBAf}`).

```@example logo
img = set_color_julia.(X, Y, wig_normalized, α1, α2, α3, Ref(cmap1), Ref(cmap2), Ref(cmap3), δ)
```

### Final Plot

Finally, we plot the Wigner function with the custom colormap.

```@example logo
fig = Figure(size = (500, 500), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1], backgroundcolor = :transparent)
image!(ax, img')
hidespines!(ax)
hidexdecorations!(ax)
hideydecorations!(ax)
fig
```

## Conclusion

This tutorial demonstrates how to generate the [QuantumToolbox.jl](https://github.com/albertomercurio/QuantumToolbox.jl) logo using the package itself and [Makie.jl](https://github.com/MakieOrg/Makie.jl) for visualization. The logo is a visualization of the Wigner function of a triangular cat state, with a custom colormap that highlights the different coherent states with colors matching the Julia logo.
