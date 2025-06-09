# [Plotting on the Bloch Sphere](@id doc:Plotting-on-the-Bloch-Sphere)

```@setup Bloch_sphere_rendering
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## Introduction

When studying the dynamics of a two-level system, it's often convenient to visualize the state of the system by plotting the state vector or density matrix on the Bloch sphere.

In [`QuantumToolbox`](https://qutip.org/QuantumToolbox.jl/), this can be done using the [`Bloch`](@ref) or [`plot_bloch`](@ref) methods that provide same syntax as [QuTiP](https://qutip.readthedocs.io/en/stable/guide/guide-bloch.html).

## Create a Bloch Sphere

In [`QuantumToolbox`](https://qutip.org/QuantumToolbox.jl/), creating a [`Bloch`](@ref) sphere is accomplished by calling either:

!!! note "Import plotting libraries"
    Remember to import plotting libraries first. Here, we demonstrate the functionalities with [`CairoMakie.jl`](https://docs.makie.org/stable/explanations/backends/cairomakie.html).

```@example Bloch_sphere_rendering
b = Bloch()
```

which will load an instance of [`Bloch`](@ref). Before getting into the details of these objects, we can simply plot the blank [`Bloch`](@ref) sphere associated with these instances via:

```@example Bloch_sphere_rendering
fig, _ = render(b)
fig
```

See the [API documentation for Bloch sphere](@ref doc-API:Bloch-Sphere) for a full list of other available functions.

## Add a single data point

As an example, we can add a single data point via [`add_points!`](@ref):

```@example Bloch_sphere_rendering
pnt = [1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)]
add_points!(b, pnt)
fig, _ = render(b)
fig
```

## Add a single vector

Add a single vector via  [`add_vectors!`](@ref):

```@example Bloch_sphere_rendering
vec = [0, 1, 0]
add_vectors!(b, vec)
fig, _ = render(b)
fig
```

## Add a single quantum state

Add another vector corresponding to the ``|0\rangle`` state:

```@example Bloch_sphere_rendering
z0 = basis(2, 0)
add_states!(b, z0)
fig, _ = render(b)
fig
```

## Add multiple data

We can also plot multiple points, vectors, and states at the same time by passing arrays instead of individual elements via [`add_points!`](@ref), [`add_vectors!`](@ref), and [`add_states!`](@ref), respectively. Before giving an example, we can use [`clear!`](@ref) to remove the current data from our [`Bloch`](@ref) sphere instead of creating a new instance:

```@example Bloch_sphere_rendering
clear!(b)
fig, _ = render(b)
fig
```

Now on the same [`Bloch`](@ref) sphere, we can plot the three states via [`add_states!`](@ref) associated with the `x`, `y`, and `z` directions:

```@example Bloch_sphere_rendering
x = basis(2, 0) + basis(2, 1)
y = basis(2, 0) + im * basis(2, 1)
z = basis(2, 0)
add_states!(b, [x, y, z])
fig, _ = render(b)
fig
```

!!! note "State normalization"
    The function [`add_states!`](@ref) will automatically normalize the given quantum state(s), while [`add_vectors!`](@ref) does not normalize the given vectors.

A similar method works for adding vectors:

```@example Bloch_sphere_rendering
clear!(b)
vecs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
add_vectors!(b, vecs)
fig, _ = render(b)
fig
```

# Add lines and arcs

You can also add lines and arcs via [`add_line!`](@ref) and [`add_arc!`](@ref) respectively:

```@example Bloch_sphere_rendering
add_line!(b, x, y)
add_arc!(b, y, z)
fig, _ = render(b)
fig
```

## Add multiple points

Adding multiple points to the [`Bloch`](@ref) sphere works slightly differently than adding multiple states or vectors. For example, lets add a set of `20` points around the equator (after calling [`clear!`](@ref)):

```@example Bloch_sphere_rendering
clear!(b)

th = LinRange(0, 2π, 20)
xp = cos.(th)
yp = sin.(th)
zp = zeros(20)
pnts = [xp, yp, zp]
add_points!(b, pnts)
fig, lscene = render(b)
fig
```

Notice that, in contrast to states or vectors, each point remains the same color as the initial point. This is because adding multiple data points using [`add_points!`](@ref) is interpreted, by default, to correspond to a single data point (single qubit state) plotted at different times. This is very useful when visualizing the dynamics of a qubit. If we want to plot additional qubit states we can call additional [`add_points!`](@ref) function:

```@example Bloch_sphere_rendering
xz = zeros(20)
yz = sin.(th)
zz = cos.(th)
add_points!(b, [xz, yz, zz])
fig, lscene = render(b)
fig
```

The color and shape of the data points is varied automatically by [`Bloch`](@ref). Notice how the color and point markers change for each set of data. Again, we have had to call [`add_points!`](@ref) twice because adding more than one set of multiple data points is not supported by the [`add_points!`](@ref) function.

What if we want to vary the color of our points. We can tell [`Bloch`](@ref) to vary the color of each point according to the colors listed in the `point_color` field (see [Configuring the Bloch sphere](@ref doc:Configuring-the-Bloch-sphere) below). Again, after [`clear!`](@ref):

```@example Bloch_sphere_rendering
clear!(b)

xp = cos.(th)
yp = sin.(th)
zp = zeros(20)
pnts = [xp, yp, zp]
add_points!(b, pnts, meth=:m) # add `meth=:m` to signify 'multi' colored points
fig, lscene = render(b)
fig
```

Now, the data points cycle through a variety of predefined colors. Now lets add another set of points, but this time we want the set to be a single color, representing say a qubit going from the ``|0\rangle`` state to the ``|1\rangle`` state in the `y-z` plane:

```@example Bloch_sphere_rendering
pnts = [xz, yz, zz]
add_points!(b, pnts) # no `meth=:m`
fig, lscene = render(b)
fig
```

## [Configuring the Bloch sphere](@id doc:Configuring-the-Bloch-sphere)

At the end of the last section we saw that the colors and marker shapes of the data plotted on the Bloch sphere are automatically varied according to the number of points and vectors added. But what if you want a different choice of color, or you want your sphere to be purple with different axes labels? Well then you are in luck as the Bloch class has many attributes which one can control. Assuming `b = Bloch()`:

| **Field** | **Description** | **Default setting** |
|:----------|:----------------|:--------------------|
| `b.points` | Points to plot on the Bloch sphere (3D coordinates) | `Vector{Matrix{Float64}}()` (empty) |
| `b.vectors` | Vectors to plot on the Bloch sphere | `Vector{Vector{Float64}}()` (empty) |
| `b.lines` | Lines to draw on the sphere (points, style, properties) | `Vector{Tuple{Vector{Vector{Float64}},String}}()` (empty) |
| `b.arcs` | Arcs to draw on the sphere | `Vector{Vector{Vector{Float64}}}()` (empty) |
| `b.font_color` | Color of axis labels and text | `"black"` |
| `b.font_size` | Font size for labels | `20` |
| `b.frame_alpha` | Transparency of the frame background | `0.1` |
| `b.frame_color` | Background color of the frame | `"gray"` |
| `b.frame_limit` | Axis limits for the 3D frame (symmetric around origin) | `1.2` |
| `b.point_default_color` | Default color cycle for points | `["blue", "red", "green", "#CC6600"]` |
| `b.point_color` | List of colors for Bloch point markers to cycle through | `Union{Nothing,String}[]` |
| `b.point_marker` | List of point marker shapes to cycle through | `[:circle, :rect, :diamond, :utriangle]` |
| `b.point_size` | List of point marker sizes (not all markers look the same size when plotted) | `[5.5, 6.2, 6.5, 7.5]` |
| `b.point_style` | List of marker styles | `Symbol[]` |
| `b.point_alpha` | List of marker transparencies | `Float64[]` |
| `b.sphere_color` | Color of Bloch sphere surface | `0.2` |
| `b.sphere_alpha` | Transparency of sphere surface | `"#FFDDDD"` |
| `b.vector_color` | Colors for vectors | `["green", "#CC6600", "blue", "red"]` |
| `b.vector_width` | Width of vectors | `0.025` |
| `b.vector_arrowsize` | Arrow size parameters as (head length, head width, stem width) | `[0.07, 0.08, 0.08]` |
| `b.view` | Azimuthal and elevation viewing angles in degrees | `[30, 30]` |
| `b.xlabel` | Labels for x-axis | `[L"x", ""]` (``+x`` and ``-x``) |
| `b.xlpos` | Positions of x-axis labels | `[1.0, -1.0]` |
| `b.ylabel` | Labels for y-axis | `[L"y", ""]` (``+y`` and ``-y``) |
| `b.ylpos` | Positions of y-axis labels | `[1.0, -1.0]` |
| `b.zlabel` | Labels for z-axis | `[L"\|0\rangle", L"\|1\rangle]"` (``+z`` and ``-z``) |
| `b.zlpos` | Positions of z-axis labels | `[1.0, -1.0]` |

These properties can also be accessed via the `print` command:

```@example Bloch_sphere_rendering
b = Bloch()
print(b)
```