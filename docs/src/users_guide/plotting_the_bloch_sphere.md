# [Plotting on the Bloch Sphere](@id doc:Plotting-on-the-Bloch-Sphere)

```@setup Bloch_sphere_rendering
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## [Introduction](@id doc:Bloch_sphere_rendering)

When studying the dynamics of a two-level system, it's often convenient to visualize the state of the system by plotting the state vector or density matrix on the Bloch sphere.

In [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/), this can be done using the [`Bloch`](@ref) or [`plot_bloch`](@ref) methods that provide same syntax as [QuTiP](https://qutip.readthedocs.io/en/stable/guide/guide-bloch.html).

## Create a Bloch Sphere

In [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/), creating a [`Bloch`](@ref) sphere is accomplished by calling either:

```@example Bloch_sphere_rendering
b = Bloch();
```

which will load an instance of [`Bloch`](@ref). Before getting into the details of these objects, we can simply plot the blank [`Bloch`](@ref) sphere associated with these instances via:

```@example Bloch_sphere_rendering
fig, _ = render(b);
fig
```

## Add a Single Data Point

As an example, we can add a single data point via [`add_points!`](@ref):

```@example Bloch_sphere_rendering
pnt = [1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)];
add_points!(b, pnt);
fig, _ = render(b);
fig
```

## Add a Single Vector

and then a single vector via  [`add_vectors!`](@ref):

```@example Bloch_sphere_rendering
vec = [0, 1, 0];
add_vectors!(b, vec)
fig, _ = render(b)
fig
```

and then add another vector corresponding to the ``|0\rangle`` state:

```@example Bloch_sphere_rendering
x = basis(2, 0)
add_states!(b, [x])
fig, _ = render(b)
fig
```

## Add Multiple Vectors

We can also plot multiple points, vectors, and states at the same time by passing arrays instead of individual elements via  [`add_vectors!](@ref). Before giving an example, we can use [`clear!`](@ref) to remove the current data from our [`Bloch`](@ref) sphere instead of creating a new instance:

```@example Bloch_sphere_rendering
clear!(b)
fig, _ = render(b)
fig
```

Now on the same [`Bloch`](@ref) sphere, we can plot the three states via [`add_states!`](@ref) associated with the `x`, `y`, and `z` directions:

```@example Bloch_sphere_rendering
x = basis(2, 0) + basis(2, 1)
y = basis(2, 0) - im * basis(2, 1)
z = basis(2, 0)
b = Bloch()
add_states!(b, [x, y, z])
fig, _ = render(b)
fig
```

A similar method works for adding vectors:

```@example Bloch_sphere_rendering
clear!(b)
vecs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
add_vectors!(b, vecs)
fig, _ = render(b)
fig
```

# Add Arc, Line, and Vector

You can also add lines and arcs via [`add_line!`](@ref) and [`add_arc!`](@ref) respectively:

```@example Bloch_sphere_rendering
clear!(b)
vec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
add_vectors!(b, vec);
add_line!(b, [1,0,0], [0,1,0])
add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
fig, _ = render(b)
fig
```

## Add Multiple Points

Adding multiple points to the [`Bloch`](@ref) sphere works slightly differently than adding multiple states or vectors. For example, lets add a set of `20` points around the equator (after calling [`clear!`](@ref)):

```@example Bloch_sphere_rendering
th = range(0, 2π; length=20);
clear!(b)
xp = cos.(th);
yp = sin.(th);
zp = zeros(20);
pnts = [xp, yp, zp];
add_points!(b, pnts);
fig, ax = render(b);
fig
```

Notice that, in contrast to states or vectors, each point remains the same color as the initial point. This is because adding multiple data points using [`add_points!`](@ref) is interpreted, by default, to correspond to a single data point (single qubit state) plotted at different times. This is very useful when visualizing the dynamics of a qubit. If we want to plot additional qubit states we can call additional [`add_points!`](@ref):

## Add Another Set of Points

```@example Bloch_sphere_rendering
xz = zeros(20);
yz = sin.(th);
zz = cos.(th);
pnts = [xz, yz, zz];
add_points!(b, pnts);
fig, ax = render(b);
fig
```

The color and shape of the data points is varied automatically by [`Bloch`](@ref). Notice how the color and point markers change for each set of data.

What if we want to vary the color of our points. We can tell [`Bloch`](@ref) to vary the color of each point according to the colors listed in the `point_color` attribute.

```@example Bloch_sphere_rendering
clear!(b)
xp = cos.(th);
yp = sin.(th);
zp = zeros(20);
pnts = [xp, yp, zp];
add_points!(b, pnts, meth=:m);
fig, ax = render(b);
fig
```

Now, the data points cycle through a variety of predefined colors. Now lets add another set of points, but this time we want the set to be a single color, representing say a qubit going from the ``|0\rangle`` state to the ``|1\rangle`` state in the `y-z` plane:

```@example Bloch_sphere_rendering
pnts = [xz, yz, zz] ;
add_points!(b, pnts);
fig, ax = render(b);
fig
```
