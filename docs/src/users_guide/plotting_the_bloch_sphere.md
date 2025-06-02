# [Plotting on the Bloch Sphere](@id doc: plotting_the_bloch_sphere.md)

```@setup Bloch_sphere_rendering
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## [Introduction](@id doc:Bloch_sphere_rendering)

When studying the dynamics of a two-level system, it's often convenient to visualize the state of the system by plotting the state vector or density matrix on the Bloch sphere.

In [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/), this can be done using the [`Bloch`](@ref) or [`plot_bloch`](@ref) methods that provide same syntax as [QuTiP](https://qutip.readthedocs.io/en/stable/guide/guide-bloch.html).

## Create a Bloch Sphere

```@example Bloch_sphere_rendering
b = Bloch();
fig, _ = render(b)
```

## Add a Single Data Point

```@example Bloch_sphere_rendering
pnt = [1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)]
add_points!(b, pnt)
fig, _ = render(b)
fig
```

## Add a Single Vector

```@example Bloch_sphere_rendering
clear!(b)
vec = [0, 1, 0];
add_vectors!(b, vec)
fig, _ = render(b)
fig
```

## Add Multiple Vectors

```@example Bloch_sphere_rendering
clear!(b)
vecs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
add_vectors!(b, vecs)
fig, _ = render(b)
fig
```

# Add Arc, Line, and Vector

```@example Bloch_sphere_rendering
clear!(b)
vec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
add_vectors!(b, vec);
add_line!(b, [1,0,0], [0,1,0])
add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
fig, _ = render(b)
fig
```

# Add Quantum States

```@example Bloch_sphere_rendering
clear!(b)
x = basis(2, 0) + basis(2, 1)
y = basis(2, 0) - im * basis(2, 1)
z = basis(2, 0)
b = Bloch()
add_states!(b, [x, y, z])
fig, _ = render(b)
fig
```

## Add Multiple Points Around the Equator

```@example Bloch_sphere_rendering
th = range(0, 2π; length=20);
clear!(b)
xp = cos.(th);
yp = sin.(th);
zp = zeros(20);
pnts = [xp, yp, zp] ;
pnts = Matrix(hcat(xp, yp, zp)');
add_points!(b, pnts);
fig, ax = render(b);
fig
```

## Add Another Set of Points

```@example Bloch_sphere_rendering
xz = zeros(20);
yz = sin.(th);
zz = cos.(th);
points = hcat(xz, yz, zz)';
add_points!(b, Matrix(points), meth=:s, color="orange");
fig, ax = render(b);
fig
```
