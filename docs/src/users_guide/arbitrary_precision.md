# [Arbitrary Precision](@id doc:Arbitrary-Precision)

Double precision (`Float64`) carries about 16 significant decimal digits. For the overwhelming majority of quantum simulations that is far more than enough, and you should keep using it.

But some of the most interesting quantities in quantum physics are *exponentially small*: tunneling splittings, Liouvillian gaps, metastable switching rates, and the lifetimes of symmetry-broken states all scale like ``e^{-S}`` for some action or system size ``S``. Push a little deeper into the parameter regime and the number you are chasing drops below ``10^{-16}`` relative to the energy scales in your Hamiltonian. At that point the answer is no longer physics — it is roundoff.

The failure is quiet, which is what makes it dangerous. The eigensolver still converges, the numbers still look plausible, and the plot still has a curve on it. It is simply the wrong curve. Below we show an example where `Float64` confidently reports a result that is wrong by nearly ten orders of magnitude, and where switching a single type argument fixes it.

`QuantumToolbox.jl` supports arbitrary precision throughout: states, operators, superoperators, the eigensolvers, and the time-evolution solvers are all generic over the number type.

## [Why Julia makes this possible](@id doc:Arbitrary-Precision:why-julia)

This capability is not something `QuantumToolbox.jl` implements by hand — it is inherited from Julia's generic programming model, in the same way that [Automatic Differentiation](@ref doc:autodiff) support was largely inherited rather than built.

Julia compiles the *same generic source code* into specialized native code for whichever number type you hand it. [`destroy(Complex{BigFloat}, N)`](@ref destroy) does not run a separate high-precision code path; it runs the same one, respecialized by the compiler for `Complex{BigFloat}`. Nothing in the library needs to know in advance which types users will care about.

It is worth being precise about why this is hard elsewhere:

- **`Python QuTiP`** cannot do this at all. Its data layer (`Dense` and `CSR`) is compiled Cython with `double complex` hardcoded, so precision is not a knob that exists.
- **`JAX`** cannot either. XLA supports only hardware float types (`bf16`, `f16`, `f32`, `f64`); arbitrary precision has nowhere to compile *to*.
- **`NumPy`** technically can, via `dtype=object` arrays holding `mpmath` scalars, but every elementary operation then becomes an interpreted per-element Python call with no BLAS behind it — a fundamentally different and far slower execution path than the fast one.

The honest caveat is that `BigFloat` is not magic: it wraps the MPFR C library, allocates on the heap, and is genuinely slow. The more interesting part of the story is types like [`Double64`](https://github.com/JuliaMath/DoubleFloats.jl) from `DoubleFloats.jl`, which represents ~32 digits as a pair of `Float64`s. It is written in pure Julia, so it inlines and compiles down to ordinary floating-point instructions with no interpreter overhead at all — while remaining a type the library has never heard of. That combination is what Python and JAX structurally cannot offer.

## [A warm-up: where double precision runs out](@id doc:Arbitrary-Precision:warm-up)

Before touching quantum objects, it helps to see the mechanism in isolation. Consider evaluating ``\sqrt{x+1} - \sqrt{x}`` for large ``x``. For ``x = 10^{16}`` the two square roots agree to ~16 digits, so subtracting them cancels away every digit you had — this is *catastrophic cancellation*. The algebraically identical form

```math
\sqrt{x+1} - \sqrt{x} = \frac{1}{\sqrt{x+1} + \sqrt{x}}
```

never subtracts nearby numbers, and so is stable.

```@example arbitrary_precision
using DoubleFloats

setprecision(BigFloat, 256)

naive(x) = sqrt(x + 1) - sqrt(x)
stable(x) = 1 / (sqrt(x + 1) + sqrt(x))

reference = stable(BigFloat("1e16"))   # 256-bit reference value

for (label, x) in (
    ("Float64 ", 1e16),
    ("Double64", Double64("1e16")),
    ("BigFloat", BigFloat("1e16")),
)
    err(v) = Float64(abs((BigFloat(v) - reference) / reference))
    println("$label   naive: ", rpad(Float64(naive(x)), 24), " (rel. err. ", err(naive(x)), ")")
    println("$label  stable: ", rpad(Float64(stable(x)), 24), " (rel. err. ", err(stable(x)), ")")
end
```

`Float64` with the naive formula returns exactly zero: not an inaccurate answer, but a complete loss of every significant digit. Note the two independent levers here — a better *algorithm* (the stable formula) and a better *number type*. Extra precision is not a substitute for numerical care, and the rest of this page is about the cases where care alone is not enough.

## [Using arbitrary precision in QuantumToolbox](@id doc:Arbitrary-Precision:usage)

Every state and operator generating function takes an optional leading type argument:

```@example arbitrary_precision
using QuantumToolbox

destroy(Complex{BigFloat}, 5)
```

which propagates through arithmetic exactly as you would expect:

```@example arbitrary_precision
normalize(fock(BigFloat, 10, 3) + fock(BigFloat, 10, 4))
```

Some functions instead infer the element type from the *value* you pass, which is often more convenient:

```@example arbitrary_precision
α = BigFloat(1) + 1im * BigFloat(1)
coherent(6, α)   # eltype follows α
```

The working precision of `BigFloat` is global and set with `setprecision` (in bits — 128 bits is ~38 decimal digits, and is usually plenty):

```@example arbitrary_precision
setprecision(BigFloat, 128)
eltype(destroy(Complex{BigFloat}, 5))
```

!!! note "Which packages do I need?"
    Core operations work out of the box, but the eigensolvers need generic linear algebra backends that are not `QuantumToolbox.jl` dependencies. Load them yourself:

    - [`GenericSchur.jl`](https://github.com/RalphAS/GenericSchur.jl) — the generic Schur decomposition used by [`eigenstates`](@ref) and [`eigsolve`](@ref) for non-`BlasFloat` element types.
    - [`Sparspak.jl`](https://github.com/PetrKryslUCSD/Sparspak.jl) — generic sparse LU, needed for the sparse shift-and-invert path (`sparse = Val(true)` with a `sigma`).
    - [`DoubleFloats.jl`](https://github.com/JuliaMath/DoubleFloats.jl) — optional, provides `Double64`.

    Forgetting these produces a `MethodError` from `LinearAlgebra`, not an obvious "unsupported precision" message.

## [Exponentially small tunneling splittings](@id doc:Arbitrary-Precision:double-well)

Now the physics. Take a particle in a quartic double well (in units ``\hbar = m = \omega = 1``),

```math
\hat{H} = \frac{\hat{p}^2}{2} + \lambda \left( \hat{x}^2 - x_0^2 \right)^2 \, ,
```

with minima at ``x = \pm x_0`` separated by a barrier of height ``\lambda x_0^4``. The two lowest eigenstates are the symmetric and antisymmetric combinations of the states localized in each well, and they are split by an energy ``\Delta E = E_1 - E_0`` that is *exponentially small* in the barrier action. That splitting is the entire physics of the problem: a particle prepared in the left well tunnels to the right one in a time ``\pi / \Delta E``.

We build ``\hat{x}`` and ``\hat{p}`` from the annihilation operator, so that the whole Hamiltonian inherits whatever type we pass in:

```@example arbitrary_precision
using LinearAlgebra
using GenericSchur   # generic Schur decomposition, needed for non-Float64 eigensolvers

function double_well(::Type{T}, N, λ, x0) where {T}
    a = destroy(T, N)
    x = (a' + a) / sqrt(T(2))
    p = im * (a' - a) / sqrt(T(2))
    W = x * x - T(x0)^2 * qeye(T, N)
    return p * p / 2 + T(λ) * W * W
end

# splitting between the two lowest levels
function splitting(::Type{T}, N, λ, x0) where {T}
    E = real.(eigenstates(double_well(T, N, λ, x0)).values)
    return E[2] - E[1]
end
nothing # hide
```

Note that `double_well` contains no reference to precision whatsoever — it is ordinary code, written the way you would write it anyway. Now sweep the well separation ``x_0`` (deepening the barrier) at three different precisions:

```@example arbitrary_precision
N = 250     # Fock cutoff, large enough to be converged at the deepest barrier
λ = 0.5
x0_list = range(1.0, 3.5, 11)

Δ64  = [splitting(ComplexF64, N, λ, x0) for x0 in x0_list]
Δd64 = [splitting(Complex{Double64}, N, λ, x0) for x0 in x0_list]

# cross-check the deepest barrier against BigFloat (much slower, so just one point)
Δbig_end = splitting(Complex{BigFloat}, N, λ, x0_list[end])

println("at the deepest barrier, x0 = ", x0_list[end])
println("  Float64  : ", Δ64[end])
println("  Double64 : ", Δd64[end])
println("  BigFloat : ", Δbig_end)
```

`Float64` reports a splitting of ``\sim 10^{-13}`` while the true value is ``\sim 10^{-23}`` — wrong by nearly ten orders of magnitude. `Double64` already agrees with `BigFloat` to several digits here, so the cheap type is enough and we can plot it alone. Sweeping the whole range shows how the failure sets in:

```@example arbitrary_precision
using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())

fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],
    title = "Tunneling splitting of a quartic double well",
    xlabel = L"well separation $x_0$",
    ylabel = L"\Delta E = E_1 - E_0",
    yscale = log10,
)
lines!(ax, x0_list, abs.(Δ64), label = "Float64", linewidth = 2)
lines!(ax, x0_list, abs.(Float64.(Δd64)), label = "Double64", linewidth = 2, linestyle = :dash)
axislegend(ax, position = :lb)

fig
```

The two curves agree perfectly while the splitting is large. Around ``\Delta E \sim 10^{-13}`` the `Float64` curve simply *stops descending* and flattens into a noise floor, wandering non-monotonically from point to point — it is now reporting the accumulated roundoff of the eigensolver rather than a physical energy. `Double64` continues down the straight line that the exponential law predicts.

Read physically, the `Float64` curve claims that beyond ``x_0 \approx 3`` the barrier stops mattering and the tunneling rate saturates. That is not a small quantitative error; it is a qualitatively wrong statement about the system. Tunneling keeps slowing down exponentially, and only the high-precision types can see it.

The same trap appears in open systems. The **Liouvillian gap** of a driven-dissipative resonator — the eigenvalue of ``\mathcal{L}`` with the smallest nonzero real part, which sets the switching time between metastable states — becomes exponentially small near a bistability, and shrinks further with system size. Computed in `Float64` it eventually stops being physics: it flattens into the same roundoff floor, can come out *negative* (implying a density matrix that grows without bound, which no Lindblad master equation permits), and can place the gap minimum at entirely the wrong drive amplitude. The remedy is identical — pass a wider element type to [`liouvillian`](@ref) and [`eigenstates`](@ref).

Note also that `Double64`, not `BigFloat`, was enough here, at roughly a tenth of the cost. Reaching for `BigFloat` on every problem is rarely the right move; it is usually worth trying the cheap type first.

## [The solvers are generic too](@id doc:Arbitrary-Precision:solvers)

The same type argument flows through the time-evolution solvers, so [`sesolve`](@ref) and [`mesolve`](@ref) work unchanged at high precision:

```@example arbitrary_precision
N_solver = 20
a_big = destroy(Complex{BigFloat}, N_solver)
ψ0_big = fock(Complex{BigFloat}, N_solver, 1)
H_big = a_big' * a_big + BigFloat(0.1) / 2 * a_big' * a_big' * a_big * a_big
c_ops_big = [sqrt(BigFloat(0.1)) * a_big]

tlist = range(0, 10, 100)
sol = mesolve(H_big, ψ0_big, tlist, c_ops_big, e_ops = [a_big' * a_big], progress_bar = Val(false))

eltype(sol.expect)
```

Note that the expectation values come back in the element type you asked for, rather than being silently narrowed to `ComplexF64`.

Be aware of what this does and does not buy you, though. Extra precision does **not** generally make a time-evolution result better: the error of an ODE integration is dominated by the solver tolerances, not by roundoff, and tightening `abstol`/`reltol` in `Float64` is both cheaper and more effective. High-precision solvers earn their cost when the evolution feeds something else that is precision-critical — for example [`eigsolve_al`](@ref), which integrates the master equation to extract Liouvillian eigenvalues, and so inherits exactly the sensitivity described above.

## [Caveats and performance](@id doc:Arbitrary-Precision:caveats)

!!! warning "Arbitrary precision is not free"
    - **It is slow.** `BigFloat` allocates on the heap and calls into MPFR for every operation; expect one to two orders of magnitude slowdown, and no multithreaded BLAS. `Double64` is typically only a few times slower than `Float64` and should be your first try — in the double-well sweep above it is roughly ten times faster than `BigFloat` and gives an identical answer.
    - **No GPU support.** The CUDA extension deliberately supports only `Float32`/`Float64` word sizes; arbitrary precision is CPU-only.
    - **Convergence is not just about precision.** In the double-well example the Fock cutoff `N` must be large enough that the *truncated* problem is converged. At `N = 150`, even `BigFloat` returns a wrong splitting — the arithmetic is exact but the basis is not. Increasing precision cannot fix a basis-truncation error, so check both.
    - **Tolerances do not scale automatically.** [`eigsolve`](@ref)'s default `tol = 1e-8` is a hardcoded `Float64` literal. If you want more digits than that out of an iterative solver, pass a tighter `tol` explicitly.
    - **Tested surface.** Arbitrary precision is covered by the test suite for [`sesolve`](@ref), [`mesolve`](@ref), [`eigenstates`](@ref), and [`eigsolve_al`](@ref). Other routines such as [`steadystate`](@ref), [`mcsolve`](@ref), and [`spectrum`](@ref) are written generically and may well work, but are not currently tested at high precision — if you rely on them, validate against a `Float64` result in a regime where both are trustworthy.
