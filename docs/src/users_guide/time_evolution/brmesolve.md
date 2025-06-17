# [Bloch-Redfield master equation](@id doc-TE:Bloch-Redfield-master-equation)

## Introduction

The Lindblad master equation introduced earlier is constructed so that it describes a physical evolution of the density matrix (i.e., trace and positivity preserving), but it does not provide a connection to any underlying microscopic physical model.
The Lindblad operators (collapse operators) describe phenomenological processes, such as for example dephasing and spin flips, and the rates of these processes are arbitrary parameters in the model.
In many situations the collapse operators and their corresponding rates have clear physical interpretation, such as dephasing and relaxation rates, and in those cases the Lindblad master equation is usually the method of choice.

However, in some cases, for example systems with varying energy biases and eigenstates and that couple to an environment in some well-defined manner (through a physically motivated system-environment interaction operator), it is often desirable to derive the master equation from more fundamental physical principles, and relate it to for example the noise-power spectrum of the environment.

The Bloch-Redfield formalism is one such approach to derive a master equation from a microscopic system.
It starts from a combined system-environment perspective, and derives a perturbative master equation for the system alone, under the assumption of weak system-environment coupling.
One advantage of this approach is that the dissipation processes and rates are obtained directly from the properties of the environment.
On the downside, it does not intrinsically guarantee that the resulting master equation unconditionally preserves the physical properties of the density matrix (because it is a perturbative method).
The Bloch-Redfield master equation must therefore be used with care, and the assumptions made in the derivation must be honored.
(The Lindblad master equation is in a sense more robust -- it always results in a physical density matrix -- although some collapse operators might not be physically justified).
For a full derivation of the Bloch Redfield master equation, see e.g. [Cohen_Tannoudji_atomphoton](@citet) or [breuer2002](@citet).
Here we present only a brief version of the derivation, with the intention of introducing the notation and how it relates to the implementation in `QuantumToolbox.jl`.


## Brief Derivation and Definitions

The starting point of the Bloch-Redfield formalism is the total Hamiltonian for the system and the environment (bath): ``H = H_{\rm S} + H_{\rm B} + H_{\rm I}``, where ``H`` is the total system+bath Hamiltonian, ``H_{\rm S}`` and ``H_{\rm B}`` are the system and bath Hamiltonians, respectively, and ``H_{\rm I}`` is the interaction Hamiltonian.

The most general form of a master equation for the system dynamics is obtained by tracing out the bath from the von-Neumann equation of motion for the combined system (``\dot\rho = -i\hbar^{-1}[H, \rho]``). In the interaction picture the result is

```math
    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^t d\tau\;  {\rm Tr}_B [H_I(t), [H_I(\tau), \rho_S(\tau)\otimes\rho_B]],
```

where the additional assumption that the total system-bath density matrix can be factorized as ``\rho(t) \approx \rho_S(t) \otimes \rho_B``.
This assumption is known as the Born approximation, and it implies that there never is any entanglement between the system and the bath, neither in the initial state nor at any time during the evolution.
*It is justified for weak system-bath interaction.*

The master equation above is non-Markovian, i.e., the change in the density matrix at a time ``t`` depends on states at all times ``\tau < t``, making it intractable to solve both theoretically and numerically.
To make progress towards a manageable master equation, we now introduce the Markovian approximation, in which ``\rho_S(\tau)`` is replaced by ``\rho_S(t)``.
The result is the Redfield equation

```math
    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^t d\tau\; {\rm Tr}_B [H_I(t), [H_I(\tau), \rho_S(t)\otimes\rho_B]],
```

which is local in time with respect the density matrix, but still not Markovian since it contains an implicit dependence on the initial state. By extending the integration to infinity and substituting ``\tau \rightarrow t-\tau``, a fully Markovian master equation is obtained:

```math
    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^\infty d\tau\; {\rm Tr}_B [H_I(t), [H_I(t-\tau), \rho_S(t)\otimes\rho_B]].
```

The two Markovian approximations introduced above are valid if the time-scale with which the system dynamics changes is large compared to the time-scale with which correlations in the bath decays (corresponding to a "short-memory" bath, which results in Markovian system dynamics).

The Markovian master equation above is still on a too general form to be suitable for numerical implementation. We therefore assume that the system-bath interaction takes the form ``H_I = \sum_\alpha A_\alpha \otimes B_\alpha`` and where ``A_\alpha`` are system operators and ``B_\alpha`` are bath operators.
This allows us to write master equation in terms of system operators and bath correlation functions:

```math
    \frac{d}{dt}\rho_S(t) =
    -\hbar^{-2}
    \sum_{\alpha\beta}
    \int_0^\infty d\tau\;
    \left\{
    g_{\alpha\beta}(\tau) \left[A_\alpha(t)A_\beta(t-\tau)\rho_S(t) - A_\beta(t-\tau)\rho_S(t)A_\alpha(t)\right]
    \right. \nonumber\\
    \left.
    g_{\alpha\beta}(-\tau) \left[\rho_S(t)A_\alpha(t-\tau)A_\beta(t) - A_\beta(t)\rho_S(t)A_\beta(t-\alpha)\right]
    \right\},
```

where ``g_{\alpha\beta}(\tau) = {\rm Tr}_B\left[B_\alpha(t)B_\beta(t-\tau)\rho_B\right] = \left<B_\alpha(\tau)B_\beta(0)\right>``,
since the bath state ``\rho_B`` is a steady state.

In the eigenbasis of the system Hamiltonian, where ``A_{mn}(t) = A_{mn} e^{i\omega_{mn}t}``,
``\omega_{mn} = \omega_m - \omega_n`` and ``\omega_m`` are the eigenfrequencies
corresponding the eigenstate ``\left|m\right>``, we obtain in matrix form in the Schrödinger picture

```math
\begin{aligned}
    \frac{d}{dt} \rho_{ab}(t) = & -i\omega_{ab}\rho_{ab}(t)\\
    &-\hbar^{-2} \sum_{\alpha,\beta} \sum_{c,d}^{\rm sec} \int_0^\infty d\tau\;
    \left\{
    g_{\alpha\beta}(\tau)
    \left[\delta_{bd}\sum_nA^\alpha_{an}A^\beta_{nc}e^{i\omega_{cn}\tau}
    -
    A^\beta_{ac} A^\alpha_{db} e^{i\omega_{ca}\tau}
    \right]
    \right. \\
    &+
    \left.
    g_{\alpha\beta}(-\tau)
    \left[\delta_{ac}\sum_n A^\alpha_{dn}A^\beta_{nb} e^{i\omega_{nd}\tau}
    -
    A^\beta_{ac}A^\alpha_{db}e^{i\omega_{bd}\tau}
    \right]
    \right\} \rho_{cd}(t),
\end{aligned}
```

where the "sec" above the summation symbol indicate summation of the secular terms which satisfy ``|\omega_{ab}-\omega_{cd}| \ll \tau_ {\rm decay}``.
This is an almost-useful form of the master equation. The final step before arriving at the form of the Bloch-Redfield master equation that is implemented in QuTiP, involves rewriting the bath correlation function ``g(\tau)`` in terms of the noise-power spectrum of the environment ``S(\omega) = \int_{-\infty}^\infty d\tau e^{i\omega\tau} g(\tau)``:

```math
    \int_0^\infty d\tau\; g_{\alpha\beta}(\tau) e^{i\omega\tau} = \frac{1}{2}S_{\alpha\beta}(\omega) + i\lambda_{\alpha\beta}(\omega),
```

where ``\lambda_{ab}(\omega)`` is an energy shift that is neglected here. The final form of the Bloch-Redfield master equation is


```math
    \frac{d}{dt}\rho_{ab}(t)
    =
    -i\omega_{ab}\rho_{ab}(t)
    +
    \sum_{c,d}^{\rm sec}R_{abcd}\rho_{cd}(t),
```

where

```math
    R_{abcd} =  -\frac{\hbar^{-2}}{2} \sum_{\alpha,\beta}
    \left\{
    \delta_{bd}\sum_nA^\alpha_{an}A^\beta_{nc}S_{\alpha\beta}(\omega_{cn})
    -
    A^\beta_{ac} A^\alpha_{db} S_{\alpha\beta}(\omega_{ca})
    \right. \nonumber\\
    +
    \left.
    \delta_{ac}\sum_n A^\alpha_{dn}A^\beta_{nb} S_{\alpha\beta}(\omega_{dn})
    -
    A^\beta_{ac}A^\alpha_{db} S_{\alpha\beta}(\omega_{db})
    \right\},
```

is the Bloch-Redfield tensor.

The Bloch-Redfield master equation in this form is suitable for numerical implementation. The input parameters are the system Hamiltonian ``H``, the system operators through which the environment couples to the system ``A_\alpha``, and the noise-power spectrum ``S_{\alpha\beta}(\omega)`` associated with each system-environment interaction term.

To simplify the numerical implementation we often assume that ``A_\alpha`` are Hermitian and that cross-correlations between different environment operators vanish, resulting in

```math
    R_{abcd} =  -\frac{\hbar^{-2}}{2} \sum_{\alpha}
    \left\{
    \delta_{bd}\sum_nA^\alpha_{an}A^\alpha_{nc}S_{\alpha}(\omega_{cn})
    -
    A^\alpha_{ac} A^\alpha_{db} S_{\alpha}(\omega_{ca})
    \right. \nonumber\\
    +
    \left.
    \delta_{ac}\sum_n A^\alpha_{dn}A^\alpha_{nb} S_{\alpha}(\omega_{dn})
    -
    A^\alpha_{ac}A^\alpha_{db} S_{\alpha}(\omega_{db})
    \right\}.
```

## Bloch-Redfield master equation in `QuantumToolbox.jl`

### Preparing the Bloch-Redfield tensor

In `QuantumToolbox.jl`, the Bloch-redfield master equation can be calculated using the function [`bloch_redfield_tensor`](@ref). 
It takes two mandatory arguments: The system Hamiltonian ``H`` and a nested list `a_ops` consist of tuples as `(A, spec)` 
with `A::QuantumObject{Operator}` being the operator ``A_\alpha`` and `spec::Function` being the spectral density functions ``S_\alpha(\omega)``.

It is possible to also get the term for a bath `\alpha` directly using [`brterm`](@ref). This function takes only one hermitian coupling operator ``A_\alpha`` and spectral response function ``S_\alpha(\omega)``.

To illustrate how to calculate the Bloch-Redfield tensor, let's consider a two-level atom

```math
    H = -\frac{1}{2}\Delta\sigma_x - \frac{1}{2}\varepsilon_0\sigma_z
```

```@example brmesolve
Δ  = 0.2 * 2π
ε0 = 1.0 * 2π
γ1 = 0.5

H = -Δ/2.0 * sigmax() - ε0/2 * sigmaz()

function ohmic_spectrum(ω)
    if ω == 0.   # dephasing noise
        return γ1
    elseif ω > 0 # relaxing noise  
        return γ1 / 2 * ω / (2π)
    else
        return o
    end
end

R, U = bloch_redfield_tensor(H, [(sigmax(), ohmic_spectrum)])

R
```

Note that it is also possible to add Lindblad dissipation superoperators in the
Bloch-Refield tensor by passing the operators via the third argument `c_ops`
like you would in the [`mesolve`](@ref) or [`mcsolve`](@ref) functions.
For convenience, when the keyword argument `fock_basis = false`, the function [`bloch_redfield_tensor`](@ref) also returns the basis
transformation operator `U`, the eigen vector matrix, since they are calculated in the
process of calculating the Bloch-Redfield tensor `R`, and the `U` are usually
needed again later when transforming operators between the laboratory basis and the eigen basis.
The tensor can be obtained in the laboratory basis by setting `fock_basis = true`,
in that case, the transformation operator is not returned.

### Time evolution

The evolution of a wave function or density matrix, according to the Bloch-Redfield
master equation, can be calculated using the function [`mesolve`](@ref)
using Bloch-Refield tensor in the laboratory basis instead of a liouvillian.
For example, to evaluate the expectation values of the ``\sigma_x``,
``\sigma_y``, and ``\sigma_z`` operators for the example above, we can use the following code:

```@example brmesolve
Δ = 0.2 * 2 * π
ϵ0 = 1.0 * 2 * π
γ1 = 0.5

H = - Δ/2.0 * sigmax() - ϵ0/2.0 * sigmaz()

ohmic_spectrum(ω) = (ω == 0.0) ? γ1 : γ1 / 2 * (ω / (2 * π)) * (ω > 0.0)

a_ops = ((sigmax(), ohmic_spectrum),)
e_ops = [sigmax(), sigmay(), sigmaz()]

# same initial random ket state in QuTiP doc 
ψ0 = Qobj([
    0.05014193+0.66000276im,
    0.67231376+0.33147603im
])

tlist = LinRange(0, 15.0, 1000)
sol = brmesolve(H, ψ0, tlist, a_ops, e_ops=e_ops)
expt_list = real(sol.expect)

# plot the evolution of state on Bloch sphere
sphere = Bloch()
add_points!(sphere, [expt_list[1,:], expt_list[2,:], expt_list[3,:]])
sphere.vector_color = ["red"]

add_vectors!(sphere, [Δ, 0, ϵ0] / √(Δ^2 + ϵ0^2))

fig, _ = render(sphere)
fig
```

The two steps of calculating the Bloch-Redfield tensor and evolving according to
the corresponding master equation can be combined into one by using the function
[`brmesolve`](@ref), in addition to the same arguments as [`mesolve`](@ref) and
[`mcsolve`](@ref), the nested list of operator-spectrum
tuple should be given under `a_ops`.

```@example brmesolve
    output = brmesolve(H, psi0, tlist, a_ops=[[sigmax(),ohmic_spectrum]], e_ops=e_ops)
```

The resulting `output` is of the `struct` [`TimeEvolutionSol`](@ref) like as from [`mesolve`](@ref).
