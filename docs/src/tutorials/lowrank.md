# [Low rank master equation](@id doc-tutor:Low-rank-master-equation)

In this tutorial, we will show how to solve the master equation using the low-rank method. For a detailed explanation of the method, we recommend to read the article [gravina2024adaptive](@cite).

As a test, we will consider the dissipative Ising model with a transverse field. The Hamiltonian is given by

```math
\hat{H} = \frac{J_x}{2} \sum_{\langle i,j \rangle} \sigma_i^x \sigma_j^x + \frac{J_y}{2} \sum_{\langle i,j \rangle} \sigma_i^y \sigma_j^y + \frac{J_z}{2} \sum_{\langle i,j \rangle} \sigma_i^z \sigma_j^z - \sum_i h_i \sigma_i^z + h_x \sum_i \sigma_i^x + h_y \sum_i \sigma_i^y + h_z \sum_i \sigma_i^z,
```

where the sums are over nearest neighbors, and the collapse operators are given by 

```math
c_i = \sqrt{\gamma} \sigma_i^-.
```

We start by importing the packages

```@example lowrank
using QuantumToolbox
using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

Define lattice

```@example lowrank
Nx, Ny = 2, 3
latt = Lattice(Nx = Nx, Ny = Ny)
```

Define lr-space dimensions

```@example lowrank
N_cut = 2         # Number of states of each mode
N_modes = latt.N  # Number of modes
N = N_cut^N_modes # Total number of states
M = latt.N + 1       # Number of states in the LR basis
```

Define lr states. Take as initial state all spins up. All other N states are taken as those with miniman Hamming distance to the initial state.

```@example lowrank
ϕ = Vector{QuantumObject{Vector{ComplexF64},KetQuantumObject,M-1}}(undef, M)
ϕ[1] = kron(fill(basis(2, 1), N_modes)...)

i = 1
for j in 1:N_modes
    global i += 1
    i <= M && (ϕ[i] = SingleSiteOperator(sigmap(), j, latt) * ϕ[1])
end
for k in 1:N_modes-1
    for l in k+1:N_modes
        global i += 1
        i <= M && (ϕ[i] = SingleSiteOperator(sigmap(), k, latt) * SingleSiteOperator(sigmap(), l, latt) * ϕ[1])
    end
end
for i in i+1:M
    ϕ[i] = QuantumObject(rand(ComplexF64, size(ϕ[1])[1]), dims = ϕ[1].dims)
    normalize!(ϕ[i])
end
nothing # hide
```

Define the initial state

```@example lowrank
z = hcat(get_data.(ϕ)...)
B = Matrix(Diagonal([1 + 0im; zeros(M - 1)]))
S = z' * z # Overlap matrix
B = B / tr(S * B) # Normalize B

ρ = QuantumObject(z * B * z', dims = ntuple(i->N_cut, Val(N_modes))); # Full density matrix
```

Define the Hamiltonian and collapse operators

```@example lowrank
# Define Hamiltonian and collapse operators
Jx = 0.9
Jy = 1.04
Jz = 1.0
hx = 0.0
hy = 0.0
hz = 0.0
γ = 1

Sx = mapreduce(i->SingleSiteOperator(sigmax(), i, latt), +, 1:latt.N)
Sy = mapreduce(i->SingleSiteOperator(sigmay(), i, latt), +, 1:latt.N)
Sz = mapreduce(i->SingleSiteOperator(sigmaz(), i, latt), +, 1:latt.N)

H, c_ops = DissipativeIsing(Jx, Jy, Jz, hx, hy, hz, γ, latt; boundary_condition = Val(:periodic_bc), order = 1)
e_ops = (Sx, Sy, Sz)

tl = range(0, 10, 100)
nothing # hide
```

### Full evolution

```@example lowrank
sol_me = mesolve(H, ρ, tl, c_ops; e_ops = [e_ops...]);
Strue = entropy_vn(sol_me.states[end], base=2) / latt.N
```

### Low Rank evolution

Define functions to be evaluated during the low-rank evolution

```@example lowrank
function f_purity(p, z, B)
    N = p.N
    M = p.M
    S = p.S
    T = p.temp_MM

    mul!(T, S, B)
    return tr(T^2)
end

function f_trace(p, z, B)
    N = p.N
    M = p.M
    S = p.S
    T = p.temp_MM

    mul!(T, S, B)
    return tr(T)
end

function f_entropy(p, z, B)
    C = p.A0
    σ = p.Bi

    mul!(C, z, sqrt(B))
    mul!(σ, C', C)
    return entropy_vn(Qobj(Hermitian(σ), type=Operator), base=2)
end
```

Define the options for the low-rank evolution

```@example lowrank
opt = (err_max = 1e-3, p0 = 0.0, atol_inv = 1e-6, adj_condition = "variational", Δt = 0.0);

sol_lr = lr_mesolve(H, z, B, tl, c_ops; e_ops = e_ops, f_ops = (f_purity, f_entropy, f_trace), opt = opt);
```

Plot the results

```@example lowrank
m_me = real(sol_me.expect[3, :]) / Nx / Ny
m_lr = real(sol_lr.expect[3, :]) / Nx / Ny

fig = Figure(size = (500, 350), fontsize = 15)
ax = Axis(fig[1, 1], xlabel = L"\gamma t", ylabel = L"M_{z}", xlabelsize = 20, ylabelsize = 20)
lines!(ax, tl, m_lr, label = L"LR $[M=M(t)]$", linewidth = 2)
lines!(ax, tl, m_me, label = "Fock", linewidth = 2, linestyle = :dash)
axislegend(ax, position = :rb)

ax2 = Axis(fig[1, 2], xlabel = L"\gamma t", ylabel = "Value", xlabelsize = 20, ylabelsize = 20)
lines!(ax2, tl, 1 .- real(sol_lr.fexpect[1, :]), label = L"$1-P$", linewidth = 2)
lines!(
    ax2,
    tl,
    1 .- real(sol_lr.fexpect[3, :]),
    label = L"$1-\mathrm{Tr}(\rho)$",
    linewidth = 2,
    linestyle = :dash,
    color = :orange,
)
lines!(ax2, tl, real(sol_lr.fexpect[2, :]) / Nx / Ny, color = :blue, label = L"S", linewidth = 2)
hlines!(ax2, [Strue], color = :blue, linestyle = :dash, linewidth = 2, label = L"S^{\,\mathrm{true}}_{\mathrm{ss}}")
axislegend(ax2, position = :rb)

fig
```
