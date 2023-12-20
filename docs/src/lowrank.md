# Low Rank Master Equation

We start by importing the packages

```@example lowrank
using Plots
using LaTeXStrings
using QuPhys;
```

Define lattice
```@example lowrank
Nx, Ny = 2, 3
latt = Lattice(Nx=Nx, Ny=Ny)
```

Define lr-space dimensions
```@example lowrank
N_cut = 2         # Number of states of each mode
N_modes = latt.N  # Number of modes
N = N_cut^N_modes # Total number of states
M = Nx*Ny+1       # Number of states in the LR basis
```

Define lr states. Take as initial state all spins up. All other N states are taken as those with miniman Hamming distance to the initial state.
```@example lowrank
ϕ = Vector{QuantumObject{Vector{ComplexF64}, KetQuantumObject}}(undef, M)
ϕ[1] = kron(repeat([basis(2,0)],N_modes)...)

global i=1
for j in 1:N_modes
    global i+=1
    i<=M && (ϕ[i] = mb(sp, j, latt) * ϕ[1])
end
for k in 1:N_modes-1
    for l=k+1:N_modes
        global i+=1
        i<=M && (ϕ[i] = mb(sp, k, latt) * mb(sp, l, latt) * ϕ[1])
    end
end
for i in i+1:M
    ϕ[i] = QuantumObject(rand(ComplexF64,size(ϕ[1])[1]), dims=ϕ[1].dims)
    normalize!(ϕ[i])
end
```

Define the initial state
```@example lowrank
z  = hcat(broadcast(x->x.data, ϕ)...)
p0 = 0. # Population of the lr states other than the initial state
B  = Matrix(Diagonal([1+0im; p0 * ones(M-1)]))
S  = z'*z # Overlap matrix
B = B / tr(S*B) # Normalize B

ρ = QuantumObject(z*B*z', dims=ones(Int,N_modes)*N_cut); # Full density matrix
```

Define the Hamiltonian and collapse operators
```@example lowrank
# Define Hamiltonian and collapse operators
Jx = 0.9 
Jy = 1.04
Jz = 1. 
hx = 0.
γ  = 1

Sx = sum([mb(sx, i, latt) for i in 1:latt.N])
Sy = sum([mb(sy, i, latt) for i in 1:latt.N])
Sz = sum([mb(sz, i, latt) for i in 1:latt.N])
SFxx = sum([mb(sx, i, latt) * mb(sx, j, latt) for i in 1:latt.N for j in 1:latt.N])

H, c_ops = TFIM(Jx, Jy, Jz, hx, γ, latt; bc=pbc, order=1)
e_ops = (Sx,Sy,Sz,SFxx)

tl = LinRange(0,10,100);
```

### Full evolution
```@example lowrank
@time mesol = mesolve(H, ρ, tl, c_ops; e_ops=[e_ops...]);
A = Matrix(mesol.states[end].data)
λ = eigvals(Hermitian(A))
Strue  = -sum(λ.*log2.(λ))/latt.N;
```

### Low Rank evolution
Define functions to be evaluated during the low-rank evolution
```@example lowrank
function f_purity(p,z,B)
    N  = p.N
    M  = p.M
    S  = p.S
    T  = p.temp_MM

    mul!(T, S, B)
    tr(T^2)
end

function f_trace(p,z,B)
    N  = p.N
    M  = p.M
    S  = p.S
    T  = p.temp_MM

    mul!(T,S,B)
    tr(T)
end

function f_entropy(p,z,B)
    C = p.A0
    σ = p.Bi

    mul!(C, z, sqrt(B))
    mul!(σ, C', C)
    λ = eigvals(Hermitian(σ))
    λ = λ[λ.>1e-10]
    return -sum(λ .* log2.(λ))
end;
```

Define the options for the low-rank evolution
```@example lowrank
opt = LRMesolveOptions( 
    alg        = Tsit5(),
    err_max    = 1e-3,
    p0         = 0.,
    atol_inv   = 1e-6,
    adj_condition="variational",
    Δt = 0.
);

@time lrsol = lr_mesolve(H, z, B, tl, c_ops; e_ops=e_ops, f_ops=(f_purity, f_entropy, f_trace,), opt=opt);
```

Plot the results
```@example lowrank
fig = plot(layout=(1,2), size=(800,400), legend=:topleft, xlabel=L"\gamma t")

m_me = real(mesol.expect[3,:])/Nx/Ny
m_lr = real(lrsol.expvals[3,:])/Nx/Ny

plot!(fig[1], tl, m_lr, label=raw"LR $[M=M(t)]$", lw=2)
plot!(fig[1], tl, m_me, ls=:dash, label="Fock", lw=2)
ylabel!(fig[1], L"M_{z}")

plot!(fig[2], tl, 1 .-real(lrsol.funvals[1,:]), label=L"$1-P$", lw=2)
plot!(fig[2], tl, 1 .-real(lrsol.funvals[3,:]), c=:orange, label=L"$1-\rm{Tr}(\rho)$", lw=2, ls=:dash)
plot!(fig[2], tl, real(lrsol.funvals[2,:])/Nx/Ny, c=:blue, label=L"S", lw=2)
hline!(fig[2], [Strue], c=:blue, ls=:dash, lw=2, label=L"S^{\rm \,true}_{\rm ss}")
ylabel!(fig[2], "value")
xlabel!(fig[2], L"\gamma t")
```