using BenchmarkTools
using QuPhys

const SUITE = BenchmarkGroup()

# Jaynes-Cummings Model
a = kron(destroy(10), eye(2))
sm = kron(eye(10), sigmam())
sp = kron(eye(10), sigmap())
sz = kron(eye(10), sigmaz())

ω = 1
g = 0.05
F = 0.005
γ = 0.01
T = 0.05
nth = n_th(ω, T)

H = ω*a'*a + ω*sz/2 + g*(a*sp + a'*sm)
c_ops = [sqrt(γ*(nth+1))*a, sqrt(γ*nth)*a', sqrt(γ)*sm]
e_ops = [a'*a, sz]
ψ0 = kron(fock(10, 4), fock(2, 1))
tlist = range(0, 10/γ, 1000)
fields = [a+a', sm+sp]
γ_list = [γ, γ]
ω_list = [ω, ω]
T_list = [T, 0]

liouvillian(H, c_ops) # precompile
steadystate(H, c_ops) # precompile
steadystate_floquet(H, c_ops, -1im * 0.5 * F * (a+a'), 1im * 0.5 * F * (a+a'), 1) # precompile
liouvillian_generalized(H, fields, γ_list, ω_list, T_list, N_trunc=10) # precompile
liouvillian_generalized(H, fields, γ_list, ω_list, T_list) # precompile
mesolve(H, ψ0, tlist, c_ops, e_ops=e_ops, progress=false, saveat=[tlist[end]]) # precompile
mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, progress=false, ensemble_method=EnsembleSerial(), n_traj=500, saveat=[tlist[end]]) # precompile
mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, progress=false, n_traj=500, saveat=[tlist[end]]) # precompile

SUITE["timeevolution"] = BenchmarkGroup()
SUITE["timeevolution"]["liouvillian"] = @benchmarkable liouvillian($H, $c_ops)
SUITE["timeevolution"]["liouvillian_generalized"] = @benchmarkable liouvillian_generalized($H, $fields, $γ_list, $ω_list, $T_list)
SUITE["timeevolution"]["liouvillian_generalized_trunc"] = @benchmarkable liouvillian_generalized($H, $fields, $γ_list, $ω_list, $T_list, N_trunc=10)
SUITE["timeevolution"]["steadystate"] = @benchmarkable steadystate($H, $c_ops)
SUITE["timeevolution"]["steadystate_floquet"] = @benchmarkable steadystate_floquet($H, $c_ops, -1im * 0.5 * F * $(a+a'), 1im * 0.5 * F * $(a+a'), $ω)
SUITE["timeevolution"]["mesolve"] = @benchmarkable mesolve($H, $ψ0, $tlist, $c_ops, e_ops=$e_ops, progress=false, saveat=[$tlist[end]])
SUITE["timeevolution"]["mcsolve_serial"] = @benchmarkable mcsolve($H, $ψ0, $tlist, $c_ops, e_ops=$e_ops, progress=false, n_traj=500, ensemble_method=EnsembleSerial(), saveat=[$tlist[end]])
SUITE["timeevolution"]["mcsolve_parallel"] = @benchmarkable mcsolve($H, $ψ0, $tlist, $c_ops, e_ops=$e_ops, progress=false, n_traj=500, saveat=[$tlist[end]])

# Entanglement
ψ = normalize( kron(fock(20, 0), fock(10, 2)) + kron(fock(20, 1), fock(10, 1)) )
ρ = ket2dm(ψ)
ρ1 = ptrace(ρ, [1])

ptrace(ψ, [1]) # precompile
ptrace(ρ, [1]) # precompile
entropy_vn(ρ1) # precompile

SUITE["entanglement"] = BenchmarkGroup()
SUITE["entanglement"]["ptrace"] = @benchmarkable ptrace($ψ, [1])
SUITE["entanglement"]["ptrace_dm"] = @benchmarkable ptrace($ρ, [1])
SUITE["entanglement"]["entropy_vn"] = @benchmarkable entropy_vn($ρ1)

# Wigner
α = 0.5 + 0.8im
ψ = coherent(30, α)
ρ = dense_to_sparse(ket2dm(ψ), 1e-6)
ψ2 = dense_to_sparse(normalize(fock(20, 1) + fock(20, 3)))
xvec = LinRange(-3, 3, 300)
yvec = LinRange(-3, 3, 300)

wigner(ψ, xvec, yvec, solver=WignerLaguerre(tol=1e-6)) # precompile
wigner(ρ, xvec, yvec, solver=WignerLaguerre(parallel=false)) # precompile
wigner(ρ, xvec, yvec, solver=WignerLaguerre(parallel=true)) # precompile
wigner(ψ, xvec, yvec, solver=WignerClenshaw()) # precompile

SUITE["wigner"] = BenchmarkGroup()
SUITE["wigner"]["wigner_clenshaw"] = @benchmarkable wigner($ψ, $xvec, $yvec, solver=WignerClenshaw())
SUITE["wigner"]["wigner_laguerre"] = @benchmarkable wigner($ψ, $xvec, $yvec, solver=WignerLaguerre(tol=1e-6))
SUITE["wigner"]["wigner_laguerre_sparse"] = @benchmarkable wigner($ρ, $xvec, $yvec, solver=WignerLaguerre(parallel=false))
SUITE["wigner"]["wigner_laguerre_sparse_parallel"] = @benchmarkable wigner($ρ, $xvec, $yvec, solver=WignerLaguerre(parallel=true))
SUITE["wigner"]["wigner_laguerre_fock"] = @benchmarkable wigner($ψ2, $xvec, $yvec, solver=WignerLaguerre())
SUITE["wigner"]["wigner_clenshaw_fock"] = @benchmarkable wigner($ψ2, $xvec, $yvec, solver=WignerClenshaw())

# Permutation
N = 50
Δ = 0
G = 5
tg = 0
θ  = atan(tg)
U  = sin(θ)
κ2 = cos(θ)
κ1  = 0.
κϕ  = 1e-3
nth = 0.

a     = destroy(N)
ad    = create(N)
H     = -Δ*ad*a + G/2*(ad^2 + a^2) + U/2*(ad^2*a^2)
c_ops = [√(κ2)*a^2, √(κ1*(nth+1))*a, √(κ1*nth)*ad, √(κϕ)*ad*a]
L     = liouvillian(H,c_ops)

P, L_bd, block_sizes = bdf(L) # precompile
blocks_list, block_indices = get_bdf_blocks(L_bd, block_sizes) # precompile

SUITE["permutation"] = BenchmarkGroup()
SUITE["permutation"]["bdf"] = @benchmarkable bdf($L)
SUITE["permutation"]["get_bdf_blocks"] = @benchmarkable get_bdf_blocks($L_bd, $block_sizes)

# Correlation and spectrum
a = destroy(10)
H = a' * a
c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']

ω_list = range(0, 3, length=1000)
spectrum(H, ω_list, a', a, c_ops, solver=FFTCorrelation(), 
        progress=false, abstol=1e-7, reltol=1e-5) # precompile
spectrum(H, ω_list, a', a, c_ops) # precompile

SUITE["spectrum"] = BenchmarkGroup()
SUITE["spectrum"]["spectrum_fft"] = @benchmarkable spectrum($H, $ω_list, $(a'), $a, $c_ops, solver=FFTCorrelation(), progress=false, abstol=1e-7, reltol=1e-5)
SUITE["spectrum"]["spectrum_exponential_series"] = @benchmarkable spectrum($H, $ω_list, $(a'), $a, $c_ops)
