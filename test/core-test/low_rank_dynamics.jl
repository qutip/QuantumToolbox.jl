@testset "Low Rank Dynamics" begin
    # Define lattice
    Nx, Ny = 2, 3
    latt = Lattice(Nx = Nx, Ny = Ny)
    ##
    N_cut = 2         # Number of states of each mode
    N_modes = latt.N  # Number of modes
    N = N_cut^N_modes # Total number of states
    M = latt.N + 1       # Number of states in the LR basis

    # Define initial state
    ϕ = Vector{QuantumObject{Vector{ComplexF64},KetQuantumObject,M - 1}}(undef, M)
    ϕ[1] = kron(fill(basis(2, 1), N_modes)...)

    i = 1
    for j in 1:N_modes
        i += 1
        i <= M && (ϕ[i] = SingleSiteOperator(sigmap(), j, latt) * ϕ[1])
    end
    for k in 1:N_modes-1
        for l in k+1:N_modes
            i += 1
            i <= M && (ϕ[i] = SingleSiteOperator(sigmap(), k, latt) * SingleSiteOperator(sigmap(), l, latt) * ϕ[1])
        end
    end
    for i in i+1:M
        ϕ[i] = QuantumObject(rand(ComplexF64, size(ϕ[1])[1]), dims = ϕ[1].dims)
        normalize!(ϕ[i])
    end

    z = hcat(get_data.(ϕ)...)
    B = Matrix(Diagonal([1 + 0im; zeros(M - 1)]))
    S = z' * z # Overlap matrix
    B = B / tr(S * B) # Normalize B
    ρ = Qobj(z * B * z', dims = ntuple(i -> N_cut, Val(N_modes))) # Full density matrix

    # Define Hamiltonian and collapse operators
    Jx = 0.9
    Jy = 1.04
    Jz = 1.0
    hx = 0.0
    hy = 0.0
    hz = 0.0
    γ = 1

    Sx = mapreduce(i -> SingleSiteOperator(sigmax(), i, latt), +, 1:latt.N)
    Sy = mapreduce(i -> SingleSiteOperator(sigmay(), i, latt), +, 1:latt.N)
    Sz = mapreduce(i -> SingleSiteOperator(sigmaz(), i, latt), +, 1:latt.N)
    SFxx = mapreduce(
        x -> SingleSiteOperator(sigmax(), x[1], latt) * SingleSiteOperator(sigmax(), x[2], latt),
        +,
        Iterators.product(1:latt.N, 1:latt.N),
    )

    H, c_ops = DissipativeIsing(Jx, Jy, Jz, hx, hy, hz, γ, latt; boundary_condition = Val(:periodic_bc), order = 1)
    e_ops = (Sz,)

    tl = range(0, 10, 100)

    # Full solution
    sol_me = mesolve(H, ρ, tl, c_ops; e_ops = [e_ops...])
    Strue = entropy_vn(sol_me.states[end], base = 2) / latt.N

    # Low rank solution
    function f_entropy(p, z, B)
        C = p.A0
        σ = p.Bi

        mul!(C, z, sqrt(B))
        mul!(σ, C', C)
        return entropy_vn(Qobj(Hermitian(σ), type = Operator), base = 2)
    end

    opt = (err_max = 1e-3, p0 = 0.0, atol_inv = 1e-6, adj_condition = "variational", Δt = 0.0, progress = false)

    sol_lr = lr_mesolve(H, z, B, tl, c_ops; e_ops = e_ops, f_ops = (f_entropy,), opt = opt)

    # Test
    S_lr = real(sol_lr.fexpect[1, end]) / latt.N

    @test real(sol_me.expect[1, :]) ≈ real(sol_lr.expect[1, :]) atol = 1e-1
    @test S_lr ≈ Strue atol = 1e-1
end
