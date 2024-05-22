@testset "Low Rank Dynamics" begin
    # Define lattice
    Nx, Ny = 2, 3
    latt = Lattice(Nx = Nx, Ny = Ny)
    N_cut = 2
    N_modes = latt.N
    N = N_cut^N_modes
    M = Nx * Ny + 1

    # Define initial state
    ϕ = Vector{QuantumObject{Vector{ComplexF64},KetQuantumObject}}(undef, M)
    ϕ[1] = tensor(repeat([basis(2, 0)], N_modes)...)
    i = 1
    for j in 1:N_modes
        i += 1
        i <= M && (ϕ[i] = mb(sp, j, latt) * ϕ[1])
    end
    for k in 1:N_modes-1
        for l in k+1:N_modes
            i += 1
            i <= M && (ϕ[i] = mb(sp, k, latt) * mb(sp, l, latt) * ϕ[1])
        end
    end
    for i in i+1:M
        ϕ[i] = Qobj(rand(ComplexF64, size(ϕ[1])[1]), dims = ϕ[1].dims)
        normalize!(ϕ[i])
    end
    z = hcat(broadcast(x -> x.data, ϕ)...)
    B = Matrix(Diagonal([1 + 0im; zeros(M - 1)]))
    S = z' * z
    B = B / tr(S * B)
    ρ = Qobj(z * B * z', dims = ones(Int, N_modes) * N_cut)

    # Define Hamiltonian and collapse operators
    Jx = 0.9
    Jy = 1.02
    Jz = 1.0
    hx = 0.0
    γ = 1
    Sz = sum([mb(sz, i, latt) for i in 1:latt.N])
    tl = LinRange(0, 10, 100)

    H, c_ops = TFIM(Jx, Jy, Jz, hx, γ, latt; bc = pbc, order = 1)
    e_ops = (Sz,)

    # Full solution
    mesol = mesolve(H, ρ, tl, c_ops; e_ops = [e_ops...], progress_bar = false)
    A = Matrix(mesol.states[end].data)
    λ = eigvals(Hermitian(A))
    Strue = -sum(λ .* log2.(λ))

    # Low rank solution
    function f_entropy(p, z, B)
        C = p.A0
        σ = p.Bi
        mul!(C, z, sqrt(B))
        mul!(σ, C', C)
        λ = eigvals(Hermitian(σ))
        λ = λ[λ.>1e-10]
        return -sum(λ .* log2.(λ))
    end

    opt = LRMesolveOptions(
        alg = Tsit5(),
        err_max = 1e-3,
        p0 = 0.0,
        atol_inv = 1e-6,
        adj_condition = "variational",
        Δt = 0.2,
        progress = false,
    )
    lrsol = lr_mesolve(H, z, B, tl, c_ops; e_ops = e_ops, f_ops = (f_entropy,), opt = opt)

    # Test
    m_me = real(mesol.expect[1, :])
    m_lr = real(lrsol.expvals[1, :])
    @test all(abs.((m_me .- m_lr) ./ m_me) .< 0.1)

    S_lr = real(lrsol.funvals[1, end])
    @test abs((S_lr - Strue) / Strue) < 0.5
end
