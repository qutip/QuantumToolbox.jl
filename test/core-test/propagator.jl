@testitem "Propagator (by solvers)" begin
    ϵ0 = 1.0 * 2π
    Δ = 0.8 * 2π
    H = (ϵ0/2) * sigmaz() + (Δ/2) * sigmax()
    L = liouvillian(H)
    ψ0 = basis(2, 0)
    ρ0 = mat2vec(ket2dm(ψ0))

    dt = π/5
    tlist = 0:dt:(2π)
    ψt = sesolve(H, ψ0, tlist; progress_bar = Val(false)).states[2:end] # ignore the initial state
    ρt = mesolve(H, ρ0, tlist; progress_bar = Val(false)).states[2:end] # ignore the initial state
    Prop_se = sesolve(H, qeye_like(H), [0, dt]; progress_bar = Val(false)).states[end]
    Prop_me = mesolve(L, qeye_like(L), [0, dt]; progress_bar = Val(false)).states[end]

    for n in 1:(length(tlist)-1)
        @test isapprox(Prop_se^n * ψ0, ψt[n]; atol = 1e-5)
        @test isapprox(Prop_me^n * ρ0, ρt[n]; atol = 1e-5)
    end
end
