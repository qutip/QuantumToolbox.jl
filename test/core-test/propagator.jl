@testitem "Propagators" begin
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()
    I = one(sz)

    H = sz + QobjEvo((0.1*sx, (p, t) -> sin(2*t)))

    p_se = propagator(H)

    c_ops = [0.1*sz]
    p_me = propagator(H, c_ops)

    sesolve_res = []
    mesolve_res = []

    t0 = 0.0
    t = 10.0
    for i in 0:1
        psi = fock(2, i)
        push!(sesolve_res, sesolve(H, psi, [t0, t]))
    end

    sup_psi_base = mat2vec(fock(2, 0)*fock(2, 0)')
    for i in 1:4
        sup_psi = 0 * sup_psi_base
        sup_psi[i] = 1
        push!(mesolve_res, mesolve(H, vec2mat(sup_psi), [t0, t], c_ops))
    end

    U_se = p_se(t)
    U_me = p_me(t)

    U_sesolve = QuantumObject(hcat([sesolve_res[i].states[end].data for i in 1:2]...))
    U_mesolve = QuantumObject(hcat([mat2vec(mesolve_res[i].states[end]).data for i in 1:4]...), type = SuperOperator())

    @test norm(U_se-U_sesolve) < 1e-8
    @test norm(U_me-U_mesolve) < 1e-8

    U_se = p_se(t, 3.0)
    U_se_JI = p_se(t, 3.0; just_interval = true)
    U_me = p_me(t, 3.0)
    U_me_JI = p_me(t, 3.0; just_interval = true)

    @test norm(U_se-U_se_JI)<1e-6
    @test norm(U_me-U_me_JI) < 1e-6

end
