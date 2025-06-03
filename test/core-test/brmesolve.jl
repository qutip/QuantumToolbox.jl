@testitem "Bloch-Redfield tensor sec_cutoff" begin
    N = 5
    H = num(N)
    a = destroy(N)
    A_op = a+a'
    spectra(x) = (x>0) * 0.5
    for sec_cutoff in [0, 0.1, 1, 3, -1]
        R = bloch_redfield_tensor(H, [(A_op, spectra)], [a^2], sec_cutoff = sec_cutoff, fock_basis = true)
        R_eig, evecs = bloch_redfield_tensor(H, [(A_op, spectra)], [a^2], sec_cutoff = sec_cutoff, fock_basis = false)
        @test isa(R, QuantumObject)
        @test isa(R_eig, QuantumObject)
        @test isa(evecs, QuantumObject)

        state = rand_dm(N) |> mat2vec
        fock_computed = R * state
        eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state
        @test isapprox(fock_computed, eig_computed, atol = 1e-15)
    end
end

@testitem "Compare brterm and Lindblad" begin
    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'
    spectra(x) = x>0

    # this test applies for limited cutoff
    lindblad = lindblad_dissipator(a)
    computation = brterm(H, A_op, spectra, sec_cutoff = 1.5, fock_basis = true)
    @test isapprox(lindblad, computation, atol = 1e-15)
end

@testitem "brterm basis" begin
    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'
    spectra(x) = x>0
    for sec_cutoff in [0, 0.1, 1, 3, -1]
        R = brterm(H, A_op, spectra, sec_cutoff = sec_cutoff, fock_basis = true)
        R_eig, evecs = brterm(H, A_op, spectra, sec_cutoff = sec_cutoff, fock_basis = false)
        @test isa(R, QuantumObject)
        @test isa(R_eig, QuantumObject)
        @test isa(evecs, QuantumObject)

        state = rand_dm(N) |> mat2vec
        fock_computed = R * state
        eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state
        @test isapprox(fock_computed, eig_computed, atol = 1e-15)
    end;
end

@testitem "brterm sprectra function" begin
    f(x) = exp(x)/10
    function g(x)
        nbar = n_thermal(abs(x), 1)
        if x > 0
            return nbar
        elseif x < 0
            return 1 + nbar
        else
            return 0.0
        end
    end

    spectra_list = [
        x -> (x>0),  # positive frequency filter
        x -> one(x), # no filter
        f, # smooth filter
        g, # thermal field
    ]

    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'
    for spectra in spectra_list
        R = brterm(H, A_op, spectra, sec_cutoff = 0.1, fock_basis = true)
        R_eig, evecs = brterm(H, A_op, spectra, sec_cutoff = 0.1, fock_basis = false)
        @test isa(R, QuantumObject)
        @test isa(R_eig, QuantumObject)
        @test isa(evecs, QuantumObject)

        state = rand_dm(N) |> mat2vec
        fock_computed = R * state
        eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state
        @test isapprox(fock_computed, eig_computed, atol = 1e-15)
    end
end

@testitem "simple qubit system" begin
    pauli_vectors = [sigmax(), sigmay(), sigmaz()]
    γ = 0.25
    spectra(x) = γ * (x>=0)
    _m_c_op = √γ * sigmam()
    _z_c_op = √γ * sigmaz()
    _x_a_op = (sigmax(), spectra)

    arg_sets = [[[_m_c_op], [], [_x_a_op]], [[_m_c_op], [_m_c_op], []], [[_m_c_op, _z_c_op], [_z_c_op], [_x_a_op]]]

    δ = 0
    ϵ = 0.5 * 2π
    e_ops = pauli_vectors
    H = δ * 0.5 * sigmax() + ϵ * 0.5 * sigmaz()
    ψ0 = unit(2basis(2, 0) + basis(2, 1))
    times = LinRange(0, 10, 100)

    for (me_c_ops, brme_c_ops, brme_a_ops) in arg_sets
        me = mesolve(H, ψ0, times, me_c_ops, e_ops = e_ops, progre)
        brme = brmesolve(H, ψ0, times, brme_a_ops, brme_c_ops, e_ops = e_ops)

        @test all(me.expect .== brme.expect)
    end
end;
