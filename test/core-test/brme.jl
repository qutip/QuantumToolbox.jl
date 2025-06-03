@testset "Bloch-Redfield tensor sec_cutoff = $sec_cutoff" for sec_cutoff in [0,0.1,1,3,-1]
    N = 5
    H = num(N)
    a = destroy(N)
    A_op = a+a'
    spectra(x) = (x>0) * 0.5

    R = bloch_redfield_tensor(
        H,
        [(A_op, spectra)],
        [a^2],
        sec_cutoff=sec_cutoff,
        fock_basis=true
    )
    R_eig, evecs = bloch_redfield_tensor(
        H,
        [(A_op, spectra)],
        [a^2],
        sec_cutoff=sec_cutoff,
        fock_basis=false
    )
    @test isa(R, QuantumObject)
    @test isa(R_eig, QuantumObject)
    @test isa(evecs, QuantumObject)
    state = rand_dm(N) |> mat2vec
    fock_computed = R * state
    eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state

    @test isapprox(fock_computed, eig_computed, atol=1e-15)
end;

@testset "brterm lindblad" begin
    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'
    spectra(x) = x>0
    
    # this test applys for limited cutoff
    lindblad = lindblad_dissipator(a)
    computation = brterm(H, A_op, spectra, sec_cutoff=1.5, fock_basis=true)
    @test isapprox(lindblad, computation, atol=1e-15)
end;

@testset "brterm basis for sec_cutoff = $sec_cutoff" for sec_cutoff in [0,0.1,1,3,-1]
    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'
    spectra(x) = x>0    
    
    R = brterm(
        H,
        A_op,
        spectra,
        sec_cutoff=sec_cutoff,
        fock_basis=true
    )
    R_eig, evecs = brterm(
        H,
        A_op,
        spectra,
        sec_cutoff=sec_cutoff,
        fock_basis=false
    )
    @test isa(R, QuantumObject)
    @test isa(R_eig, QuantumObject)
    @test isa(evecs, QuantumObject)
    state = rand_dm(N) |> mat2vec
    fock_computed = R * state
    eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state

    @test isapprox(fock_computed, eig_computed, atol=1e-15)
end;

f(x) = exp(x)/10
function g(x) 
    nbar = n_thermal(abs(x), 1)
    if x > 0
        return nbar
    elseif x < 0
        return 1 + nbar
    else
        return 0.
    end
end

spectras = [
    (x -> (x>0), "positive frequency filter"),
    (x -> one(x), "no filter"),
    (f, "smooth filter"),
    (g, "thermal field")
]

@testset "brterms sprectra with $description" for (spectra, description) in spectras
    N = 5
    H = num(N)
    a = destroy(N) + destroy(N)^2/2
    A_op = a+a'

    R = brterm(
        H,
        A_op,
        spectra,
        sec_cutoff=0.1,
        fock_basis=true
    )
    R_eig, evecs = brterm(
        H,
        A_op,
        spectra,
        sec_cutoff=0.1,
        fock_basis=false
    )
    @test isa(R, QuantumObject)
    @test isa(R_eig, QuantumObject)
    @test isa(evecs, QuantumObject)
    state = rand_dm(N) |> mat2vec
    fock_computed = R * state
    eig_computed = (sprepost(evecs, evecs') * R_eig * sprepost(evecs', evecs)) * state

    @test isapprox(fock_computed, eig_computed, atol=1e-15)
end;    