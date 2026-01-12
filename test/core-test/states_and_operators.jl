@testitem "States and Operators" begin
    import QuantumToolbox: position, momentum
    using LinearAlgebra
    using SparseArrays

    @testset "zero state" begin
        v1 = zero_ket(4)
        v2 = zero_ket((2, 2))
        @test size(v1) == (4,)
        @test size(v2) == (4,)
        @test v1.data == v2.data
        @test v1 != v2
        @test isket(v1)
        @test isket(v2)
        @test v1.dims == [4]
        @test v2.dims == [2, 2]
    end

    @testset "fock state" begin
        # fock, basis, and fock_dm
        @test fock_dm(4; dims = (2, 2), sparse = Val(true)) ≈ ket2dm(basis(4; dims = (2, 2)))
        @test_throws DimensionMismatch fock(4; dims = 2)
        @test_throws ArgumentError fock(4, 4)
    end

    @testset "coherent state" begin
        # coherent and coherent_dm
        N = 4
        α = 0.25im
        ψα = coherent(N, α)
        ρα = coherent_dm(N, α)
        @test isket(ψα)
        @test isoper(ρα)
        @test ket2dm(ψα) ≈ ρα
        @test tr(ρα) ≈ 1.0
    end

    @testset "thermal state" begin
        N = 5

        # extreme cases
        ρTd0 = thermal_dm(N, 0.0)
        ρTs0 = thermal_dm(N, 0.0; sparse = Val(true))
        ρTd∞ = thermal_dm(N, Inf)
        ρTs∞ = thermal_dm(N, Inf; sparse = Val(true))
        @test tr(ρTd0) ≈ tr(ρTs0) ≈ tr(ρTd∞) ≈ tr(ρTs∞) ≈ 1.0
        @test ρTd0 ≈ ρTs0 ≈ fock_dm(N, 0)
        @test ρTd∞ ≈ ρTs∞ ≈ maximally_mixed_dm(N)

        # general case (also test BigFloat)
        ρTd = thermal_dm(N, big(0.123))
        ρTs = thermal_dm(N, big(0.123); sparse = Val(true))
        @test isoper(ρTd)
        @test ρTd.dims == [N]
        @test tr(ρTd) ≈ tr(ρTs) ≈ 1.0
        @test diag(ρTd) ≈ Float64[
            0.8904859864731106,
            0.09753319353178326,
            0.010682620484781245,
            0.0011700465891612583,
            0.00012815292116369966,
        ]
        @test ρTs.data isa AbstractSparseMatrix
        @test ρTd ≈ ρTs
    end

    @testset "maximally mixed state" begin
        ρ1 = maximally_mixed_dm(4)
        ρ2 = maximally_mixed_dm((2, 2))
        @test size(ρ1) == (4, 4)
        @test size(ρ2) == (4, 4)
        @test tr(ρ1) ≈ 1.0
        @test ρ1.data == ρ2.data
        @test ρ1 != ρ2
        @test isoper(ρ1)
        @test isoper(ρ2)
        @test ρ1.dims == [4]
        @test ρ2.dims == [2, 2]
        @test entropy_vn(ρ1, base = 2) ≈ log(2, 4)
    end

    @testset "random state" begin
        # rand_ket and rand_dm
        ψ = rand_ket(10)
        ρ_AB = rand_dm((2, 2))
        ρ_A = ptrace(ρ_AB, 1)
        ρ_B = ptrace(ρ_AB, 2)
        rank = 5
        ρ_low_rank = rand_dm(10, rank = rank)
        eig_val = eigenenergies(ρ_low_rank)
        @test ψ' * ψ ≈ 1.0
        @test tr(ρ_AB) ≈ 1.0
        @test tr(ρ_A) ≈ 1.0
        @test tr(ρ_B) ≈ 1.0
        @test tr(ρ_low_rank) ≈ 1.0
        @test isposdef(ρ_AB) == true
        @test ishermitian(ρ_AB) == true
        @test ishermitian(ρ_A) == true
        @test ishermitian(ρ_B) == true
        @test ishermitian(ρ_low_rank) == true
        @test all(eigenenergies(ρ_AB) .>= 0)
        @test all(eigenenergies(ρ_A) .>= 0)
        @test all(eigenenergies(ρ_B) .>= 0)
        @test all(isapprox.(eig_val[1:rank], 0.0, atol = 1.0e-10))
        @test all(eig_val[(rank + 1):10] .>= 0)
        @test_throws DomainError rand_dm(4, rank = rank)
        @test_throws DomainError rand_dm(4, rank = 0)
    end

    @testset "entanglement states" begin
        # bell_state, singlet_state, triplet_states, w_state, and ghz_state
        e0 = basis(2, 0)
        e1 = basis(2, 1)
        d0 = basis(3, 0)
        d1 = basis(3, 1)
        d2 = basis(3, 2)
        B00 = (e0 ⊗ e0 + e1 ⊗ e1) / √2
        B01 = (e0 ⊗ e0 - e1 ⊗ e1) / √2
        B10 = (e0 ⊗ e1 + e1 ⊗ e0) / √2
        B11 = (e0 ⊗ e1 - e1 ⊗ e0) / √2
        S = singlet_state()
        T = triplet_states()
        @test bell_state(0, 0) == B00 == ghz_state(2)
        @test bell_state(0, 1) == B01
        @test bell_state(1, 0) == B10 == T[2] == w_state(2)
        @test bell_state(1, 1) == B11 == S
        @test T[1] == e1 ⊗ e1
        @test T[3] == e0 ⊗ e0
        @test w_state(3) == (e0 ⊗ e0 ⊗ e1 + e0 ⊗ e1 ⊗ e0 + e1 ⊗ e0 ⊗ e0) / √3
        @test ghz_state(3) == (e0 ⊗ e0 ⊗ e0 + e1 ⊗ e1 ⊗ e1) / √2
        @test ghz_state(3; d = 3) == (d0 ⊗ d0 ⊗ d0 + d1 ⊗ d1 ⊗ d1 + d2 ⊗ d2 ⊗ d2) / √3
        @test_throws ArgumentError bell_state(0, 2)
        @test_throws ArgumentError bell_state(3, 1)
        @test_throws ArgumentError bell_state(2, 3)
        @test_throws ArgumentError w_state(1)
        @test_throws ArgumentError ghz_state(1)
        @test_throws ArgumentError ghz_state(2; d = 1)
    end

    @testset "bosonic operators" begin
        # destroy, create, num, position, momentum
        n = 10
        i, j = rand(0:(n - 1), 2)
        a = destroy(n)
        ad = create(n)
        N = num(n)
        x = position(n)
        p = momentum(n)
        Pij = projection(n, i, j)
        @test isoper(x)
        @test isoper(p)
        @test a.dims == ad.dims == N.dims == x.dims == p.dims == [n]
        @test eigenenergies(ad * a) ≈ 0:(n - 1)
        @test commutator(N, a) ≈ -a
        @test commutator(N, ad) ≈ ad
        @test all(diag(commutator(x, p))[1:(n - 1)] .≈ 1.0im)
        @test fock(n, i) == Pij * fock(n, j)
        @test_throws ArgumentError projection(n, n, 0)
        @test_throws ArgumentError projection(n, 0, n)
    end

    @testset "displacement and squeezing operators" begin
        N = 10
        α = rand(ComplexF64)
        z = rand(ComplexF64)
        II = qeye(N)
        D = displace(N, α)
        S = squeeze(N, z)
        @test D * D' ≈ II
        @test S * S' ≈ II
    end

    @testset "phase" begin
        # verify Equation (33) in:
        # Michael Martin Nieto, QUANTUM PHASE AND QUANTUM PHASE OPERATORS: Some Physics and Some History
        s = 10
        ϕ0 = 2 * π * rand()
        II = qeye(s + 1)
        ket = [basis(s + 1, i) for i in 0:s]
        SUM = Qobj(zeros(ComplexF64, s + 1, s + 1))
        for j in 0:s
            for k in 0:s
                (j != k) && (
                    SUM +=
                        (exp(1im * (j - k) * ϕ0) / (exp(1im * (j - k) * 2 * π / (s + 1)) - 1)) * ket[j + 1] * ket[k + 1]'
                )
            end
        end
        @test isapprox((ϕ0 + s * π / (s + 1)) * II + (2 * π / (s + 1)) * SUM, phase(s + 1, ϕ0); atol = 1.0e-8)
    end

    @testset "tunneling" begin
        @test tunneling(10, 2) == tunneling(10, 2; sparse = Val(true))
        @test_throws ArgumentError tunneling(10, 0)
    end

    @testset "quantum Fourier transform" begin
        N = 9
        dims = (3, 3)
        ω = exp(2.0im * π / N)
        x = Qobj(rand(ComplexF64, N))
        ψx = basis(N, 0; dims = dims)
        ψk = unit(Qobj(ones(ComplexF64, N); dims = dims))
        F_9 = qft(N)
        F_3_3 = qft(dims)
        y = F_9 * x
        for k in 0:(N - 1)
            nk = collect(0:(N - 1)) * k
            @test y[k + 1] ≈ sum(x.data .* (ω .^ nk)) / sqrt(N)
        end
        @test ψk ≈ F_3_3 * ψx
        @test ψx ≈ F_3_3' * ψk
    end

    @testset "random unitary" begin
        U1 = rand_unitary(20)
        U2 = rand_unitary(20, :exp)
        U3 = rand_unitary((5, 5))
        U4 = rand_unitary((5, 5), :exp)
        @test isunitary(U1)
        @test isunitary(U2)
        @test isunitary(U3)
        @test isunitary(U4)
        @test U1.dims == U2.dims == [20]
        @test U3.dims == U4.dims == [5, 5]

        @test_throws ArgumentError rand_unitary(20, :wrong)
    end

    @testset "Spin-j operators" begin
        # Pauli matrices and general Spin-j operators
        J0 = Qobj(spdiagm(0 => [0.0im]))
        Jx, Jy, Jz = spin_J_set(0.5)
        Sx = sigmax()
        Sy = sigmay()
        Sz = sigmaz()
        @test spin_Jx(0) == spin_Jy(0) == spin_Jz(0) == spin_Jp(0) == spin_Jm(0) == J0
        @test Sx ≈ 2 * spin_Jx(0.5) ≈ 2 * Jx ≈ Qobj(sparse([0.0im 1; 1 0]))
        @test Sy ≈ 2 * spin_Jy(0.5) ≈ 2 * Jy ≈ Qobj(sparse([0.0im -1im; 1im 0]))
        @test Sz ≈ 2 * spin_Jz(0.5) ≈ 2 * Jz ≈ Qobj(sparse([1 0.0im; 0 -1]))
        @test sigmap() ≈ spin_Jp(0.5) ≈ Qobj(sparse([0.0im 1; 0 0]))
        @test sigmam() ≈ spin_Jm(0.5) ≈ Qobj(sparse([0.0im 0; 1 0]))
        @test isunitary(Sx)
        @test isunitary(Sy)
        @test isunitary(Sz)
        for which in [:x, :y, :z, :+, :-]
            @test jmat(2.5, which).dims == [6] # 2.5 * 2 + 1 = 6

            for j_wrong in [-1, 8.7]  # Invalid j
                @test_throws ArgumentError jmat(j_wrong, which)
            end
        end
        @test_throws ArgumentError jmat(0.5, :wrong)

        # test commutation relations for general Spin-j operators
        # the negative signs follow the rule of Levi-Civita symbol (ϵ_ijk)
        j = 2
        Jx, Jy, Jz = spin_J_set(j)
        @test commutator(Jx, Jy) ≈ 1im * Jz
        @test commutator(Jy, Jz) ≈ 1im * Jx
        @test commutator(Jz, Jx) ≈ 1im * Jy
        @test commutator(Jx, Jz) ≈ -1im * Jy
        @test commutator(Jy, Jx) ≈ -1im * Jz
        @test commutator(Jz, Jy) ≈ -1im * Jx
    end

    @testset "spin state" begin
        j = 3.5
        Jz = spin_Jz(j)
        for m in (-j):1:j
            ψm = spin_state(j, m)
            @test Jz * ψm ≈ m * ψm  # check eigenvalue
        end
        @test_throws ArgumentError spin_state(8.7, 8.7)
        @test_throws ArgumentError spin_state(-1, 0)
        @test_throws ArgumentError spin_state(2.5, 2)
        @test_throws ArgumentError spin_state(2.5, 3.5)
        @test_throws ArgumentError spin_state(2.5, -3.5)
    end

    @testset "spin coherent state" begin
        # spin-half (state in Bloch sphere)
        s = 0.5
        θ = π * rand()
        ϕ = 2 * π * rand()
        @test spin_coherent(s, θ, ϕ) ≈ cos(θ / 2) * basis(2, 0) + exp(1im * ϕ) * sin(θ / 2) * basis(2, 1)

        # also verify the equation in:
        # Robert Jones, Spin Coherent States and Statistical Physics, page 3
        s = 1.5
        θ1 = π * rand()
        θ2 = π * rand()
        ϕ1 = 2 * π * rand()
        ϕ2 = 2 * π * rand()
        γ = rand()
        Sz = spin_Jz(s)
        n1 = spin_coherent(s, θ1, ϕ1)
        n2 = spin_coherent(s, θ2, ϕ2)
        @test dot(n1, n2) ≈ (cos(θ1 / 2) * cos(θ2 / 2) + exp(1im * (ϕ2 - ϕ1)) * sin(θ1 / 2) * sin(θ2 / 2))^(2 * s)
        @test dot(n1, exp(-γ * Sz), n2) ≈
            (
            exp(-γ / 2) * cos(θ1 / 2) * cos(θ2 / 2) + exp(1im * (ϕ2 - ϕ1) + γ / 2) * sin(θ1 / 2) * sin(θ2 / 2)
        )^(2 * s)

        # test commutation relations for fermionic creation and annihilation operators
        sites = 4
        SIZE = 2^sites
        dims = ntuple(i -> 2, Val(sites))
        Q_iden = Qobj(sparse((1.0 + 0.0im) * I, SIZE, SIZE); dims = dims)
        Q_zero = Qobj(spzeros(ComplexF64, SIZE, SIZE); dims = dims)
        for i in 1:sites
            d_i = fdestroy(sites, i)
            @test d_i' ≈ fcreate(sites, i)

            for j in 1:sites
                d_j = fdestroy(sites, j)

                if i == j
                    @test commutator(d_i, d_j'; anti = true) ≈ Q_iden
                else
                    @test commutator(d_i, d_j'; anti = true) ≈ Q_zero
                end
                @test commutator(d_i', d_j'; anti = true) ≈ Q_zero
                @test commutator(d_i, d_j; anti = true) ≈ Q_zero
            end
        end
        zero = zero_ket((2, 2))
        vac = tensor(basis(2, 0), basis(2, 0))
        d1 = sigmap() ⊗ qeye(2)
        d2 = sigmaz() ⊗ sigmap()
        @test d1 * vac == zero
        @test d2 * vac == zero
        @test d1' * vac == tensor(basis(2, 1), basis(2, 0))
        @test d2' * vac == tensor(basis(2, 0), basis(2, 1))
        @test d1' * d2' * vac == tensor(basis(2, 1), basis(2, 1))
        @test d1' * d1' * d2' * vac == zero
        @test d2' * d1' * d2' * vac == zero
        @test_throws ArgumentError fdestroy(0, 0)
        @test_throws ArgumentError fdestroy(sites, 0)
        @test_throws ArgumentError fdestroy(sites, sites + 1)
    end

    @testset "identity operator" begin
        I_op1 = qeye(4)
        I_op2 = qeye(4, dims = (2, 2))
        I_su1 = qeye(4, type = SuperOperator())
        I_su2 = qeye(4, type = SuperOperator(), dims = 2)
        @test isunitary(I_op1) == true
        @test isunitary(I_op2) == true
        @test isunitary(I_su1) == false
        @test isunitary(I_su2) == false
        @test I_op1.data == I_op2.data == I_su1.data == I_su2.data
        @test (I_op1 == I_op2) == false
        @test (I_op1 == I_su1) == false
        @test (I_op1 == I_su2) == false
        @test (I_op2 == I_su1) == false
        @test (I_op2 == I_su2) == false
        @test (I_su1 == I_su2) == true
        @test_throws DimensionMismatch qeye(4, dims = 2)
        @test_throws DimensionMismatch qeye(2, type = SuperOperator())
        @test_throws DimensionMismatch qeye(4, type = SuperOperator(), dims = (2, 2))
    end

    @testset "superoperators" begin
        # spre, spost, and sprepost
        Xs = sigmax()
        Xd = to_dense(Xs)
        A_wrong1 = Qobj(rand(4, 4), dims = 4)
        A_wrong2 = Qobj(rand(4, 4), dims = (2, 2))
        A_wrong3 = Qobj(rand(3, 3))
        @test (typeof(spre(Xd).data) <: SparseMatrixCSC) == true
        @test (typeof(spre(Xs).data) <: SparseMatrixCSC) == true
        @test (typeof(spost(Xd).data) <: SparseMatrixCSC) == true
        @test (typeof(spost(Xs).data) <: SparseMatrixCSC) == true
        @test (typeof(sprepost(Xd, Xd).data) <: SparseMatrixCSC) == true
        @test (typeof(sprepost(Xs, Xs).data) <: SparseMatrixCSC) == true
        @test (typeof(sprepost(Xs, Xd).data) <: SparseMatrixCSC) == true
        @test (typeof(sprepost(Xd, Xs).data) <: SparseMatrixCSC) == true
        @test_throws DimensionMismatch sprepost(A_wrong1, A_wrong2)
        @test_throws DimensionMismatch sprepost(A_wrong1, A_wrong3)
    end

    @testset "Type Inference (States)" begin
        @inferred zero_ket(4)
        @inferred zero_ket((2, 2))

        @inferred fock(4, 1, dims = (2, 2), sparse = Val(true))
        @inferred fock(4, 1, dims = (2, 2), sparse = Val(false))
        @inferred basis(4, 1, dims = (2, 2))

        @inferred coherent(5, 0.25im)

        @inferred rand_ket(10)
        @inferred rand_ket((2, 5))

        @inferred fock_dm(4; dims = (2, 2), sparse = Val(true))
        @inferred fock_dm(4; dims = (2, 2), sparse = Val(false))

        @inferred coherent_dm(5, 0.25im)

        @inferred thermal_dm(10, 0.123)
        @inferred thermal_dm(10, 0.123; sparse = Val(true))

        @inferred maximally_mixed_dm(4)
        @inferred maximally_mixed_dm((2, 2))

        @inferred rand_dm((2, 2))
        @inferred rand_dm(10, rank = 5)

        @inferred spin_state(3.5, 1.5)
        @inferred spin_coherent(0.5, π, 2π)

        @inferred bell_state(Val(0), Val(0))
        @inferred bell_state(Val(0), Val(1))
        @inferred bell_state(Val(1), Val(0))
        @inferred bell_state(Val(1), Val(1))
        @inferred singlet_state()
        @inferred triplet_states()
        @inferred w_state(Val(2))

        @inferred ghz_state(Val(2))
        @inferred ghz_state(Val(3))
        @inferred ghz_state(Val(3); d = 3)
    end

    @testset "Type Inference (Operators)" begin
        @inferred rand_unitary(10, Val(:haar))
        @inferred rand_unitary(10, Val(:exp))

        a = destroy(20)
        a_d = create(20)
        @inferred commutator(a, a_d)

        @inferred destroy(20)
        @inferred create(20)
        @inferred num(20)
        @inferred displace(20, 0.5 + 0.5im)
        @inferred squeeze(20, 0.5 + 0.5im)
        @inferred position(20)
        @inferred momentum(20)

        @inferred phase(20, 0.5)

        @inferred jmat(2.5)
        @inferred jmat(2.5, Val(:x))
        @inferred jmat(2.5, Val(:y))
        @inferred jmat(2.5, Val(:z))

        @inferred sigmap()
        @inferred sigmam()
        @inferred sigmax()
        @inferred sigmay()
        @inferred sigmaz()

        @inferred eye(20)

        @inferred fdestroy(Val(10), 4)
        @inferred fcreate(Val(10), 4)

        @inferred projection(20, 5, 3)

        @inferred tunneling(20, 1, sparse = Val(false))
        @inferred tunneling(20, 1, sparse = Val(true))

        @inferred qft(20)
        @inferred qft((2, 10))
    end

    @testset "Element types" begin
        N = 2
        type_list = [ComplexF32, Complex{BigFloat}]
        for T in type_list
            @test T == eltype(zero_ket(T, N))
            @test T == eltype(fock(T, N, 0; sparse = Val(true)))
            @test T == eltype(fock(T, N, 0; sparse = Val(false)))
            @test T == eltype(rand_ket(T, N))
        end
    end
end
