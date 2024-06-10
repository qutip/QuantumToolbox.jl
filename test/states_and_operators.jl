@testset "States and Operators" begin
    # zero_ket
    v1 = zero_ket(4)
    v2 = zero_ket([2, 2])
    @test size(v1) == (4,)
    @test size(v2) == (4,)
    @test v1.data == v2.data
    @test v1 != v2
    @test isket(v1)
    @test isket(v2)
    @test v1.dims == [4]
    @test v2.dims == [2, 2]

    # fock, basis, fock_dm
    @test fock_dm(4; dims = [2, 2], sparse = true) ≈ ket2dm(basis(4; dims = [2, 2]))
    @test_throws DimensionMismatch fock(4; dims = [2])

    # coherent, coherent_dm
    N = 4
    α = 0.25im
    ψα = coherent(N, α)
    ρα = coherent_dm(N, α)
    @test isket(ψα)
    @test isoper(ρα)
    @test ket2dm(ψα) ≈ ρα

    # thermal_dm
    ρT = thermal_dm(5, 0.123)
    @test isoper(ρT)
    @test ρT.dims == [5]
    @test ρT.data ≈ spdiagm(
        0 => Float64[
            0.8904859864731106,
            0.09753319353178326,
            0.010682620484781245,
            0.0011700465891612583,
            0.00012815292116369966,
        ],
    )

    # maximally_mixed_dm
    ρ1 = maximally_mixed_dm(4)
    ρ2 = maximally_mixed_dm([2, 2])
    @test size(ρ1) == (4, 4)
    @test size(ρ2) == (4, 4)
    @test ρ1.data == ρ2.data
    @test ρ1 != ρ2
    @test isoper(ρ1)
    @test isoper(ρ2)
    @test ρ1.dims == [4]
    @test ρ2.dims == [2, 2]

    # rand_dm
    ρ_AB = rand_dm(4, dims = [2, 2])
    ρ_A = ptrace(ρ_AB, 1)
    ρ_B = ptrace(ρ_AB, 2)
    @test tr(ρ_AB) ≈ 1.0
    @test tr(ρ_A) ≈ 1.0
    @test tr(ρ_B) ≈ 1.0
    @test ishermitian(ρ_AB) == true
    @test ishermitian(ρ_A) == true
    @test ishermitian(ρ_B) == true
    @test all(eigenenergies(ρ_AB) .>= 0)
    @test all(eigenenergies(ρ_A) .>= 0)
    @test all(eigenenergies(ρ_B) .>= 0)
    @test_throws DimensionMismatch rand_dm(4, dims = [2])

    # Pauli matrices and general Spin-j operators
    J0 = Qobj(spdiagm(0 => [0.0im]))
    Jx, Jy, Jz = spin_J_set(0.5)
    @test spin_Jx(0) == spin_Jy(0) == spin_Jz(0) == spin_Jp(0) == spin_Jm(0) == J0
    @test sigmax() ≈ 2 * spin_Jx(0.5) ≈ 2 * Jx ≈ Qobj(sparse([0.0im 1; 1 0]))
    @test sigmay() ≈ 2 * spin_Jy(0.5) ≈ 2 * Jy ≈ Qobj(sparse([0.0im -1im; 1im 0]))
    @test sigmaz() ≈ 2 * spin_Jz(0.5) ≈ 2 * Jz ≈ Qobj(sparse([1 0.0im; 0 -1]))
    @test sigmap() ≈ spin_Jp(0.5) ≈ Qobj(sparse([0.0im 1; 0 0]))
    @test sigmam() ≈ spin_Jm(0.5) ≈ Qobj(sparse([0.0im 0; 1 0]))
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

    # test commutation relations for fermionic creation and annihilation operators
    sites = 4
    SIZE = 2^sites
    dims = fill(2, sites)
    Q_iden = Qobj(sparse((1.0 + 0.0im) * LinearAlgebra.I, SIZE, SIZE); dims = dims)
    Q_zero = Qobj(spzeros(ComplexF64, SIZE, SIZE); dims = dims)
    for i in 0:(sites-1)
        d_i = fdestroy(sites, i)
        @test d_i' ≈ fcreate(sites, i)

        for j in 0:(sites-1)
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
    @test_throws ArgumentError fdestroy(0, 0)
    @test_throws ArgumentError fdestroy(sites, -1)
    @test_throws ArgumentError fdestroy(sites, sites)

    # identity operator
    I_op1 = qeye(4)
    I_op2 = qeye(4, dims = [2, 2])
    I_su1 = qeye(4, type = SuperOperator)
    I_su2 = qeye(4, type = SuperOperator, dims = [2])
    @test I_op1.data == I_op2.data == I_su1.data == I_su2.data
    @test (I_op1 == I_op2) == false
    @test (I_op1 == I_su1) == false
    @test (I_op1 == I_su2) == false
    @test (I_op2 == I_su1) == false
    @test (I_op2 == I_su2) == false
    @test (I_su1 == I_su2) == true
    @test_throws DimensionMismatch qeye(4, dims = [2])
    @test_throws DimensionMismatch qeye(2, type = SuperOperator)
    @test_throws DimensionMismatch qeye(4, type = SuperOperator, dims = [2, 2])
end
