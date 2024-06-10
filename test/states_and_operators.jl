@testset "States and Operators" begin
    # fock, basis
    @test fock(4; dims = [2, 2], sparse = true) ≈ basis(4; dims = [2, 2])
    @test_throws DimensionMismatch fock(4; dims = [2])

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

    # spre, spost, and sprepost
    Xs = sigmax()
    Xd = sparse_to_dense(Xs)
    A_wrong1 = Qobj(rand(4, 4), dims = [4])
    A_wrong2 = Qobj(rand(4, 4), dims = [2, 2])
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
