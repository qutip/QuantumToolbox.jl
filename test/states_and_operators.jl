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
end
