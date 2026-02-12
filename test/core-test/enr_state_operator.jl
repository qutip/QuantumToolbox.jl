@testitem "Excitation number restricted state space" begin
    using StaticArraysCore
    using SparseArrays

    @testset "EnrSpace" begin
        s_enr = EnrSpace((2, 2, 3), 3)

        # check if the idx2state is the same as qutip
        qutip_idx2state = Dict(
            1 => SVector{3}(0, 0, 0),
            2 => SVector{3}(0, 0, 1),
            3 => SVector{3}(0, 0, 2),
            4 => SVector{3}(0, 1, 0),
            5 => SVector{3}(0, 1, 1),
            6 => SVector{3}(0, 1, 2),
            7 => SVector{3}(1, 0, 0),
            8 => SVector{3}(1, 0, 1),
            9 => SVector{3}(1, 0, 2),
            10 => SVector{3}(1, 1, 0),
            11 => SVector{3}(1, 1, 1),
        )
        @test s_enr.idx2state == qutip_idx2state
    end

    @testset "Element type" begin
        s_enr = EnrSpace((2, 2, 3), 3)
        float_type_list = [Float32, BigFloat]
        for FT in float_type_list
            CT = Complex{FT}
            @test CT == eltype(enr_fock(CT, s_enr, zeros(Int, 3)))
            @test FT == eltype(enr_thermal_dm(s_enr, rand(FT); sparse = Val(true)))
            @test FT == eltype(enr_thermal_dm(s_enr, rand(FT, 3); sparse = Val(false)))
            @test all(==(CT), eltype.(enr_destroy(CT, s_enr)))
            @test CT == eltype(enr_identity(CT, s_enr))
        end
    end

    @testset "kron" begin
        # normal Hilbert space
        D1 = 4
        D2 = 5
        dims_s = (D1, D2)
        ρ_s = rand_dm(dims_s)
        I_s = qeye(D1) ⊗ qeye(D2)
        size_s = prod(dims_s)
        space_s = TensorSpace(Space(D1), Space(D2))

        # EnrSpace
        dims_enr = (2, 3, 2)
        excitations = 4
        space_enr = EnrSpace(dims_enr, excitations)
        I_enr = enr_identity(space_enr)
        size_enr = space_enr.size

        # enr_thermal_dm (extreme cases)
        ρTd0 = enr_thermal_dm(space_enr, 0.0)
        ρTs0 = enr_thermal_dm(space_enr, 0.0; sparse = Val(true))
        ρTd∞ = enr_thermal_dm(space_enr, Inf)
        ρTs∞ = enr_thermal_dm(space_enr, Inf; sparse = Val(true))
        @test tr(ρTd0) ≈ tr(ρTs0) ≈ tr(ρTd∞) ≈ tr(ρTs∞) ≈ 1.0
        @test ρTd0.data ≈ ρTs0.data ≈ fock_dm(size_enr, 0).data
        @test ρTd∞.data ≈ ρTs∞.data ≈ maximally_mixed_dm(size_enr).data

        # general case (also test Int and BigFloat)
        nvec = BigFloat[0.123, 0.456, 0.789]
        ρTI = enr_thermal_dm(space_enr, Int64[1, 2, 3]; sparse = Val(false))
        ρTd = enr_thermal_dm(space_enr, nvec)
        ρTs = enr_thermal_dm(space_enr, nvec; sparse = Val(true))
        @test eltype(ρTI.data) == Float64
        @test isoper(ρTd)
        @test tr(ρTd) ≈ tr(ρTs) ≈ 1.0
        @test diag(ρTd) ≈ Float64[
            0.44317797863426783,
            0.19545412249437527,
            0.13879749880303996,
            0.06121365374823841,
            0.0434695463284246,
            0.019171309140931812,
            0.04854041974355738,
            0.021407708875163092,
            0.015202219370235007,
            0.006704612120243388,
            0.004761134637930744,
            0.0020997961035927096,
        ]
        @test ρTs.data isa AbstractSparseMatrix
        @test ρTd ≈ ρTs

        # tensor between normal and ENR space
        ρ_enr = enr_thermal_dm(space_enr, rand(3))
        ρ_tot = tensor(ρ_s, ρ_enr)
        opstring = sprint((t, s) -> show(t, "text/plain", s), ρ_tot)
        datastring = sprint((t, s) -> show(t, "text/plain", s), ρ_tot.data)
        ρ_tot_dims = ρ_tot.dims  # Now returns tuple format
        ρ_tot_size = size_s * size_enr
        ρ_tot_isherm = isherm(ρ_tot)
        @test opstring ==
            "\nQuantum Object:   type=Operator()   dims=$ρ_tot_dims   size=$((ρ_tot_size, ρ_tot_size))   ishermitian=$ρ_tot_isherm\n$datastring"

        # use non-square Dimensions to do partial trace
        new_dims1 = Dimensions(
            TensorSpace(Space(1), Space(1), space_enr),
            TensorSpace(Space(1), Space(1), space_enr),
        )
        ρ_enr_compound = Qobj(zeros(ComplexF64, size_enr, size_enr), dims = new_dims1)
        basis_list = [tensor(basis(D1, i), basis(D2, j)) for i in 0:(D1 - 1) for j in 0:(D2 - 1)]
        for b in basis_list
            ρ_enr_compound += tensor(b', I_enr) * ρ_tot * tensor(b, I_enr)
        end
        new_dims2 = Dimensions(
            TensorSpace(space_s..., Space(1)),
            TensorSpace(space_s..., Space(1)),
        )
        ρ_s_compound = Qobj(zeros(ComplexF64, size_s, size_s), dims = new_dims2)
        basis_list = [enr_fock(space_enr, space_enr.idx2state[idx]) for idx in 1:space_enr.size]
        for b in basis_list
            ρ_s_compound += tensor(I_s, b') * ρ_tot * tensor(I_s, b)
        end
        @test ρ_enr.data ≈ ρ_enr_compound.data
        @test ρ_s.data ≈ ρ_s_compound.data
    end

    @testset "mesolve, steadystate, and eigenstates" begin
        ε = 2π
        ωc = 2π
        g = 0.1ωc
        γ = 0.01ωc
        tlist = range(0, 20, 100)
        N_cut = 2

        # normal mesolve and steadystate
        sz = sigmaz() ⊗ qeye(N_cut)
        sm = sigmam() ⊗ qeye(N_cut)
        a = qeye(2) ⊗ destroy(N_cut)
        H_JC = 0.5ε * sz + ωc * a' * a + g * (sm' * a + a' * sm)
        c_ops_JC = (√γ * a,)
        ψ0_JC = basis(2, 0) ⊗ fock(N_cut, 0)
        sol_JC = mesolve(H_JC, ψ0_JC, tlist, c_ops_JC; e_ops = [sz], progress_bar = Val(false))
        ρ_ss_JC = steadystate(H_JC, c_ops_JC)

        # ENR mesolve and steadystate
        N_exc = 1
        dims = (2, N_cut)
        sm_enr, a_enr = enr_destroy(dims, N_exc)
        sz_enr = 2 * sm_enr' * sm_enr - 1
        ψ0_enr = enr_fock(dims, N_exc, [1, 0])
        H_enr = ε * sm_enr' * sm_enr + ωc * a_enr' * a_enr + g * (sm_enr' * a_enr + a_enr' * sm_enr)
        c_ops_enr = (√γ * a_enr,)
        sol_enr = mesolve(H_enr, ψ0_enr, tlist, c_ops_enr; e_ops = [sz_enr], progress_bar = Val(false))
        ρ_ss_enr = steadystate(H_enr, c_ops_enr)

        # check mesolve result
        @test all(isapprox.(sol_JC.expect, sol_enr.expect, atol = 1.0e-4))

        # check steadystate result
        @test expect(sz, ρ_ss_JC) ≈ expect(sz_enr, ρ_ss_enr) atol = 1.0e-4

        # check eigenstates
        λ, v = eigenstates(H_enr)
        @test all([H_enr * v[k] ≈ λ[k] * v[k] for k in eachindex(λ)])
    end

    @testset "Type Inference" begin
        N = 3
        dims = (2, 2, 3)
        excitations = 3
        @inferred enr_identity(dims, excitations)
        @inferred enr_fock(dims, excitations, zeros(Int, N))
        @inferred enr_destroy(dims, excitations)
        @inferred enr_thermal_dm(dims, excitations, rand(N))
    end
end
