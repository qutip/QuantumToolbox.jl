@testset "Metal Extension" verbose = true begin
    ψdi = Qobj(Int64[1, 0])
    ψdf = Qobj(Float64[1, 0])
    ψdc = Qobj(ComplexF64[1, 0])
    ψsi = dense_to_sparse(ψdi)
    ψsf = dense_to_sparse(ψdf)
    ψsc = dense_to_sparse(ψdc)

    Xdi = Qobj(Int64[0 1; 1 0])
    Xdf = Qobj(Float64[0 1; 1 0])
    Xdc = Qobj(ComplexF64[0 1; 1 0])
    Xsi = dense_to_sparse(Xdi)
    Xsf = dense_to_sparse(Xdf)
    Xsc = dense_to_sparse(Xdc)

    # type conversion of dense arrays
    @test typeof(mtl(ψdi).data) == typeof(MtlArray{Int32}(ψdi).data) <: MtlVector{Int32}
    @test typeof(mtl(ψdf).data) ==
          typeof(MtlArray(ψdf).data) ==
          typeof(MtlArray{Float32}(ψdf).data) <:
          MtlVector{Float32}
    @test typeof(mtl(ψdc).data) ==
          typeof(MtlArray(ψdc).data) ==
          typeof(MtlArray{ComplexF32}(ψdc).data) <:
          MtlVector{ComplexF32}
    @test typeof(mtl(Xdi).data) == typeof(MtlArray{Int32}(Xdi).data) <: MtlMatrix{Int32}
    @test typeof(mtl(Xdf).data) ==
          typeof(MtlArray(Xdf).data) ==
          typeof(MtlArray{Float32}(Xdf).data) <:
          MtlMatrix{Float32}
    @test typeof(mtl(Xdc).data) ==
          typeof(MtlArray(Xdc).data) ==
          typeof(MtlArray{ComplexF32}(Xdc).data) <:
          MtlMatrix{ComplexF32}

    # type conversion of sparse arrays
    @test typeof(mtl(ψsi).data) == typeof(MtlArray{Int32}(ψsi).data) <: MtlVector{Int32}
    @test typeof(mtl(ψsf).data) ==
          typeof(MtlArray(ψsf).data) ==
          typeof(MtlArray{Float32}(ψsf).data) <:
          MtlVector{Float32}
    @test typeof(mtl(ψsc).data) ==
          typeof(MtlArray(ψsc).data) ==
          typeof(MtlArray{ComplexF32}(ψsc).data) <:
          MtlVector{ComplexF32}
    @test typeof(mtl(Xsi).data) == typeof(MtlArray{Int32}(Xsi).data) <: MtlMatrix{Int32}
    @test typeof(mtl(Xsf).data) ==
          typeof(MtlArray(Xsf).data) ==
          typeof(MtlArray{Float32}(Xsf).data) <:
          MtlMatrix{Float32}
    @test typeof(mtl(Xsc).data) ==
          typeof(MtlArray(Xsc).data) ==
          typeof(MtlArray{ComplexF32}(Xsc).data) <:
          MtlMatrix{ComplexF32}

    # brief example in README and documentation
    N = 5 # cannot be too large since Metal.jl does not support sparse matrix
    ω = 1.0f0  # Float32
    γ = 0.1f0  # Float32
    tlist = range(0, 10, 100)

    ## calculate by CPU
    a_cpu = destroy(N)
    ψ0_cpu = fock(N, 3)
    H_cpu = ω * a_cpu' * a_cpu
    sol_cpu = mesolve(H_cpu, ψ0_cpu, tlist, [sqrt(γ) * a_cpu], e_ops = [a_cpu' * a_cpu], progress_bar = Val(false))

    ## calculate by GPU
    a_gpu = mtl(destroy(N))
    ψ0_gpu = mtl(fock(N, 3))
    H_gpu = ω * a_gpu' * a_gpu
    sol_gpu = mesolve(H_gpu, ψ0_gpu, tlist, [sqrt(γ) * a_gpu], e_ops = [a_gpu' * a_gpu], progress_bar = Val(false))

    @test all(isapprox.(sol_cpu.expect, sol_gpu.expect; atol = 1e-6))
end
