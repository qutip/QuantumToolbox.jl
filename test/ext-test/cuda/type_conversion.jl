using Test
using QuantumToolbox
using CUDA
using CUDA.cuSPARSE

@testset "CUDA (type conversion)" verbose = true begin
    # Test that scalar indexing is disallowed
    @test_throws ErrorException CUDA.rand(1)[1]

    ψdi = Qobj(Int64[1, 0])
    ψdf = Qobj(Float64[1, 0])
    ψdc = Qobj(ComplexF64[1, 0])
    ψsi = to_sparse(ψdi)
    ψsf = to_sparse(ψdf)
    ψsc = to_sparse(ψdc)

    Xdi = Qobj(Int64[0 1; 1 0])
    Xdf = Qobj(Float64[0 1; 1 0])
    Xdc = Qobj(ComplexF64[0 1; 1 0])
    Xsi = to_sparse(Xdi)
    Xsf = to_sparse(Xdf)
    Xsc = to_sparse(Xdc)

    @test_throws DomainError cu(ψdi; word_size = 16)

    # type conversion of CUDA dense arrays
    @test typeof(cu(ψdi; word_size = 64).data) == typeof(CuArray(ψdi).data) <: CuArray{Int64, 1}
    @test typeof(cu(ψdi; word_size = 32).data) == typeof(CuArray{Int32}(ψdi).data) <: CuArray{Int32, 1}
    @test typeof(cu(ψdf; word_size = 64).data) == typeof(CuArray(ψdf).data) <: CuArray{Float64, 1}
    @test typeof(cu(ψdf; word_size = 32).data) == typeof(CuArray{Float32}(ψdf).data) <: CuArray{Float32, 1}
    @test typeof(cu(ψdc; word_size = 64).data) == typeof(CuArray(ψdc).data) <: CuArray{ComplexF64, 1}
    @test typeof(cu(ψdc; word_size = 32).data) == typeof(CuArray{ComplexF32}(ψdc).data) <: CuArray{ComplexF32, 1}
    @test typeof(cu(Xdi; word_size = 64).data) == typeof(CuArray(Xdi).data) <: CuArray{Int64, 2}
    @test typeof(cu(Xdi; word_size = 32).data) == typeof(CuArray{Int32}(Xdi).data) <: CuArray{Int32, 2}
    @test typeof(cu(Xdf; word_size = 64).data) == typeof(CuArray(Xdf).data) <: CuArray{Float64, 2}
    @test typeof(cu(Xdf; word_size = 32).data) == typeof(CuArray{Float32}(Xdf).data) <: CuArray{Float32, 2}
    @test typeof(cu(Xdc; word_size = 64).data) == typeof(CuArray(Xdc).data) <: CuArray{ComplexF64, 2}
    @test typeof(cu(Xdc; word_size = 32).data) == typeof(CuArray{ComplexF32}(Xdc).data) <: CuArray{ComplexF32, 2}

    # type conversion of CUDA sparse arrays
    @test typeof(cu(ψsi; word_size = 64).data) == typeof(CuSparseVector(ψsi).data) == CuSparseVector{Int64, Int32}
    @test typeof(cu(ψsi; word_size = 32).data) == typeof(CuSparseVector{Int32}(ψsi).data) == CuSparseVector{Int32, Int32}
    @test typeof(cu(ψsf; word_size = 64).data) == typeof(CuSparseVector(ψsf).data) == CuSparseVector{Float64, Int32}
    @test typeof(cu(ψsf; word_size = 32).data) ==
        typeof(CuSparseVector{Float32}(ψsf).data) ==
        CuSparseVector{Float32, Int32}
    @test typeof(cu(ψsc; word_size = 64).data) == typeof(CuSparseVector(ψsc).data) == CuSparseVector{ComplexF64, Int32}
    @test typeof(cu(ψsc; word_size = 32).data) ==
        typeof(CuSparseVector{ComplexF32}(ψsc).data) ==
        CuSparseVector{ComplexF32, Int32}
    @test typeof(cu(Xsi; word_size = 64).data) == typeof(CuSparseMatrixCSC(Xsi).data) == CuSparseMatrixCSC{Int64, Int32}
    @test typeof(cu(Xsi; word_size = 32).data) ==
        typeof(CuSparseMatrixCSC{Int32}(Xsi).data) ==
        CuSparseMatrixCSC{Int32, Int32}
    @test typeof(cu(Xsf; word_size = 64).data) ==
        typeof(CuSparseMatrixCSC(Xsf).data) ==
        CuSparseMatrixCSC{Float64, Int32}
    @test typeof(cu(Xsf; word_size = 32).data) ==
        typeof(CuSparseMatrixCSC{Float32}(Xsf).data) ==
        CuSparseMatrixCSC{Float32, Int32}
    @test typeof(cu(Xsc; word_size = 64).data) ==
        typeof(CuSparseMatrixCSC(Xsc).data) ==
        CuSparseMatrixCSC{ComplexF64, Int32}
    @test typeof(cu(Xsc; word_size = 32).data) ==
        typeof(CuSparseMatrixCSC{ComplexF32}(Xsc).data) ==
        CuSparseMatrixCSC{ComplexF32, Int32}
    @test typeof(CuSparseMatrixCSR(Xsi).data) == CuSparseMatrixCSR{Int64, Int32}
    @test typeof(CuSparseMatrixCSR{Int32}(Xsi).data) == CuSparseMatrixCSR{Int32, Int32}
    @test typeof(CuSparseMatrixCSR(Xsf).data) == CuSparseMatrixCSR{Float64, Int32}
    @test typeof(CuSparseMatrixCSR{Float32}(Xsf).data) == CuSparseMatrixCSR{Float32, Int32}
    @test typeof(CuSparseMatrixCSR(Xsc).data) == CuSparseMatrixCSR{ComplexF64, Int32}
    @test typeof(CuSparseMatrixCSR{ComplexF32}(Xsc).data) == CuSparseMatrixCSR{ComplexF32, Int32}

    # Sparse To Dense
    # @test to_dense(cu(ψsi; word_size = 64)).data isa CuVector{Int64} # TODO: Fix this in CUDA.jl
    @test to_dense(cu(ψsf; word_size = 64)).data isa CuVector{Float64}
    @test to_dense(cu(ψsc; word_size = 64)).data isa CuVector{ComplexF64}
    # @test to_dense(cu(Xsi; word_size = 64)).data isa CuMatrix{Int64} # TODO: Fix this in CUDA.jl
    @test to_dense(cu(Xsf; word_size = 64)).data isa CuMatrix{Float64}
    @test to_dense(cu(Xsc; word_size = 64)).data isa CuMatrix{ComplexF64}

    # @test to_dense(Int32, cu(ψsf; word_size = 64)).data isa CuVector{Int32} # TODO: Fix this in CUDA.jl
    # @test to_dense(Float32, cu(ψsf; word_size = 64)).data isa CuVector{Float32} # TODO: Fix this in CUDA.jl
    # @test to_dense(ComplexF32, cu(ψsf; word_size = 64)).data isa CuVector{ComplexF32} # TODO: Fix this in CUDA.jl
    # @test to_dense(Int64, cu(Xsf; word_size = 32)).data isa CuMatrix{Int64} # TODO: Fix this in CUDA.jl
    # @test to_dense(Float64, cu(Xsf; word_size = 32)).data isa CuMatrix{Float64} # TODO: Fix this in CUDA.jl
    # @test to_dense(ComplexF64, cu(Xsf; word_size = 32)).data isa CuMatrix{ComplexF64} # TODO: Fix this in CUDA.jl
end
