using Revise
using LinearAlgebra
using SparseArrays
using SciMLOperators
using QuantumToolbox
using QuantumToolboxOperators
using Chairmarks

function to_profile(y, A, x)
    for _ in 1:100
        mul!(y, A, x)
    end
    return y
end

# %% -------- TEST with N SPINS -----------

T = ComplexF32
N = 20
dims = ntuple(_ -> 2, Val(N))

σx_data = T[0 1; 1 0]
σz_data = T[1 0; 0 -1]

σx = sigmax()
σz = sigmaz()

ψ = randn(T, 2^N) |> normalize
dψ = similar(ψ)

H_full = sum(1:(N - 1)) do n
    multisite_operator(dims, n => σx, n + 1 => σx).data
end

H_tensor = sum(1:(N - 1)) do n
    prova = ntuple(Val(N)) do j
        if j == n || j == n + 1
            MatrixOperator(σx_data)
        else
            IdentityOperator(dims[j])
        end
    end
    TensorProductOperator(prova...)
end
H_tensor = cache_operator(H_tensor, ψ)

H_local = sum(1:(N - 1)) do n
    LocalTensorProductOperator(dims, n => MatrixOperator(σx_data), n + 1 => MatrixOperator(σx_data))
end
H_local_cache = cache_operator(H_local.ops[1], ψ).cache  # Cache the first term (all are the same)
H_local = SciMLOperators.AddedOperator((LocalTensorProductOperator(op.dims, op.indices, op.ops, H_local_cache) for op in H_local.ops)...)

mul!(dψ, H_full, ψ)
mul!(dψ, H_tensor, ψ)
mul!(dψ, H_local, ψ)

memory_ratio_tensor = Base.summarysize(H_full) / Base.summarysize(H_tensor)
memory_ratio_local = Base.summarysize(H_full) / Base.summarysize(H_local)

# %%

@be mul!($dψ, $H_full, $ψ)
@be mul!($dψ, $H_tensor, $ψ)
@be mul!($dψ, $H_local, $ψ)


# %% -------- TEST with N CAVITIES -----------

T = ComplexF32
N = 5
dims = ntuple(_ -> 30, Val(N))

ψ = randn(T, prod(dims)) |> normalize
dψ = similar(ψ)

H_full = sum(1:(N - 1)) do n
    multisite_operator(dims, n => create(dims[n]), n + 1 => destroy(dims[n + 1])).data
end

H_tensor = sum(1:(N - 1)) do n
    prova = ntuple(Val(N)) do j
        if j == n
            DestroyOperator{T}(dims[j])'
        elseif j == n + 1
            DestroyOperator{T}(dims[j])
        else
            IdentityOperator(dims[j])
        end
    end
    TensorProductOperator(prova...)
end
H_tensor = cache_operator(H_tensor, ψ)

H_local = sum(1:(N - 1)) do n
    LocalTensorProductOperator(dims, n => DestroyOperator{T}(dims[n])', n + 1 => DestroyOperator{T}(dims[n + 1]))
end
H_local_cache = cache_operator(H_local.ops[1], ψ).cache  # Cache the first term (all are the same)
H_local = SciMLOperators.AddedOperator((LocalTensorProductOperator(op.dims, op.indices, op.ops, H_local_cache) for op in H_local.ops)...)

mul!(dψ, H_full, ψ)
mul!(dψ, H_tensor, ψ)
mul!(dψ, H_local, ψ)

Base.summarysize(H_local) / Base.summarysize(ψ)

memory_ratio_tensor = Base.summarysize(H_full) / Base.summarysize(H_tensor)
memory_ratio_local = Base.summarysize(H_full) / Base.summarysize(H_local)

# %%

@be mul!($dψ, $H_full, $ψ)
@be mul!($dψ, $H_tensor, $ψ)
@be mul!($dψ, $H_local, $ψ)

# %%

using Reactant

ψ_reactant = Reactant.to_rarray(ψ)
dψ_reactant = similar(ψ_reactant)

H_local_cache_reactant = Reactant.to_rarray(H_local_cache) # cache_operator(H_local.ops[1], ψ).cache  # Cache the first term (all are the same)
H_local_reactant = SciMLOperators.AddedOperator((LocalTensorProductOperator(op.dims, op.indices, op.ops, H_local_cache_reactant) for op in H_local.ops)...)

mul_compiled! = @compile mul!(dψ_reactant, H_local_reactant, ψ_reactant)

mul_compiled!(dψ_reactant, H_local_reactant, ψ_reactant)

# %%

using Revise
using SciMLOperators

m1, n1 = 3, 5
m2, n2 = 7, 11
m3, n3 = 13, 17

Id1 = IdentityOperator(m1)
Id2 = IdentityOperator(m2)
Id3 = IdentityOperator(m3)
A1 = MatrixOperator(rand(n1, n1))
A2 = MatrixOperator(rand(n2, n2))

op1 = kron(A1, Id1, Id2, Id3)
op2 = kron(Id1, A1, Id2, Id3)
op3 = kron(Id1, Id2, A1, Id3)
op4 = kron(Id1, Id2, Id3, A1)

op5 = kron(A1, A2, Id1, Id2)
op6 = kron(Id1, A1, A2, Id2)
op7 = kron(Id1, Id2, A1, A2)

# Test the structure of the resulting operators
# The nesting structure should be deep 2 at most
op1.ops[1] isa MatrixOperator
op1.ops[2] isa IdentityOperator

op2.ops[1] isa TensorProductOperator
op2.ops[2] isa IdentityOperator
op2.ops[1].ops[1] isa IdentityOperator
op2.ops[1].ops[2] isa MatrixOperator

op3.ops[1] isa TensorProductOperator
op3.ops[2] isa IdentityOperator
op3.ops[1].ops[1] isa IdentityOperator
op3.ops[1].ops[2] isa MatrixOperator

op4.ops[1] isa IdentityOperator
op4.ops[2] isa MatrixOperator

op5.ops[1] isa MatrixOperator
op5.ops[2] isa TensorProductOperator
op5.ops[2].ops[1] isa MatrixOperator
op5.ops[2].ops[2] isa IdentityOperator

op6.ops[1] isa TensorProductOperator
op6.ops[2] isa TensorProductOperator
op6.ops[1].ops[1] isa IdentityOperator
op6.ops[1].ops[2] isa MatrixOperator
op6.ops[2].ops[1] isa MatrixOperator
op6.ops[2].ops[2] isa IdentityOperator

op7.ops[1] isa TensorProductOperator
op7.ops[2] isa MatrixOperator
op7.ops[1].ops[1] isa IdentityOperator
op7.ops[1].ops[2] isa MatrixOperator
