function spre(O::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T}
    A = spdiagm( ones(T, size(O, 1)) )
    QuantumObject(kron(A, O.data), SuperOperatorQuantumObject, O.dims)
end

function spost(O::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T}
    B = spdiagm( ones(T, size(O, 1)) )
    QuantumObject(kron(transpose(O.data), B), SuperOperatorQuantumObject, O.dims)
end

sprepost(A::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, B::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T} = spre(A) * spost(B)

function lindblad_dissipator(O::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T}
    Od_O = O' * O
    return sprepost(O, O') - spre(Od_O) / 2 - spost(Od_O) / 2
end

destroy(N::Int) = QuantumObject(spdiagm(1 => Array{ComplexF64}(sqrt.(1:N - 1))), OperatorQuantumObject, [N])
create(N::Int) = QuantumObject(spdiagm(-1 => Array{ComplexF64}(sqrt.(1:N - 1))), OperatorQuantumObject, [N])

sigmap() = destroy(2)
sigmam() = create(2)
sigmax() = sigmam() + sigmap()
sigmay() = 1im * (sigmam() - sigmap())
sigmaz() = sigmap() * sigmam() - sigmam() * sigmap()

eye(N::Int) = QuantumObject(spdiagm(ones(ComplexF64, N)), OperatorQuantumObject, [N])

function fock(N::Int, pos::Int)
    array = zeros(N)
    array[pos + 1] = 1
    QuantumObject(Array{ComplexF64}(array), KetQuantumObject, [N])
end

basis(N::Int, pos::Int) = fock(N, pos)

function coherent(N::Real, α::T) where T <: Number
    a = destroy(N)
    ad = a'
    return exp(α * ad - α' * a) * fock(N, 0)
end

projection(N::Int, i::Int, j::Int) = fock(N, i) * fock(N, j)'

sinm(O::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T} = -0.5im * (exp(1im * O) - exp(-1im * O))

cosm(O::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T} = 0.5 * (exp(1im * O) + exp(-1im * O))