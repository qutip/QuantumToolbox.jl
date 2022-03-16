function spre(O::AbstractArray)
    N = size(O, 1)
    A = spdiagm( ones(ComplexF64, N) )
    B = O
    return kron(A, B)
end
 
function spost(O::AbstractArray)
    N = size(O, 1)
    A = sparse( transpose(O) )
    B = spdiagm( ones(ComplexF64, N) )
    return kron(A, B)
end
 
function sprepost(A, B::AbstractArray)
    return sparse( spre(A) * spost(B) )
end
 
function lindblad_dissipator(O::AbstractArray)
    Od_O = adjoint(O) * O
    return sprepost(O, adjoint(O)) - 0.5 * spre(Od_O) - 0.5 * spost(Od_O)
end

function destroy(N::Number)
    return spdiagm(1 => Array{ComplexF64}(sqrt.(1:N - 1)))
end

function eye(N::Number)
    return spdiagm(ones(ComplexF64, N))
end

function fock(N::Number, pos::Number)
    array = zeros(N)
    array[pos + 1] = 1
    return Array{ComplexF64}( array )
end

function projection(N::Number, i::Number, j::Number, shift = false)
    if shift
        return sparse( fock(N, i - 1) * adjoint(fock(N, j - 1)) )
    else
        return sparse( fock(N, i) * adjoint(fock(N, j)) )
    end
end

# function exp_gpu(M)
#     M_tmp = copy(M)
#     E_gpu, U_gpu = CUDA.CUSOLVER.heevd!('V','U', M_tmp)
#     return U_gpu * Diagonal(exp.(E_gpu)) * adjoint(U_gpu)
# end

function sinm(O::AbstractArray)
    M = collect(O)
    return sparse( 0.5im * (exp(1im * M) - exp(-1im * M)) )
end

function cosm(O::AbstractArray)
    M = collect(O)
    return sparse( 0.5 * (exp(1im * M) + exp(-1im * M)) )
end

function expect(op::AbstractArray, state::AbstractArray)
    if ishermitian(op)
        return real(adjoint(state) * op * state)
    else
        return adjoint(state) * op * state
    end
end