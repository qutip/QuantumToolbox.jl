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
 
sprepost(A::AbstractArray, B::AbstractArray) = spre(A) * spost(B)
 
function lindblad_dissipator(O::AbstractArray)
    Od_O = adjoint(O) * O
    return sprepost(O, adjoint(O)) - 0.5 * spre(Od_O) - 0.5 * spost(Od_O)
end

destroy(N::Number) = spdiagm(1 => Array{ComplexF64}(sqrt.(1:N - 1)))

create(N::Number) = spdiagm(-1 => Array{ComplexF64}(sqrt.(1:N - 1)))

sigmap() = destroy(2)

sigmam() = create(2)

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
        if length(size(state)) == 1
            return real(adjoint(state) * op * state)
        else
            return real(tr(op * state))
        end
    else
        if length(size(state)) == 1
            return adjoint(state) * op * state
        else
            return tr(op * state)
        end
    end
end

function wigner(state, xvec, yvec; g = sqrt(2))
    if length(size(state)) == 1
        ρ = state * state'
    else
        ρ = state
    end
    M = size(ρ, 1)
    X, Y = meshgrid(xvec, yvec)
    A = @. 0.5 * g * (X + 1.0im * Y)

    Wlist = [zeros(ComplexF64, size(A)) for k in 1:M]
    Wlist[1] = @. exp(-2.0 * abs(A)^2) / pi

    W = @. real(ρ[1,1]) * real(Wlist[1])
    for n in 2:M
        Wlist[n] = @. (2.0 * A * Wlist[n - 1]) / sqrt(n)
        W += @. 2 * real(ρ[1, n] * Wlist[n])
    end

    for m in 2:M
        temp = deepcopy(Wlist[m])
        Wlist[m] = @. (2 * conj(A) * temp - sqrt(m) * Wlist[m - 1]) / sqrt(m)

        # Wlist[m] = Wigner function for |m><m|
        W += @. real(ρ[m, m] * Wlist[m])

        for n in m+1:M
            temp2 = @. (2 * A * Wlist[n - 1] - sqrt(m) * temp) / sqrt(n)
            temp = deepcopy(Wlist[n])
            Wlist[n] = temp2

            # Wlist[n] = Wigner function for |m><n|
            W += @. 2 * real(ρ[m, n] * Wlist[n])
        end
    end

    return 0.5 * W * g^2
end