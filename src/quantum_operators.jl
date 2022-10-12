function spre(O::AbstractArray)
    A = spdiagm( ones(ComplexF64, size(O, 1)) )
    return kron(A, O)
end
 
function spost(O::AbstractArray)
    B = spdiagm( ones(ComplexF64, size(O, 1)) )
    return kron(O', B)
end
 
sprepost(A::AbstractArray, B::AbstractArray) = spre(A) * spost(B)
 
function lindblad_dissipator(O::AbstractArray)
    Od_O = O' * O
    return sprepost(O, O') - spre(Od_O) / 2 - spost(Od_O) / 2
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

function _wig_laguerre_val(L, x, c)
    if length(c) == 1
        y0 = c[1]
        y1 = 0
    elseif  length(c) == 2
        y0 = c[1]
        y1 = c[2]
    else
        k = length(c)
        y0 = c[end-1]
        y1 = c[end]
        for i in range(3, length(c), step = 1)
            k -= 1
            y0, y1 = @. c[end+1-i] - y1 * ( (k - 1)*(L+k-1)/((L + k)*k) )^0.5, y0 - y1 * ((L + 2*k -1) - x) * ((L+k)*k)^(-0.5)
        end
    end
    return @. y0 - y1 * ((L + 1) - x) * (L + 1)^(-0.5)
end

function wigner(state, xvec::Union{Vector{T1}, LinRange{T2}}, 
    yvec::Union{Vector{T1}, LinRange{T2}}; g::Real = √2) where {T1,T2}
    if length(size(state)) == 1
        ρ = state * state'
    else
        ρ = state
    end
    M = size(ρ, 1)
    X,Y = meshgrid(xvec, yvec)
    A2 = g * (X + 1im * Y)

    B = abs.(A2)
    B .*= B
    w0 = (2*ρ[1, end]) .* ones(eltype(A2), size(A2)...)
    L = M-1

    ρ = ρ .* (2 * ones(M,M) - diagm(ones(M)))
    while L > 0
        L -= 1
        w0 = _wig_laguerre_val(L, B, diag(ρ, L)) .+ w0 .* A2 .* (L+1)^(-0.5)
    end

    return @. real(w0) * exp(-B * 0.5) * (g*g*0.5 / π)
end