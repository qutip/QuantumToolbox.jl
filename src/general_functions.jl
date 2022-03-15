function chop_op(O, tol = 1e-8)
    tmp_r = (abs.(real.(O)) .> tol) .* real.(O)
    tmp_i = (abs.(imag.(O)) .> tol) .* imag.(O)
    return tmp_r + 1im .* tmp_i
end

function gaussian(x, mu, sig)
    return exp.(- 0.5 * (x .- mu).^2 / sig^2)
end

function gaussian_derivative(x, mu, sig)
    return - (x .- mu) ./ sig^2 .* exp.(- 0.5 * (x .- mu).^2 / sig^2)
end

function trunc_op(op, states)
    N_trunc = size(states)[1]
    # qstates = [qtp.Qobj(states[i], dims = [[N_s, N_s], [1, 1]]) for i in range(N_trunc)]
    res = spzeros(N_trunc, N_trunc)
    for i in range(1, N_trunc, step = 1)
        for j in range(1, N_trunc, step = 1)
            res += (adjoint(states[i, :]) * op * states[j, :]) * projection(N_trunc, i, j, true)
        end
    end
    return chop_op(res)
end