function driven_dissipative_harmonic_oscillator(; Δ = 1, γ = 0.1, nth = 0, F = 0, N = 10)
    a = destroy(N)
    H = Δ * a' * a + F * (a + a')
    c_ops = (sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a')

    return H, c_ops, a
end

function driven_dissipative_kerr()
    N = 20
    Δ = 0
    G = 5
    tg = 0
    θ = atan(tg)
    U = sin(θ)
    κ2 = cos(θ)
    κϕ = 1e-3

    a = destroy(N)
    ad = create(N)
    H = -Δ * ad * a + G / 2 * (ad^2 + a^2) + U / 2 * (ad^2 * a^2)
    c_ops = (√(κ2) * a^2, √(κϕ) * ad * a)

    return H, c_ops, a
end

#=
    Dynamical Fock Dimension
=#
function H_dfd1(dims, p)
    Δ = p.Δ
    F = p.F
    a = destroy(dims[1])
    return Δ * a' * a + F * (a + a')
end
function c_ops_dfd1(dims, p)
    κ = p.κ
    a = destroy(dims[1])
    return [√κ * a]
end
function e_ops_dfd1(dims, p)
    a = destroy(dims[1])
    return [a' * a]
end
function H_dfd2(dims, p)
    Δ = p.Δ
    F = p.F
    J = p.J
    a = kron(destroy(dims[1]), qeye(dims[2]))
    b = kron(qeye(dims[1]), destroy(dims[2]))
    return Δ * a' * a + F * (a + a') + Δ * b' * b + J * (a' * b + a * b')
end
function c_ops_dfd2(dims, p)
    κ = p.κ
    a = kron(destroy(dims[1]), qeye(dims[2]))
    b = kron(qeye(dims[1]), destroy(dims[2]))
    return [√κ * a, √κ * b]
end
function e_ops_dfd2(dims, p)
    a = kron(destroy(dims[1]), qeye(dims[2]))
    b = kron(qeye(dims[1]), destroy(dims[2]))
    return [a' * a, b' * b]
end

#=
    Dynamical Shifted Fock
=#
function H_dsf(op_list, p)
    Δ = p.Δ
    F = p.F
    U = p.U
    a = op_list[1]
    return Δ * a' * a + F * (a + a') + U * a'^2 * a^2
end
function c_ops_dsf(op_list, p)
    κ = p.κ
    a = op_list[1]
    return [√κ * a]
end
function e_ops_dsf(op_list, p)
    a = op_list[1]
    return [a' * a, a]
end
function H_dsf2(op_list, p)
    Δ = p.Δ
    F = p.F
    U = p.U
    J = p.J
    a1, a2 = op_list
    return Δ * a1' * a1 +
           Δ * a2' * a2 +
           U * a1'^2 * a1^2 +
           U * a2'^2 * a2^2 +
           F * (a1 + a1') +
           J * (a1' * a2 + a1 * a2')
end
function c_ops_dsf2(op_list, p)
    κ = p.κ
    a1, a2 = op_list
    return [√κ * a1, √κ * a2]
end
function e_ops_dsf2(op_list, p)
    a1, a2 = op_list
    return [a1' * a1, a2' * a2]
end
