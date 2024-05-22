### DYNAMICAL FOCK DIMENSION ###
@testset "Dynamical Fock Dimension" begin
    F, Δ, κ = 5, 0.25, 1
    t_l = range(0, 15, length = 100)

    N0 = 140
    a0 = destroy(N0)
    H0 = Δ * a0' * a0 + F * (a0 + a0')
    c_ops0 = [√κ * a0]
    e_ops0 = [a0' * a0]
    ψ00 = fock(N0, 0)
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops = e_ops0, progress_bar = false)

    function H_dfd0(dims, p)
        Δ = p.Δ
        F = p.F
        a = destroy(dims[1])
        return Δ * a' * a + F * (a + a')
    end
    function c_ops_dfd0(dims, p)
        κ = p.κ
        a = destroy(dims[1])
        return [√κ * a]
    end
    function e_ops_dfd0(dims, p)
        a = destroy(dims[1])
        return [a' * a]
    end
    maxdims = [150]
    ψ0 = fock(3, 0)
    dfd_params = (Δ = Δ, F = F, κ = κ)
    sol = dfd_mesolve(H_dfd0, ψ0, t_l, c_ops_dfd0, maxdims, dfd_params, e_ops = e_ops_dfd0, progress_bar = false)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) < 0.01

    ######################

    F = 0
    H0 = Δ * a0' * a0 + F * (a0 + a0')
    c_ops0 = [√κ * a0]
    e_ops0 = [a0' * a0]
    ψ00 = fock(N0, 50)
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops = e_ops0, progress_bar = false)

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
    maxdims = [150]
    ψ0 = fock(70, 50)
    dfd_params = (Δ = Δ, F = F, κ = κ)
    sol = dfd_mesolve(H_dfd1, ψ0, t_l, c_ops_dfd1, maxdims, dfd_params, e_ops = e_ops_dfd1, progress_bar = false)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) < 0.01

    ######################

    F, Δ, κ, J = 1.5, 0.25, 1, 0.05
    N0 = 25
    N1 = 20
    a0 = kron(destroy(N0), qeye(N1))
    a1 = kron(qeye(N0), destroy(N1))
    H0 = Δ * a0' * a0 + F * (a0 + a0') + Δ * a1' * a1 + J * (a0' * a1 + a0 * a1')
    c_ops0 = [√κ * a0, √κ * a1]
    e_ops0 = [a0' * a0, a1' * a1]
    ψ00 = kron(fock(N0, 0), fock(N1, 15))
    sol0 = mesolve(H0, ψ00, t_l, c_ops0, e_ops = e_ops0, progress_bar = false)

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
    maxdims = [50, 50]
    ψ0 = kron(fock(3, 0), fock(20, 15))
    dfd_params = (Δ = Δ, F = F, κ = κ, J = J)
    sol = dfd_mesolve(H_dfd2, ψ0, t_l, c_ops_dfd2, maxdims, dfd_params, e_ops = e_ops_dfd2, progress_bar = false)

    @test sum(abs.((sol.expect[1, :] .- sol0.expect[1, :]) ./ (sol0.expect[1, :] .+ 1e-16))) +
          sum(abs.((sol.expect[2, :] .- sol0.expect[2, :]) ./ (sol0.expect[2, :] .+ 1e-16))) < 0.01
end
