function run_all_functions()
    # block diagonal form
    run_block_diagonal_form()

    # correlations and spectrum
    run_correlations_and_spectrum()

    # dynamical fock dimension
    # run_dynamical_fock_dimension() # TODO: fix type instabilities here

    # dynamical shifted fock
    # run_dynamical_shifted_fock() # TODO: fix type instabilities here

    return nothing
end

#=
    Block Diagonal Form
=#
function run_block_diagonal_form()
    H, c_ops, a = driven_dissipative_kerr()
    L = liouvillian(H, c_ops)

    bdf = block_diagonal_form(L)

    return nothing
end

#=
    Correlations and Spectrum
=#
function run_correlations_and_spectrum()
    N = 10
    H, c_ops, a = driven_dissipative_harmonic_oscillator(nth = 0.01, N = N)

    t_l = range(0, 333 * π, length = 1000)
    corr1 = correlation_2op_1t(H, nothing, t_l, c_ops, a', a; progress_bar = Val(false))
    corr2 = correlation_3op_1t(H, nothing, t_l, c_ops, qeye(a.dims[1]), a', a; progress_bar = Val(false))
    ω_l1, spec1 = spectrum_correlation_fft(t_l, corr1)

    ω_l2 = range(0, 3, length = 1000)
    spec2 = spectrum(H, ω_l2, c_ops, a', a)
    # spec3 = spectrum(H, ω_l2, c_ops, a', a; solver = PseudoInverse())

    return nothing
end

#=
    Dynamical Fock Dimension
=#
function run_dynamical_fock_dimension()
    t_l = range(0, 15, length = 100)

    # Single mode
    F, Δ, κ = 5, 0.25, 1

    maxdims = [150]
    ψ0 = fock(3, 0)
    dfd_params = (Δ = Δ, F = F, κ = κ)
    sol = dfd_mesolve(H_dfd1, ψ0, t_l, c_ops_dfd1, maxdims, dfd_params, e_ops = e_ops_dfd1, progress_bar = Val(false))

    # Two modes
    F, Δ, κ, J = 1.5, 0.25, 1, 0.05

    maxdims = [50, 50]
    ψ0 = kron(fock(3, 0), fock(20, 15))
    dfd_params = (Δ = Δ, F = F, κ = κ, J = J)
    sol = dfd_mesolve(H_dfd2, ψ0, t_l, c_ops_dfd2, maxdims, dfd_params, e_ops = e_ops_dfd2, progress_bar = Val(false))

    return nothing
end

#=
    Dynamical Shifted Fock
=#
function run_dynamical_shifted_fock()
    # Single mode
    F = 3
    Δ = 0.25
    κ = 1
    U = 0.01
    α0 = 1.5

    tlist = range(0, 25, 300)

    N = 5
    a = destroy(N)

    op_list = [a]
    ψ0 = fock(N, 0)
    α0_l = [α0]
    dsf_params = (Δ = Δ, F = F, κ = κ, U = U)

    sol_dsf_me = dsf_mesolve(
        H_dsf,
        ψ0,
        tlist,
        c_ops_dsf,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf,
        progress_bar = Val(false),
    )
    sol_dsf_mc = dsf_mcsolve(
        H_dsf,
        ψ0,
        tlist,
        c_ops_dsf,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf,
        progress_bar = Val(false),
        ntraj = 500,
    )

    # Two modes
    F = 2
    Δ = 0.25
    κ = 1
    U = 0.01
    J = 0.5
    tlist = range(0, 15, 300)

    N = 5
    a1 = kron(destroy(N), qeye(N))
    a2 = kron(qeye(N), destroy(N))

    op_list = [a1, a2]
    ψ0 = kron(fock(N, 0), fock(N, 0))
    α0_l = [α0, α0]
    dsf_params = (Δ = Δ, F = F, κ = κ, U = U, J = J)

    sol_dsf_me = dsf_mesolve(
        H_dsf2,
        ψ0,
        tlist,
        c_ops_dsf2,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf2,
        progress_bar = Val(false),
    )
    sol_dsf_mc = dsf_mcsolve(
        H_dsf2,
        ψ0,
        tlist,
        c_ops_dsf2,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf2,
        progress_bar = Val(false),
        ntraj = 500,
    )

    return nothing
end
