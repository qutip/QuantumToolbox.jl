function benchmark_autodiff!(SUITE)
    # Use harmonic oscillator system for both sesolve and mesolve
    N = 20
    a = destroy(N)
    ψ0 = fock(N, 0)
    tlist = range(0, 40, 100)

    # ---- SESOLVE ----
    # For direct Forward differentiation
    function my_f_sesolve_direct(p)
        H = p[1] * a' * a + p[2] * (a + a')
        sol = sesolve(H, ψ0, tlist, progress_bar = Val(false))
        return real(expect(a' * a, sol.states[end]))
    end

    # For SciMLSensitivity.jl (reverse mode with Zygote and Enzyme)
    coef_Δ(p, t) = p[1]
    coef_F(p, t) = p[2]
    H_evo = QobjEvo(a' * a, coef_Δ) + QobjEvo(a + a', coef_F)

    function my_f_sesolve(p)
        sol = sesolve(
            H_evo,
            ψ0,
            tlist,
            progress_bar = Val(false),
            params = p,
            sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
        )
        return real(expect(a' * a, sol.states[end]))
    end

    # ---- MESOLVE ----
    # For direct Forward differentiation
    function my_f_mesolve_direct(p)
        H = p[1] * a' * a + p[2] * (a + a')
        c_ops = [sqrt(p[3]) * a]
        sol = mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
        return real(expect(a' * a, sol.states[end]))
    end

    # For SciMLSensitivity.jl (reverse mode with Zygote and Enzyme)
    coef_γ(p, t) = sqrt(p[3])
    c_ops = [QobjEvo(a, coef_γ)]
    L = liouvillian(H_evo, c_ops)

    function my_f_mesolve(p)
        sol = mesolve(
            L,
            ψ0,
            tlist,
            progress_bar = Val(false),
            params = p,
            sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
        )
        return real(expect(a' * a, sol.states[end]))
    end

    # Parameters for benchmarks
    params_sesolve = [1.0, 1.0]
    params_mesolve = [1.0, 1.0, 1.0]

    # Benchmark sesolve - Forward
    SUITE["Autodiff"]["sesolve"]["Forward"] = @benchmarkable ForwardDiff.gradient($my_f_sesolve_direct, $params_sesolve)

    # Benchmark sesolve - Reverse (Zygote)
    SUITE["Autodiff"]["sesolve"]["Reverse (Zygote)"] = @benchmarkable Zygote.gradient($my_f_sesolve, $params_sesolve)

    # Benchmark sesolve - Reverse (Enzyme)
    SUITE["Autodiff"]["sesolve"]["Reverse (Enzyme)"] = @benchmarkable Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Const($my_f_sesolve),
        Active,
        Duplicated($params_sesolve, dparams_sesolve),
    ) setup=(dparams_sesolve = Enzyme.make_zero($params_sesolve))

    # Benchmark mesolve - Forward
    SUITE["Autodiff"]["mesolve"]["Forward"] = @benchmarkable ForwardDiff.gradient($my_f_mesolve_direct, $params_mesolve)

    # Benchmark mesolve - Reverse (Zygote)
    SUITE["Autodiff"]["mesolve"]["Reverse (Zygote)"] = @benchmarkable Zygote.gradient($my_f_mesolve, $params_mesolve)

    # Benchmark mesolve - Reverse (Enzyme)
    SUITE["Autodiff"]["mesolve"]["Reverse (Enzyme)"] = @benchmarkable Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Const($my_f_mesolve),
        Active,
        Duplicated($params_mesolve, dparams_mesolve),
    ) setup=(dparams_mesolve = Enzyme.make_zero($params_mesolve))

    return nothing
end
