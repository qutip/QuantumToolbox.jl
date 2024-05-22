function benchmark_steadystate!(SUITE)
    N = 50
    Δ = 0.1
    F = 2
    γ = 1
    nth = 5

    a = destroy(N)
    H = Δ * a' * a + F * (a + a')
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a']

    SUITE["Steadystate"]["Direct"] = @benchmarkable steadystate($H, $c_ops)

    return nothing
end
