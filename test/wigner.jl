@testset "Wigner" begin
    α = 0.5 + 0.8im
    ψ = coherent(30, α)
    ρ = dense_to_sparse(ket2dm(ψ), 1e-6)
    xvec = LinRange(-3, 3, 300)
    yvec = LinRange(-3, 3, 300)

    wig = wigner(ψ, xvec, yvec, solver = WignerLaguerre(tol = 1e-6))
    wig2 = wigner(ρ, xvec, yvec, solver = WignerLaguerre(parallel = false))
    wig3 = wigner(ρ, xvec, yvec, solver = WignerLaguerre(parallel = true))
    wig4 = wigner(ψ, xvec, yvec, solver = WignerClenshaw())

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 1e-3
    @test sqrt(sum(abs.(wig3 .- wig)) / length(wig)) < 1e-3
    @test sqrt(sum(abs.(wig4 .- wig)) / length(wig)) < 1e-3

    X, Y = meshgrid(xvec, yvec)
    wig_tmp1 = gaussian.(xvec / √2, real(α), 1 / 2)
    wig_tmp2 = gaussian.(yvec / √2, imag(α), 1 / 2)
    wig2 = maximum(wig) * reshape(kron(wig_tmp1, wig_tmp2), 300, 300)

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 0.1
end
