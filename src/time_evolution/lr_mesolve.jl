export lr_mesolve, lr_mesolveProblem, LRTimeEvolutionSol, LRMesolveOptions

#=======================================================#
#                   STRUCT DEFINITIONS
#=======================================================#

struct LRTimeEvolutionSol{TT<:Vector{<:Real},TS<:AbstractVector,TE<:Matrix{ComplexF64},TM<:Vector{<:Integer}}
    times::TT
    z::TS
    B::TS
    expvals::TE
    funvals::TE
    M::TM
end

struct LRMesolveOptions{AlgType<:OrdinaryDiffEqAlgorithm}
    alg::AlgType
    progress::Bool
    err_max::Real
    p0::Real
    atol_inv::Real
    M_max::Integer
    compute_Si::Bool
    is_dynamical::Bool
    adj_condition::String
    Δt::Real
end

function LRMesolveOptions(;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    progress::Bool = true,
    err_max::Real = 0.0,
    p0::Real = 0.0,
    atol_inv::Real = 1e-4,
    M_max::Integer = typemax(Int),
    compute_Si::Bool = true,
    _is_dynamical::Bool = err_max > 0,
    adj_condition::String = "variational",
    Δt::Real = 0.0,
)
    return LRMesolveOptions{typeof(alg)}(
        alg,
        progress,
        err_max,
        p0,
        atol_inv,
        M_max,
        compute_Si,
        _is_dynamical,
        adj_condition,
        Δt,
    )
end

#=======================================================#
#                  ADDITIONAL FUNCTIONS
#=======================================================#

select(x::Real, xarr::AbstractArray, retval = false) = retval ? xarr[argmin(abs.(x .- xarr))] : argmin(abs.(x .- xarr))

@doc raw"""
    _pinv!(A, T1, T2; atol::Real=0.0, rtol::Real=(eps(real(float(oneunit(T))))*min(size(A)...))*iszero(atol)) where T
    Computes the pseudo-inverse of a matrix A, and stores it in T1. If T2 is provided, it is used as a temporary matrix. 
    The algorithm is based on the SVD decomposition of A, and is taken from the Julia package LinearAlgebra.
    The difference with respect to the original function is that the cutoff is done with a smooth function instead of a step function.

    Parameters
    ----------
    A : AbstractMatrix{T}
        The matrix to be inverted.
    T1 : AbstractMatrix{T}
    T2 : AbstractMatrix{T}
        Temporary matrices used in the calculation.
    atol : Real
        Absolute tolerance for the calculation of the pseudo-inverse.   
    rtol : Real
        Relative tolerance for the calculation of the pseudo-inverse.
"""
function _pinv!(
    A::AbstractMatrix{T},
    T1::AbstractMatrix{T},
    T2::AbstractMatrix{T};
    atol::Real = 0.0,
    rtol::Real = (eps(real(float(oneunit(T)))) * min(size(A)...)) * iszero(atol),
) where {T}
    if isdiag(A)
        idxA = diagind(A)
        diagA = view(A, idxA)
        maxabsA = maximum(abs, diagA)
        λ = max(rtol * maxabsA, atol)
        return Matrix(Diagonal(pinv.(diagA) .* 1 ./ (1 .+ (λ ./ real(diagA)) .^ 6)))
    end

    SVD = svd(A)
    λ = max(rtol * maximum(SVD.S), atol)
    SVD.S .= pinv.(SVD.S) .* 1 ./ (1 .+ (λ ./ SVD.S) .^ 6)
    mul!(T2, Diagonal(SVD.S), SVD.U')
    return mul!(T1, SVD.Vt', T2)
end

@doc raw"""
    _calculate_expectation!(p,z,B,idx) where T
    Calculates the expectation values and function values of the operators and functions in p.e_ops and p.f_ops, respectively, and stores them in p.expvals and p.funvals.
    The function is called by the callback _save_affect_lr_mesolve!.

    Parameters
    ----------
    p : NamedTuple
        The parameters of the problem.
    z : AbstractMatrix{T}
        The z matrix.
    B : AbstractMatrix{T}
        The B matrix.
    idx : Integer
        The index of the current time step.
"""
function _calculate_expectation!(p, z, B, idx)
    e_ops = p.e_ops
    f_ops = p.f_ops
    expvals = p.expvals
    funvals = p.funvals
    p.Ml[idx] = p.M

    #Expectation values
    Oz = p.A0
    zOz = p.temp_MM
    @inbounds for (i, e) in enumerate(e_ops)
        mul!(Oz, e, z)
        mul!(zOz, z', Oz)
        expvals[i, idx] = dot(B, zOz)
    end

    #Function values
    @inbounds for (i, f) in enumerate(f_ops)
        funvals[i, idx] = f(p, z, B)
    end
end

#=======================================================#
#                   SAVING FUNCTIONS
#=======================================================#

function _periodicsave_func(integrator)
    ip = integrator.p
    ip.u_save .= integrator.u
    ip.scalars[2] = integrator.t
    return u_modified!(integrator, false)
end

_save_control_lr_mesolve(u, t, integrator) = t in integrator.p.times

function _save_affect_lr_mesolve!(integrator)
    ip = integrator.p
    N, M = ip.N, ip.M
    idx = select(integrator.t, ip.times)

    @views z = reshape(integrator.u[1:N*M], N, M)
    @views B = reshape(integrator.u[N*M+1:end], M, M)
    _calculate_expectation!(ip, z, B, idx)

    if integrator.p.opt.progress
        print("\rProgress: $(round(Int, 100*idx/length(ip.times)))%")
        flush(stdout)
    end
    return u_modified!(integrator, false)
end

#=======================================================#
#                CALLBACK FUNCTIONS
#=======================================================#

@doc raw"""
    _adjM_condition_ratio(u, t, integrator) where T
    Condition for the dynamical rank adjustment based on the ratio between the smallest and largest eigenvalues of the density matrix.
    The spectrum of the density matrix is calculated efficiently using the properties of the SVD decomposition of the matrix.

    Parameters
    ----------
    u : AbstractVector{T}
        The current state of the system.
    t : Real
        The current time.
    integrator : ODEIntegrator
        The integrator of the problem.
"""
function _adjM_condition_ratio(u, t, integrator)
    ip = integrator.p
    opt = ip.opt
    N, M = ip.N, ip.M

    C = ip.A0
    σ = ip.temp_MM
    @views z = reshape(u[1:N*M], N, M)
    @views B = reshape(u[N*M+1:end], M, M)
    mul!(C, z, sqrt(B))
    mul!(σ, C', C)
    p = abs.(eigvals(σ))
    err = p[1] / p[end]

    return (err >= opt.err_max && M < N && M < opt.M_max)
end

@doc raw"""
    _adjM_condition_variational(u, t, integrator) where T
    Condition for the dynamical rank adjustment based on the leakage out of the low-rank manifold.

    Parameters
    ----------
    u : AbstractVector{T}
        The current state of the system.
    t : Real
        The current time.
    integrator : ODEIntegrator
        The integrator of the problem.
"""
function _adjM_condition_variational(u, t, integrator)
    ip = integrator.p
    opt = ip.opt
    N, M = ip.N, ip.M
    err = abs(ip.scalars[1] * sqrt(ip.M))

    return (err >= opt.err_max && M < N && M < opt.M_max)
end

@doc raw"""
    _adjM_affect!(integrator)
    Affect function for the dynamical rank adjustment. It increases the rank of the low-rank manifold by one, and updates the matrices accordingly.
    If Δt>0, it rewinds the integrator to the previous time step.

    Parameters
    ----------
    integrator : ODEIntegrator
        The integrator of the problem.
"""
function _adjM_affect!(integrator)
    ip = integrator.p
    opt = ip.opt
    Δt = opt.Δt
    N, M = ip.N, ip.M

    @views begin
        z = Δt > 0 ? reshape(ip.u_save[1:N*M], N, M) : reshape(integrator.u[1:N*M], N, M)
        B = Δt > 0 ? reshape(ip.u_save[N*M+1:end], M, M) : reshape(integrator.u[N*M+1:end], M, M)
        ψ = ip.L_tilde[:, 1]
        normalize!(ψ)

        z = hcat(z, ψ)
        B = cat(B, opt.p0, dims = (1, 2))
        resize!(integrator, length(z) + length(B))
        integrator.u[1:length(z)] .= z[:]
        integrator.u[length(z)+1:end] .= B[:]
    end

    integrator.p = merge(
        integrator.p,
        (
            M = ip.M + 1,
            L_tilde = similar(z),
            A0 = similar(z),
            Bi = similar(B),
            L = similar(B),
            temp_MM = similar(B),
            S = similar(B),
            Si = similar(B),
        ),
    )
    mul!(integrator.p.S, z', z)
    !(opt.compute_Si) &&
        (integrator.p.Si .= _pinv!(Hermitian(integrator.p.S), integrator.temp_MM, integrator.L, atol = opt.atol_inv))

    if Δt > 0
        integrator.p = merge(integrator.p, (u_save = copy(integrator.u),))
        t0 = ip.scalars[2]
        reinit!(integrator, integrator.u; t0 = t0, erase_sol = false)

        if length(integrator.sol.t) > 1
            idx = findlast(integrator.sol.t .<= t0)
            resize!(integrator.sol.t, idx)
            resize!(integrator.sol.u, idx)
            integrator.saveiter = idx
        end
    end
end

#=======================================================#
#            DYNAMICAL EVOLUTION EQUATIONS
#=======================================================#

@doc raw"""
    dBdz!(du, u, p, t) where T
    Dynamical evolution equations for the low-rank manifold. The function is called by the ODEProblem.

    Parameters
    ----------
    du : AbstractVector{T}
        The derivative of the state of the system.
    u : AbstractVector{T}
        The current state of the system.
    p : NamedTuple
        The parameters of the problem.
    t : Real
        The current time.
"""
function dBdz!(du, u, p, t)
    #NxN
    H, Γ = p.H, p.Γ
    #NxM
    L_tilde, A0 = p.L_tilde, p.A0
    #MxM
    Bi, L, temp_MM, S, Si = p.Bi, p.L, p.temp_MM, p.S, p.Si

    #Get z, B, dz, and dB by reshaping u
    N, M = p.N, p.M
    opt = p.opt

    @views z = reshape(u[1:N*M], N, M)
    @views dz = reshape(du[1:N*M], N, M)
    @views B = reshape(u[N*M+1:end], M, M)
    @views dB = reshape(du[N*M+1:end], M, M)

    #Assign A0 and S
    mul!(S, z', z)
    B .= (B + B') / 2
    mul!(temp_MM, S, B)
    B .= B ./ tr(temp_MM)
    mul!(A0, z, B)

    # Calculate inverse
    opt.compute_Si && (Si .= _pinv!(Hermitian(S), temp_MM, L, atol = opt.atol_inv))
    Bi .= _pinv!(Hermitian(B), temp_MM, L, atol = opt.atol_inv)

    # Calculate the effective Hamiltonian part of L_tilde
    mul!(dz, H, A0)
    mul!(L_tilde, dz, S)
    mul!(L, dz', z)
    mul!(L_tilde, z, L, +1im, -1im)

    # Calculate the jump operators part of L_tilde
    @inbounds for Γi in Γ
        mul!(dz, Γi, z)
        mul!(L, dz', z)
        mul!(temp_MM, B, L)
        mul!(L_tilde, dz, temp_MM, 1, 1)
    end

    # Calculate L
    mul!(L, z', L_tilde)
    mul!(temp_MM, Si, L)
    p.scalars[1] = real(tr(temp_MM) / M)

    # Calculate dz
    mul!(L_tilde, z, temp_MM, -1, 1) #Ltilde is now (Ltilde - z*Si*L)
    mul!(dB, Si, Bi)
    mul!(dz, L_tilde, dB)

    # Calculate dB
    mul!(dB, temp_MM, Si)
    temp_MM .= Si
    lmul!(p.scalars[1], temp_MM)
    return dB .-= temp_MM
end

#=======================================================#
#                   PROBLEM FORMULATION
#=======================================================#

@doc raw"""
    lr_mesolveProblem(H, z, B, tlist, c_ops; e_ops=(), f_ops=(), opt=LRMesolveOptions(), kwargs...) where T
    Formulates the ODEproblem for the low-rank time evolution of the system. The function is called by lr_mesolve.

    Parameters
    ----------
    H : QuantumObject
        The Hamiltonian of the system.
    z : AbstractMatrix{T}
        The initial z matrix.
    B : AbstractMatrix{T}
        The initial B matrix.
    tlist : AbstractVector{T}
        The time steps at which the expectation values and function values are calculated.
    c_ops : AbstractVector{QuantumObject}
        The jump operators of the system.
    e_ops : Tuple{QuantumObject}
        The operators whose expectation values are calculated.
    f_ops : Tuple{Function}
        The functions whose values are calculated.
    opt : LRMesolveOptions
        The options of the problem.
    kwargs : NamedTuple
        Additional keyword arguments for the ODEProblem.
"""
function lr_mesolveProblem(
    H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    z::AbstractArray{T2,2},
    B::AbstractArray{T2,2},
    tlist::AbstractVector,
    c_ops::AbstractVector = [];
    e_ops::Tuple = (),
    f_ops::Tuple = (),
    opt::LRMesolveOptions{AlgType} = LRMesolveOptions(),
    kwargs...,
) where {T1,T2,AlgType<:OrdinaryDiffEqAlgorithm}

    # Formulation of problem
    H -= 0.5im * sum([Γ' * Γ for Γ in c_ops])
    H = get_data(H)
    c_ops = get_data.(c_ops)
    e_ops = get_data.(e_ops)

    t_l = convert(Vector{_FType(H)}, tlist)

    # Initialization of Arrays
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    funvals = Array{ComplexF64}(undef, length(f_ops), length(t_l))
    Ml = Array{Int64}(undef, length(t_l))

    # Initialization of parameters. Scalars represents in order: Tr(S^{-1}L), t0
    p = (
        N = size(z, 1),
        M = size(z, 2),
        H = H,
        Γ = c_ops,
        e_ops = e_ops,
        f_ops = f_ops,
        opt = opt,
        times = t_l,
        expvals = expvals,
        funvals = funvals,
        Ml = Ml,
        L_tilde = similar(z),
        A0 = similar(z),
        Bi = similar(B),
        L = similar(B),
        temp_MM = similar(B),
        S = similar(B),
        Si = similar(B),
        u_save = vcat(vec(z), vec(B)),
        scalars = [0.0, t_l[1]],
    )

    mul!(p.S, z', z)
    p.Si .= pinv(Hermitian(p.S), atol = opt.atol_inv)

    # Initialization of ODEProblem's kwargs
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = [t_l[end]])
    kwargs2 = merge(default_values, kwargs)

    # Initialization of Callbacks
    if !isempty(e_ops) || !isempty(f_ops)
        _calculate_expectation!(p, z, B, 1)
        cb_save = DiscreteCallback(_save_control_lr_mesolve, _save_affect_lr_mesolve!, save_positions = (false, false))
        kwargs2 = merge(
            kwargs2,
            haskey(kwargs2, :callback) ? Dict(:callback => CallbackSet(cb_save, kwargs2[:callback])) :
            Dict(:callback => cb_save),
        )
    end

    if opt.is_dynamical
        if opt.Δt > 0
            cb_periodicsave =
                PeriodicCallback(_periodicsave_func, opt.Δt, final_affect = true, save_positions = (false, false))
            kwargs2 = merge(
                kwargs2,
                haskey(kwargs2, :callback) ? Dict(:callback => CallbackSet(cb_periodicsave, kwargs2[:callback])) :
                Dict(:callback => cb_periodicsave),
            )
        end

        if opt.adj_condition == "variational"
            adj_condition_function = _adjM_condition_variational
        elseif opt.adj_condition == "ratio"
            adj_condition_function = _adjM_condition_ratio
        else
            error("adj_condition must be either 'variational' or 'ratio'")
        end
        cb_adjM = DiscreteCallback(adj_condition_function, _adjM_affect!, save_positions = (false, false))
        kwargs2 = merge(
            kwargs2,
            haskey(kwargs2, :callback) ? Dict(:callback => CallbackSet(cb_adjM, kwargs2[:callback])) :
            Dict(:callback => cb_adjM),
        )
    end

    # Initialization of ODEProblem
    tspan = (t_l[1], t_l[end])
    return ODEProblem(dBdz!, p.u_save, tspan, p; kwargs2...)
end

function lr_mesolve(
    H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    z::AbstractArray{T2,2},
    B::AbstractArray{T2,2},
    tlist::AbstractVector,
    c_ops::AbstractVector = [];
    e_ops::Tuple = (),
    f_ops::Tuple = (),
    opt::LRMesolveOptions{AlgType} = LRMesolveOptions(),
    kwargs...,
) where {T1,T2,AlgType<:OrdinaryDiffEqAlgorithm}
    prob = lr_mesolveProblem(H, z, B, tlist, c_ops; e_ops = e_ops, f_ops = f_ops, opt = opt, kwargs...)
    return lr_mesolve(prob; kwargs...)
end

#=======================================================#
#                  OUTPUT GENNERATION
#=======================================================#

get_z(u::AbstractArray{T}, N::Integer, M::Integer) where {T} = reshape(view(u, 1:M*N), N, M)

get_B(u::AbstractArray{T}, N::Integer, M::Integer) where {T} = reshape(view(u, (M*N+1):length(u)), M, M)

@doc raw"""
    lr_mesolve(prob::ODEProblem; kwargs...)
    Solves the ODEProblem formulated by lr_mesolveProblem. The function is called by lr_mesolve.

    Parameters
    ----------
    prob : ODEProblem
        The ODEProblem formulated by lr_mesolveProblem.
    kwargs : NamedTuple
        Additional keyword arguments for the ODEProblem.
"""
function lr_mesolve(prob::ODEProblem; kwargs...)
    sol = solve(prob, prob.p.opt.alg, tstops = prob.p.times)
    prob.p.opt.progress && print("\n")

    N = prob.p.N
    Ll = length.(sol.u)
    Ml = @. Int((sqrt(N^2 + 4 * Ll) - N) / 2)

    if !haskey(kwargs, :saveat)
        Bt = map(x -> get_B(x[1], N, x[2]), zip(sol.u, Ml))
        zt = map(x -> get_z(x[1], N, x[2]), zip(sol.u, Ml))
    else
        Bt = get_B(sol.u, N, Ml)
        zt = get_z(sol.u, N, Ml)
    end

    return LRTimeEvolutionSol(sol.prob.p.times, zt, Bt, prob.p.expvals, prob.p.funvals, prob.p.Ml)
end
