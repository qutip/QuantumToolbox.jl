export ssesolveProblem

#TODO: Check if works in GPU
function _ssesolve_update_coefficients!(ψ, coefficients, c_ops)
    _get_en = op -> real(dot(ψ, op, ψ)) #this is en/2: <Sn + Sn'>/2 = Re<Sn>
    @. coefficients[2:end-1] = _get_en(c_ops) #coefficients of the OperatorSum: Σ Sn * en/2
    coefficients[end] = - sum(x->x^2, coefficients[2:end-1]) / 2 #this last coefficient is -Σen^2/8
    return nothing
end

function ssesolve_ti_drift!(du, u, p, t)
    _ssesolve_update_coefficients!(u, p.K.coefficients, p.c_ops)
    
    mul!(du, p.K, u)
    return nothing
end

function ssesolve_ti_diffusion!(du, u, p, t)
    D = p.D
    @views en = p.K.coefficients[2:end-1]

    # du:(H,W). du_reshaped:(H*W,). 
    # H:Hilbert space dimension, W: number of c_ops
    du_reshaped = reshape(du, :)
    mul!(du_reshaped, D, u) #du[:,i] = D[i] * u

    du .-= u .* reshape(en, 1, :) #du[:,i] -= en[i] * u
    return nothing
end


function ssesolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix,TJC<:LindbladJumpCallbackType}

    H.dims != ψ0.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    is_time_dependent = !(H_t isa Nothing)
    progress_bar_val = makeVal(progress_bar)

    t_l = convert(Vector{Float64}, tlist) # Convert it into Float64 to avoid type instabilities for OrdinaryDiffEq.jl

    ϕ0 = get_data(ψ0)

    H_eff = get_data(H - T2(0.5im) * mapreduce(op -> op' * op, +, c_ops))
    c_ops2 = get_data.(c_ops)

    coefficients = [1.0, fill(0.0, length(c_ops)+1)...]
    operators = [-1im * H_eff, c_ops2..., I(prod(H.dims))]
    K = OperatorSum(coefficients, operators)
    _ssesolve_update_coefficients!(ϕ0, K.coefficients, c_ops2)

    D = vcat(c_ops2...)

    progr = ProgressBar(length(t_l), enable = getVal(progress_bar_val))

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(t_l))
        e_ops2 = MT1[]
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        e_ops2 = get_data.(e_ops)
        is_empty_e_ops = isempty(e_ops)
    end

    p = (
        U = -1im * get_data(H),
        K = K,
        D = D,
        e_ops = e_ops2,
        c_ops = c_ops2,
        expvals = expvals,
        progr = progr,
        Hdims = H.dims,
        H_t = H_t,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = e_ops isa Nothing ? t_l : [t_l[end]]
    default_values = (saveat = saveat, )
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, progress_bar_val, t_l, kwargs2)

    tspan = (t_l[1], t_l[end])
    A = similar(ϕ0, length(ϕ0), length(c_ops))
    return SDEProblem(ssesolve_ti_drift!, ssesolve_ti_diffusion!, ϕ0, tspan, p; noise_rate_prototype=A, kwargs3...)
end
