#=
Helper functions for the mcsolve callbacks.
=#

struct SaveFuncMCSolve{TE, IT, TEXPV} <: AbstractSaveFunc
    e_ops::TE
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncMCSolve)(u, t, integrator) = _save_func_mcsolve(u, integrator, f.e_ops, f.iter, f.expvals)

##
struct LindbladJump{
        T1,
        RNGType <: AbstractRNG,
        RandT,
        CT <: AbstractVector,
        WT <: AbstractVector,
        JTT <: AbstractVector,
        JWT <: AbstractVector,
        JTWIT,
    }
    c_ops::T1
    traj_rng::RNGType
    random_n::RandT
    cache_mc::CT
    weights_mc::WT
    cumsum_weights_mc::WT
    col_times::JTT
    col_which::JWT
    col_times_which_idx::JTWIT
end

(f::LindbladJump)(integrator) = _lindblad_jump_affect!(
    integrator,
    f.c_ops,
    f.traj_rng,
    f.random_n,
    f.cache_mc,
    f.weights_mc,
    f.cumsum_weights_mc,
    f.col_times,
    f.col_which,
    f.col_times_which_idx,
)

##

function _save_func_mcsolve(u, integrator, e_ops, iter, expvals)
    _mcsolve_expect!(view(expvals, :, iter[]), e_ops, u, integrator)
    iter[] += 1

    derivative_discontinuity!(integrator, false)
    return nothing
end

# Write explicit functions for cleaner code and allow other packages to use them through multiple dispatch.
function _mcsolve_expect!(expvals, e_ops, u, integrator)
    norm2 = real(dot(u, u))
    expect_f = op -> dot(u, op, u) / norm2

    # @. expvals = expect_f(e_ops) allocates memory when e_ops is a tuple
    if e_ops isa Tuple
        expvals .= map(expect_f, e_ops)
    else
        @. expvals = expect_f(e_ops)
    end
    return expvals
end

function _mcsolve_jump_weights!(weights, c_ops, cache_mc, integrator, ψ = integrator.u, t = integrator.t)
    p = integrator.p
    @inbounds for i in eachindex(weights)
        c_ops[i](cache_mc, ψ, nothing, p, t)
        weights[i] = real(dot(cache_mc, cache_mc))
    end
    return weights
end

function _mcsolve_jump!(integrator, c_ops, i, cache_mc)
    c_ops[i](cache_mc, integrator.u, nothing, integrator.p, integrator.t)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)
    return nothing
end

function _generate_mcsolve_kwargs(ψ0, T, e_ops, tlist, c_ops, rng, kwargs; jump_derivative = Val(false), jump_log = Val(false))
    cache_mc = similar(ψ0.data, T)

    c_ops_data = map(op -> get_data(cache_operator(QobjEvo(op), cache_mc)), c_ops)

    weights_mc = Vector{Float64}(undef, length(c_ops))
    cumsum_weights_mc = similar(weights_mc)

    col_times = Vector{Float64}(undef, COL_TIMES_WHICH_INIT_SIZE)
    col_which = Vector{Int}(undef, COL_TIMES_WHICH_INIT_SIZE)
    col_times_which_idx = Ref(1)

    random_n = Ref(zero(_float_type(T)))

    _affect! = LindbladJump(
        c_ops_data,
        rng,
        random_n,
        cache_mc,
        weights_mc,
        cumsum_weights_mc,
        col_times,
        col_which,
        col_times_which_idx,
    )

    # `interp_points = 0` is exact: the norm of the state decreases monotonically between jumps, so it crosses the random threshold
    # at most once per step, and the crossing is visible from the sign of the condition at the two step endpoints.
    cb1 = ContinuousCallback(
        _mcsolve_jump_condition(jump_derivative, jump_log),
        _affect!,
        nothing,
        initialize = _lindblad_jump_initialize!,
        interp_points = 0,
        save_positions = (false, false),
    )

    if e_ops isa Nothing
        # We are implicitly saying that we don't have a `Progress`
        kwargs2 = _merge_kwargs_with_callback(kwargs, cb1)
    else
        expvals = Array{_complex_float_type(T)}(undef, length(e_ops), length(tlist))

        _save_func = SaveFuncMCSolve(get_data.(e_ops), Ref(1), expvals)
        cb2 = FunctionCallingCallback(_save_func, funcat = tlist)
        kwargs2 = _merge_kwargs_with_callback(kwargs, CallbackSet(cb1, cb2))
    end
    return kwargs2
end

function _lindblad_jump_initialize!(cb, u, t, integrator)
    affect! = cb.affect!
    affect!.random_n[] = rand(affect!.traj_rng)
    return nothing
end

function _lindblad_jump_affect!(
        integrator,
        c_ops,
        traj_rng,
        random_n,
        cache_mc,
        weights_mc,
        cumsum_weights_mc,
        col_times,
        col_which,
        col_times_which_idx,
    )
    _mcsolve_jump_weights!(weights_mc, c_ops, cache_mc, integrator)
    cumsum!(cumsum_weights_mc, weights_mc)
    r = rand(traj_rng) * last(cumsum_weights_mc)
    collapse_idx = something(findfirst(>(r), cumsum_weights_mc), lastindex(cumsum_weights_mc))
    _mcsolve_jump!(integrator, c_ops, collapse_idx, cache_mc)

    random_n[] = rand(traj_rng)

    idx = col_times_which_idx[]
    @inbounds col_times[idx] = integrator.t
    @inbounds col_which[idx] = collapse_idx
    col_times_which_idx[] += 1
    if col_times_which_idx[] > length(col_times)
        resize!(col_times, length(col_times) + COL_TIMES_WHICH_INIT_SIZE)
        resize!(col_which, length(col_which) + COL_TIMES_WHICH_INIT_SIZE)
    end
    derivative_discontinuity!(integrator, true)
    return nothing
end

function _mcsolve_continuous_condition(u, t, integrator, ::Val{Log} = Val(false)) where {Log}
    r = _mc_get_jump_callback(integrator).affect!.random_n[]
    s = real(dot(u, u))
    iszero(r) && return -one(s)   # a zero threshold has no finite crossing
    return Log ? log(r) - log(s) : r - s
end

function _mcsolve_continuous_derivative(u, t, integrator, ::Val{Log}) where {Log}
    s = real(dot(u, u))
    iszero(s) && return zero(s)   # bisect at an underflowed endpoint
    jump = _mc_get_jump_callback(integrator).affect!
    _mcsolve_jump_weights!(jump.weights_mc, jump.c_ops, jump.cache_mc, integrator, u, t)
    rate = sum(jump.weights_mc)
    return Log ? rate / s : rate
end

function _mcsolve_jump_condition(::Val{Derivative}, logarithm::Val{Log}) where {Derivative, Log}
    condition = Log ? (u, t, integrator) -> _mcsolve_continuous_condition(u, t, integrator, logarithm) :
        _mcsolve_continuous_condition
    derivative(u, t, integrator) = _mcsolve_continuous_derivative(u, t, integrator, logarithm)
    return Derivative ? ConditionWithDerivative(condition, derivative) : condition
end

##

function _mc_get_jump_callback(sol::AbstractODESolution)
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    return _mc_get_jump_callback(kwargs.callback) # There is always the Jump callback
end
_mc_get_jump_callback(integrator::AbstractODEIntegrator) = _mc_get_jump_callback(integrator.opts.callback)
_mc_get_jump_callback(cb::CallbackSet) = cb.continuous_callbacks[1] # The jump callback is always the first continuous callback
_mc_get_jump_callback(cb::ContinuousCallback) = cb

##

#=
    With this function we extract the c_ops from the LindbladJump `affect!` function of the callback of the integrator.
=#
function _mcsolve_get_c_ops(integrator::AbstractODEIntegrator)
    cb = _mc_get_jump_callback(integrator)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.c_ops
    end
end

#=
    _mcsolve_initialize_callbacks(prob, tlist)

Return the same callbacks of the `prob`, but with the `iter` variable reinitialized to 1 and the `expvals` variable reinitialized to a new matrix.
=#
function _mcsolve_initialize_callbacks(prob, tlist, traj_rng)
    cb = prob.kwargs[:callback]
    return _mcsolve_initialize_callbacks(cb, tlist, traj_rng)
end
function _mcsolve_initialize_callbacks(cb::CallbackSet, tlist, traj_rng)
    cb_continuous = cb.continuous_callbacks
    cb_discrete = cb.discrete_callbacks

    if cb_discrete[1].affect!.func isa SaveFuncMCSolve
        e_ops = cb_discrete[1].affect!.func.e_ops
        expvals = similar(cb_discrete[1].affect!.func.expvals)
        _save_func = SaveFuncMCSolve(e_ops, Ref(1), expvals)
        cb_save = (FunctionCallingCallback(_save_func, funcat = tlist),)
    else
        cb_save = ()
    end

    _jump_affect! = _similar_affect!(cb_continuous[1].affect!, traj_rng)
    cb_jump = _modify_field(cb_continuous[1], :affect!, _jump_affect!)

    return CallbackSet((cb_jump, cb_continuous[2:end]...), (cb_save..., cb_discrete[2:end]...))
end
function _mcsolve_initialize_callbacks(cb::ContinuousCallback, tlist, traj_rng)
    _jump_affect! = _similar_affect!(cb.affect!, traj_rng)
    return _modify_field(cb, :affect!, _jump_affect!)
end

#=
    _similar_affect!

Return a new LindbladJump with the same fields as the input LindbladJump but with new memory.
=#
function _similar_affect!(affect::LindbladJump, traj_rng)
    random_n = Ref(zero(eltype(affect.random_n))) # drawn from `traj_rng` by `_lindblad_jump_initialize!`
    cache_mc = similar(affect.cache_mc)
    weights_mc = similar(affect.weights_mc)
    cumsum_weights_mc = similar(affect.cumsum_weights_mc)
    col_times = similar(affect.col_times)
    col_which = similar(affect.col_which)
    col_times_which_idx = Ref(1)

    c_ops = map(op -> isconstant(op) ? op : deepcopy(op), affect.c_ops)

    return LindbladJump(
        c_ops,
        traj_rng,
        random_n,
        cache_mc,
        weights_mc,
        cumsum_weights_mc,
        col_times,
        col_which,
        col_times_which_idx,
    )
end

Base.@constprop :aggressive function _modify_field(obj::T, field_name::Symbol, field_val) where {T}
    # Create a NamedTuple of fields, deepcopying only the selected ones
    fields = (name != field_name ? (getfield(obj, name)) : field_val for name in fieldnames(T))
    # Reconstruct the struct with the updated fields
    return Base.typename(T).wrapper(fields...)
end
