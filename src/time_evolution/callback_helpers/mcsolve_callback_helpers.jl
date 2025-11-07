#=
Helper functions for the mcsolve callbacks.
=#

struct SaveFuncMCSolve{TE,IT,TEXPV} <: AbstractSaveFunc
    e_ops::TE
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncMCSolve)(u, t, integrator) = _save_func_mcsolve(u, integrator, f.e_ops, f.iter, f.expvals)

_get_save_callback_idx(cb, ::Type{SaveFuncMCSolve}) = _mcsolve_has_continuous_jump(cb) ? 1 : 2

##
struct LindbladJump{
    T1,
    T2,
    RNGType<:AbstractRNG,
    RandT,
    CT<:AbstractVector,
    WT<:AbstractVector,
    JTT<:AbstractVector,
    JWT<:AbstractVector,
    JTWIT,
}
    c_ops::T1
    c_ops_herm::T2
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
    f.c_ops_herm,
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
    cache_mc = _mc_get_jump_callback(integrator).affect!.cache_mc

    copyto!(cache_mc, u)
    normalize!(cache_mc)
    ψ = cache_mc
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    u_modified!(integrator, false)
    return nothing
end

function _generate_mcsolve_kwargs(ψ0, T, e_ops, tlist, c_ops, jump_callback, rng, kwargs)
    c_ops_data = get_data.(c_ops)
    c_ops_herm_data = map(op -> op' * op, c_ops_data)

    cache_mc = similar(ψ0.data, T)
    weights_mc = Vector{Float64}(undef, length(c_ops))
    cumsum_weights_mc = similar(weights_mc)

    col_times = Vector{Float64}(undef, COL_TIMES_WHICH_INIT_SIZE)
    col_which = Vector{Int}(undef, COL_TIMES_WHICH_INIT_SIZE)
    col_times_which_idx = Ref(1)

    random_n = Ref(rand(rng))

    _affect! = LindbladJump(
        c_ops_data,
        c_ops_herm_data,
        rng,
        random_n,
        cache_mc,
        weights_mc,
        cumsum_weights_mc,
        col_times,
        col_which,
        col_times_which_idx,
    )

    if jump_callback isa DiscreteLindbladJumpCallback
        cb1 = DiscreteCallback(_mcsolve_discrete_condition, _affect!, save_positions = (false, false))
    else
        cb1 = ContinuousCallback(
            _mcsolve_continuous_condition,
            _affect!,
            nothing,
            interp_points = jump_callback.interp_points,
            save_positions = (false, false),
        )
    end

    kwargs2 = _kwargs_set_tstops(kwargs, tlist)
    if e_ops isa Nothing
        # We are implicitly saying that we don't have a `Progress`
        kwargs3 =
            haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, kwargs2.callback),)) :
            merge(kwargs2, (callback = cb1,))
        return kwargs3
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(tlist))

        _save_func = SaveFuncMCSolve(get_data.(e_ops), Ref(1), expvals)
        cb2 = FunctionCallingCallback(_save_func, funcat = tlist)
        kwargs3 =
            haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, cb2, kwargs2.callback),)) :
            merge(kwargs2, (callback = CallbackSet(cb1, cb2),))
        return kwargs3
    end
end

function _lindblad_jump_affect!(
    integrator,
    c_ops,
    c_ops_herm,
    traj_rng,
    random_n,
    cache_mc,
    weights_mc,
    cumsum_weights_mc,
    col_times,
    col_which,
    col_times_which_idx,
)
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        weights_mc[i] = real(dot(ψ, c_ops_herm[i], ψ))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    r = rand(traj_rng) * sum(weights_mc)
    collapse_idx = getindex(1:length(weights_mc), findfirst(>(r), cumsum_weights_mc))
    mul!(cache_mc, c_ops[collapse_idx], ψ)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)

    random_n[] = rand(traj_rng)

    idx = col_times_which_idx[]
    @inbounds col_times[idx] = integrator.t
    @inbounds col_which[idx] = collapse_idx
    col_times_which_idx[] += 1
    if col_times_which_idx[] > length(col_times)
        resize!(col_times, length(col_times) + COL_TIMES_WHICH_INIT_SIZE)
        resize!(col_which, length(col_which) + COL_TIMES_WHICH_INIT_SIZE)
    end
    u_modified!(integrator, true)
    return nothing
end

_mcsolve_continuous_condition(u, t, integrator) =
    @inbounds _mc_get_jump_callback(integrator).affect!.random_n[] - real(dot(u, u))

_mcsolve_discrete_condition(u, t, integrator) =
    @inbounds real(dot(u, u)) < _mc_get_jump_callback(integrator).affect!.random_n[]

##

function _mc_get_jump_callback(sol::AbstractODESolution)
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    return _mc_get_jump_callback(kwargs.callback) # There is always the Jump callback
end
_mc_get_jump_callback(integrator::AbstractODEIntegrator) = _mc_get_jump_callback(integrator.opts.callback)
_mc_get_jump_callback(cb::CallbackSet) =
    if _mcsolve_has_continuous_jump(cb)
        return cb.continuous_callbacks[1]
    else
        return cb.discrete_callbacks[1]
    end
_mc_get_jump_callback(cb::ContinuousCallback) = cb
_mc_get_jump_callback(cb::DiscreteCallback) = cb

##

#=
    With this function we extract the c_ops and c_ops_herm from the LindbladJump `affect!` function of the callback of the integrator.
    This callback can be a DiscreteLindbladJumpCallback or a ContinuousLindbladJumpCallback.
=#
function _mcsolve_get_c_ops(integrator::AbstractODEIntegrator)
    cb = _mc_get_jump_callback(integrator)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.c_ops, cb.affect!.c_ops_herm
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

    if _mcsolve_has_continuous_jump(cb)
        idx = 1
        if cb_discrete[idx].affect!.func isa SaveFuncMCSolve
            e_ops = cb_discrete[idx].affect!.func.e_ops
            expvals = similar(cb_discrete[idx].affect!.func.expvals)
            _save_func = SaveFuncMCSolve(e_ops, Ref(1), expvals)
            cb_save = (FunctionCallingCallback(_save_func, funcat = tlist),)
        else
            cb_save = ()
        end

        _jump_affect! = _similar_affect!(cb_continuous[1].affect!, traj_rng)
        cb_jump = _modify_field(cb_continuous[1], :affect!, _jump_affect!)

        return CallbackSet((cb_jump, cb_continuous[2:end]...), (cb_save..., cb_discrete[2:end]...))
    else
        idx = 2
        if cb_discrete[idx].affect!.func isa SaveFuncMCSolve
            e_ops = cb_discrete[idx].affect!.func.e_ops
            expvals = similar(cb_discrete[idx].affect!.func.expvals)
            _save_func = SaveFuncMCSolve(e_ops, Ref(1), expvals)
            cb_save = (FunctionCallingCallback(_save_func, funcat = tlist),)
        else
            cb_save = ()
        end

        _jump_affect! = _similar_affect!(cb_discrete[1].affect!, traj_rng)
        cb_jump = _modify_field(cb_discrete[1], :affect!, _jump_affect!)

        return CallbackSet(cb_continuous, (cb_jump, cb_save..., cb_discrete[3:end]...))
    end
end
function _mcsolve_initialize_callbacks(cb::CBT, tlist, traj_rng) where {CBT<:Union{ContinuousCallback,DiscreteCallback}}
    _jump_affect! = _similar_affect!(cb.affect!, traj_rng)
    return _modify_field(cb, :affect!, _jump_affect!)
end

#=
    _similar_affect!

Return a new LindbladJump with the same fields as the input LindbladJump but with new memory.
=#
function _similar_affect!(affect::LindbladJump, traj_rng)
    random_n = Ref(rand(traj_rng))
    cache_mc = similar(affect.cache_mc)
    weights_mc = similar(affect.weights_mc)
    cumsum_weights_mc = similar(affect.cumsum_weights_mc)
    col_times = similar(affect.col_times)
    col_which = similar(affect.col_which)
    col_times_which_idx = Ref(1)

    return LindbladJump(
        affect.c_ops,
        affect.c_ops_herm,
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

_mcsolve_has_continuous_jump(cb::CallbackSet) =
    (length(cb.continuous_callbacks) > 0) && (cb.continuous_callbacks[1].affect! isa LindbladJump)
_mcsolve_has_continuous_jump(cb::ContinuousCallback) = true
_mcsolve_has_continuous_jump(cb::DiscreteCallback) = false
