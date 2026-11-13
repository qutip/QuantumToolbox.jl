#=
Analytic vector-Jacobian products (VJPs) of the linear ODEs solved by `sesolve` and `mesolve`:

    du/dt = L(p, t)[u],    L(p, t)[u] = Σₖ cₖ(p, t) Aₖ u Bₖ,

with constant matrices Aₖ and Bₖ (Bₖ = I for vectorized states, while the matrix form of `mesolve`
also multiplies the density matrix from the right) and coefficients cₖ depending on (p, t) only.
For a cotangent λ,

    ∂ℓ/∂u = Σₖ conj(cₖ) Aₖᴴ λ Bₖᴴ,    ∂ℓ/∂p = Re Σₖ ⟨λ, Aₖ u Bₖ⟩ ∇ₚcₖ.

They are given to SciMLSensitivity.jl as the `vjp` and `vjp_p` of the `ODEFunction`, which all its
continuous adjoints use instead of differentiating the right-hand side with the `autojacvec`
backend. The matrix products are then never differentiated, and only the scalar coefficients are
(with ForwardDiff.jl).
=#

# Split L into a tuple of terms (cₖ, Aₖ, Bₖ), or return `nothing` if it is not of this form.
_linear_terms(A::AbstractMatrix) = ((true, A, I),)
_linear_terms(L::MatrixOperator) = isconstant(L) ? ((true, L.A, I),) : nothing
_linear_terms(::IdentityOperator) = ((true, I, I),)
_linear_terms(L::ScaledOperator) = _map_terms(((c, A, B),) -> (_mul_coef(L.λ, c), A, B), _linear_terms(L.L))
_linear_terms(L::AddedOperator) = _concat_terms(map(_linear_terms, L.ops))
_linear_terms(L::ComposedOperator) = reduce(_compose_terms, map(_linear_terms, L.ops))
_linear_terms(L::SpostSuperOperator) = _map_terms(((c, B, _),) -> (c, I, B), _left_terms(L.R))
function _linear_terms(L::SprePostSuperOperator)
    left, right = _left_terms(L.L), _left_terms(L.R)
    (isnothing(left) || isnothing(right)) && return nothing
    return _concat_terms(map(((c, A, _),) -> map(((d, B, _),) -> (_mul_coef(c, d), A, B), right), left))
end
_linear_terms(::Any) = nothing

_map_terms(_, ::Nothing) = nothing
_map_terms(f, terms) = (mapped = map(f, terms); any(isnothing, mapped) ? nothing : mapped)
_concat_terms(terms) = any(isnothing, terms) ? nothing : reduce((a, b) -> (a..., b...), terms)

# The terms of an operator multiplying only from the left (the factors of `spost` and `sprepost`)
_left_terms(L) = _map_terms(((c, A, B),) -> B isa UniformScaling ? (c, A, B) : nothing, _linear_terms(L))

# (L₁ ∘ L₂)[u] = c₁ c₂ A₁ A₂ u B₂ B₁
_compose_terms(::Nothing, _) = nothing
_compose_terms(_, ::Nothing) = nothing
_compose_terms(t₁, t₂) = _concat_terms(
    map(((c₁, A₁, B₁),) -> map(((c₂, A₂, B₂),) -> (_mul_coef(c₁, c₂), _mul_matrix(A₁, A₂), _mul_matrix(B₂, B₁)), t₂), t₁),
)
_mul_coef(c, ::Bool) = c
_mul_coef(::Bool, d) = d
_mul_coef(::Bool, ::Bool) = true
_mul_coef(c, d) = c * d
_mul_matrix(A, B) = A * B
_mul_matrix(A, ::UniformScaling) = A
_mul_matrix(::UniformScaling, B) = B
_mul_matrix(::UniformScaling, ::UniformScaling) = I

# cₖ(p, t), computed without updating the operator, so that `p` can also hold dual numbers.
_coef_value(c::Number, u, p, t) = c
_coef_value(c::ScalarOperator, u, p, t) = c.update_func(c.val, u, p, t)
_coef_value(c::SciMLOperators.AddedScalarOperator, u, p, t) = sum(map(op -> _coef_value(op, u, p, t), c.ops))
_coef_value(c::SciMLOperators.ComposedScalarOperator, u, p, t) = prod(map(op -> _coef_value(op, u, p, t), c.ops))
_coef_value(c::SciMLOperators.InvertedScalarOperator, u, p, t) = inv(_coef_value(c.λ, u, p, t))
_coef_value(c::SciMLOperators.AbstractSciMLScalarOperator, u, p, t) =
    SciMLOperators.concretize(SciMLOperators.update_coefficients(c, u, p, t))

# The adjoint methods pass λ and dλ as vectors, also for matrix states (propagators and matrix form)
_reshape_as(λ, u) = size(λ) == size(u) ? λ : reshape(λ, size(u))

# The state VJP: dλ = Σₖ conj(cₖ) Aₖᴴ λ Bₖᴴ
struct LinearOperatorVJP{TT <: Tuple, CT}
    terms::TT
    cache::CT # for the terms with both Aₖ and Bₖ
end

function (f::LinearOperatorVJP)(dλ, λ, u, p, t)
    dΛ, Λ = _reshape_as(dλ, u), _reshape_as(λ, u)
    fill!(dΛ, zero(eltype(dΛ)))
    foreach(((c, A, B),) -> _add_adjoint_term!(dΛ, A, Λ, B, conj(_coef_value(c, u, p, t)), f.cache), f.terms)
    return dλ
end

_add_adjoint_term!(dΛ, A, Λ, ::UniformScaling, α, _) = mul!(dΛ, A', Λ, α, true)
_add_adjoint_term!(dΛ, ::UniformScaling, Λ, B, α, _) = mul!(dΛ, Λ, B', α, true)
_add_adjoint_term!(dΛ, ::UniformScaling, Λ, ::UniformScaling, α, _) = axpy!(α, Λ, dΛ)
_add_adjoint_term!(dΛ, A, Λ, B, α, cache) = mul!(dΛ, mul!(cache, A', Λ), B', α, true)

# The parameter VJP: dp = ∇ₚ Re Σₖ cₖ(p, t) ⟨λ, Aₖ u Bₖ⟩, only for the terms depending on (p, t)
struct LinearOperatorParamVJP{TT <: Tuple, CT}
    terms::TT
    cache::CT # for the overlaps of matrix states
end

function (f::LinearOperatorParamVJP)(dp, λ, u, p, t)
    Λ = _reshape_as(λ, u)
    overlaps = map(((c, A, B),) -> _overlap(Λ, A, u, B, f.cache), f.terms)
    g = x -> real(sum(map(((c, _, _), w) -> _coef_value(c, u, x, t) * w, f.terms, overlaps)))
    return ForwardDiff.gradient!(dp, g, p)
end

_overlap(Λ::AbstractVector, A, u::AbstractVector, ::UniformScaling, _) = dot(Λ, A, u)
_overlap(Λ, A, u, ::UniformScaling, cache) = dot(Λ, mul!(cache[1], A, u))
_overlap(Λ, ::UniformScaling, u, B, cache) = dot(Λ, mul!(cache[1], u, B))
_overlap(Λ, ::UniformScaling, u, ::UniformScaling, _) = dot(Λ, u)
_overlap(Λ, A, u, B, cache) = dot(mul!(cache[1], A', Λ), mul!(cache[2], u, B))

# The operator does not depend on the parameters
_zero_param_vjp(dp, λ, u, p, t) = fill!(dp, zero(eltype(dp)))

# They are all in-place. Saying it avoids the method reflection of the generic `isinplace`,
# which the `ODEFunction` constructor would run at every problem construction.
SciMLBase.isinplace(
    ::Union{LinearOperatorVJP, LinearOperatorParamVJP, typeof(_zero_param_vjp)},
    inplace_param_number, fname = "f", iip_preferred = true; kwargs...,
) = true

# The `ODEFunction` of `sesolve` and `mesolve`, with the analytic VJPs when the operator has the
# form above. Otherwise (or for out-of-place problems, or parameters that are not a vector), the
# adjoint differentiates the right-hand side as usual.
function _time_evolution_ode_function(U::AbstractSciMLOperator, u0, params, ::Val{inplace}) where {inplace}
    terms = inplace ? _linear_terms(U) : nothing
    isnothing(terms) && return ODEFunction{inplace, FullSpecialize}(U)

    two_sided = any(((c, A, B),) -> !(A isa UniformScaling || B isa UniformScaling), terms)
    vjp = LinearOperatorVJP(terms, two_sided ? similar(u0) : nothing)

    terms_p = filter(((c, _, _),) -> !isconstant(c), terms)
    vjp_p = if isempty(terms_p)
        _zero_param_vjp
    elseif params isa AbstractVector{<:Real}
        LinearOperatorParamVJP(terms_p, u0 isa AbstractVector ? nothing : (similar(u0), similar(u0)))
    else
        nothing
    end
    return ODEFunction{inplace, FullSpecialize}(U; vjp = vjp, vjp_p = vjp_p)
end
