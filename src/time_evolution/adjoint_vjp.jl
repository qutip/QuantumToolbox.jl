#=
Analytic vector-Jacobian products (VJPs) of the linear ODEs solved by `sesolve` and `mesolve`:

    du/dt = L(p, t) u,    L(p, t) = Σₖ cₖ(p, t) Aₖ,

with constant matrices Aₖ and coefficients cₖ depending on (p, t) only. For a cotangent λ,

    ∂ℓ/∂u = Σₖ conj(cₖ) Aₖᴴ λ,    ∂ℓ/∂p = Re Σₖ ⟨λ, Aₖ u⟩ ∇ₚcₖ.

They are given to SciMLSensitivity.jl as the `vjp` and `vjp_p` of the `ODEFunction`, which all its
continuous adjoints use instead of differentiating the right-hand side with the `autojacvec`
backend. The matrix-vector products are then never differentiated, and only the scalar
coefficients are (with ForwardDiff.jl).
=#

# Split L into a tuple of terms (cₖ, Aₖ), or return `nothing` if it is not of this form.
_linear_terms(L::MatrixOperator) = isconstant(L) ? ((true, L.A),) : nothing
_linear_terms(::IdentityOperator) = ((true, I),)
function _linear_terms(L::ScaledOperator)
    terms = _linear_terms(L.L)
    return isnothing(terms) ? nothing : map(((c, A),) -> (c isa Bool ? L.λ : L.λ * c, A), terms)
end
function _linear_terms(L::AddedOperator)
    terms = map(_linear_terms, L.ops)
    return any(isnothing, terms) ? nothing : reduce((a, b) -> (a..., b...), terms)
end
_linear_terms(::Any) = nothing

# cₖ(p, t), computed without updating the operator, so that `p` can also hold dual numbers.
_coef_value(c::Number, u, p, t) = c
_coef_value(c::ScalarOperator, u, p, t) = c.update_func(c.val, u, p, t)
_coef_value(c::SciMLOperators.AddedScalarOperator, u, p, t) = sum(map(op -> _coef_value(op, u, p, t), c.ops))
_coef_value(c::SciMLOperators.ComposedScalarOperator, u, p, t) = prod(map(op -> _coef_value(op, u, p, t), c.ops))
_coef_value(c::SciMLOperators.InvertedScalarOperator, u, p, t) = inv(_coef_value(c.λ, u, p, t))
_coef_value(c::SciMLOperators.AbstractSciMLScalarOperator, u, p, t) =
    SciMLOperators.concretize(SciMLOperators.update_coefficients(c, u, p, t))

# The state VJP: dλ = Σₖ conj(cₖ) Aₖᴴ λ
struct LinearOperatorVJP{TT <: Tuple}
    terms::TT
end

function (f::LinearOperatorVJP)(dλ, λ, u, p, t)
    fill!(dλ, zero(eltype(dλ)))
    foreach(((c, A),) -> mul!(dλ, A', λ, conj(_coef_value(c, u, p, t)), true), f.terms)
    return dλ
end

# The parameter VJP: dp = ∇ₚ Re Σₖ cₖ(p, t) ⟨λ, Aₖ u⟩, only for the terms depending on (p, t)
struct LinearOperatorParamVJP{TT <: Tuple}
    terms::TT
end

function (f::LinearOperatorParamVJP)(dp, λ, u, p, t)
    overlaps = map(((c, A),) -> dot(λ, A, u), f.terms)
    g = x -> real(sum(map(((c, _), w) -> _coef_value(c, u, x, t) * w, f.terms, overlaps)))
    return ForwardDiff.gradient!(dp, g, p)
end

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
function _time_evolution_ode_function(U::AbstractSciMLOperator, params, ::Val{inplace}) where {inplace}
    terms = inplace ? _linear_terms(U) : nothing
    isnothing(terms) && return ODEFunction{inplace, FullSpecialize}(U)

    terms_p = filter(((c, _),) -> !isconstant(c), terms)
    vjp_p = if isempty(terms_p)
        _zero_param_vjp
    elseif params isa AbstractVector{<:Real}
        LinearOperatorParamVJP(terms_p)
    else
        nothing
    end
    return ODEFunction{inplace, FullSpecialize}(U; vjp = LinearOperatorVJP(terms), vjp_p = vjp_p)
end
