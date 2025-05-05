# ======== vectorial forces ==========
struct Ψ1{f_a}
    n_a::Int
    uvec::Symbol
end
Ψ1(f_a, n::Int, uvec::Symbol) = Ψ1{f_a}(n, uvec)
@inline term_mag(trm::Ψ1{f_a}, a) where f_a = f_a(trm.n_a, a)
@inline term_mag(trm::Ψ1{f_a}, a, a0) where f_a = f_a(trm.n_a, a, a0)


# force terms
struct Ψ2{f_a, f_b, f_ab}
    n_a::Int
    n_b::Int
    n_ab::Int
    uvec::Symbol
end
Ψ2(f_a, f_b, f_ab, n_a::Int, n_b::Int, n_ab::Int, uvec::Symbol) = Ψ2{f_a, f_b, f_ab}(n_a, n_b, n_ab, uvec)


function term_mag(trm::Ψ2{f_a, f_b, f_ab}, a, b, ab, a0, b0) where {f_a, f_b, f_ab}
    fmag = 1.0
    fmag *= f_a(trm.n_a, a, a0)
    fmag *= f_b(trm.n_b, b, b0)
    fmag *= f_ab(trm.n_ab, ab)
    return fmag
end


function term_mag(trm::Ψ2{f_a, f_b, f_ab}, a, b, ab, a0b0) where {f_a, f_b, f_ab}
    fmag = 1.0
    fmag *= f_a(trm.n_a, a, a0b0[1])
    fmag *= f_b(trm.n_b, b, a0b0[2])
    fmag *= f_ab(trm.n_ab, ab)
    return fmag
end

# ===================================



# ============ geometry =============

abstract type Geometry end
struct FlatND <: Geometry end
struct PbcND <: Geometry
    lbox::Vector{Float64}
    ispbc::BitVector
end
struct SphericalSurface <: Geometry
    R::Float64
end
# ===================================



# =========== Basis functions ==========

abstract type OrthogonalFunctionBasis end
abstract type OrthogonalPolynomialBasis{T} <: OrthogonalFunctionBasis end
initialvalue0(::OrthogonalPolynomialBasis{T}, x) where {T} = one(T) # This is often the case! 
(basis::OrthogonalPolynomialBasis)(n::Int, x::Float64) = recursive_eval(basis, n, x)
(basis::OrthogonalPolynomialBasis)(n::Int, x::Float64, r0::Float64) = recursive_eval(basis, n, x/r0)
(basis::OrthogonalPolynomialBasis)(indices::Vector{Int}, x::Float64) = recursive_eval_array(basis, indices, x)
(basis::OrthogonalPolynomialBasis)(indices::Vector{Int}, x::Float64, r0::Float64) = recursive_eval_array(basis, indices, x/r0)


# use the legacy recursive_eval
function recursive_eval(op::OrthogonalPolynomialBasis{T}, n::Int, x::T) where {T} 
    # Intialization 
    n <  0 && return one(T)
    n == 0 && return initialvalue0(op, x)
    n == 1 && return initialvalue1(op, x)

    x0, x1 = initialvalue0(op, x), initialvalue1(op, x)

    # Recursive evaluation
    for i = 1:n - 1
        a1, a2, a3, a4 = recurrence_coefficients(op, i) # typo? n=>i
        x1, x0 = muladd(x, a3, a2)*x1 - a4*x0, x1 
        x1 /= a1
    end

    return x1 
end


struct GeneralizedLaguerreFunction{T, α} <: OrthogonalPolynomialBasis{T} end # Scaled by e^(-r/2) not technically a polynomial but statistfies the recurrence
@inline recurrence_coefficients(::GeneralizedLaguerreFunction{T, α}, n) where {T, α} = n + 1, 2n + α + 1, -1, n + α
initialvalue1(::GeneralizedLaguerreFunction{T, α}, x) where {T, α} = (1 - x + α)*exp(-x / 2)
initialvalue0(::GeneralizedLaguerreFunction{T, α}, x) where {T, α} = exp(- x / 2)
domain(::GeneralizedLaguerreFunction) = (0, Inf)

const LaguerreFunction{T} = GeneralizedLaguerreFunction{T, 0}
(basis::LaguerreFunction)(n::Int, x::Float64, params::Array{Float64}) = basis(n, x, params[1])



struct LegendreP{T} <: OrthogonalPolynomialBasis{T} end
@inline recurrence_coefficients(::LegendreP{T}, n) where T = n + 1, zero(T), 2n + 1, n
initialvalue1(::LegendreP{T}, x) where T = x
domain(::LegendreP) = (-1, 1)
# (basis::LegendreP)(n::Int, x::Float64, ::Any) = basis(n,x)
# (basis::LegendreP)(n::Int, x::Float64, params::Array{Float64}) = basis(n,x)


struct ChebyshevT{T} <: OrthogonalPolynomialBasis{T} end # This is not quite UltraSpherical{1} because of initalvalue1
@inline recurrence_coefficients(::ChebyshevT{T}, n) where T = one(T), zero(T), 2*one(T), one(T)
initialvalue1(::ChebyshevT{T}, x) where {T} = x
domain(::ChebyshevT) = (-1, 1)
(basis::ChebyshevT)(n::Int, x::Float64, ::Any) = basis(n,x)

struct ChebyshevU{T} <: OrthogonalPolynomialBasis{T} end
@inline recurrence_coefficients(::ChebyshevU{T}, n) where T = one(T), zero(T), 2*one(T), one(T)
initialvalue1(::ChebyshevU{T}, x) where {T} = 2x
domain(::ChebyshevU) = (-1, 1)
(basis::ChebyshevU)(n::Int, x::Float64, ::Any) = basis(n,x)

struct Polynomial end
(basis::Polynomial)(n::Int, x::T) where T = x^n
(basis::Polynomial)(n::Int, x::T, x0::T) where T = (x/x0)^n
(basis::Polynomial)(n::Int, x::T, params::Array{T}) where T = basis(n, x, params[1])
