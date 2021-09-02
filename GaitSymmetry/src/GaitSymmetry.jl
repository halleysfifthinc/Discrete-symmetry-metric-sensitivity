module GaitSymmetry

export limits, reflect

export SymmetryFunc, Sel86, Rob87, Vag92, Plo05, Zif08, Roc14, Que20, Alv20, Alv20b

"""
    (::Type{<:SymmetryFunc})(x, y)

A functor-style type to calculate the symmetry between `x` and `y`.

# Examples
```julia-repl
julia> Zif08(1., 1.1)
3.029234437673622

julia> Rob87(1., 1.1)
9.52380952380953
```
"""
abstract type SymmetryFunc end

"""
    reflect(x)

Return an array where the first half is the reverse of `x` and inverted, and the last half
is `x`, unmodified. The first element of `x` will only occur once.

# Example
```julia-repl
julia> reflect(1//1:3//1)
5-element Array{Rational{Int64},1}:
 1//3
 1//2
 1//1
 2//1
 3//1
```
"""
function reflect(x)
    return [reverse(inv.(x[2:end])); x]
end

function (f::Type{<:SymmetryFunc})(x::T, y::U, args...) where {T, U}
    f(promote(x, y)..., args...)
end

"""
    inv(::Type{<:SymmetryFunc})(s, baseline=1)

Return the inverse of the given symmetry function. The inverse symmetry function calculates
two numbers that produce a similar symmetry, `s`, when used as inputs for the given symmetry
index. `baseline` is an assumption of the value of one of the two numbers.

# Examples
```julia-repl
julia> s = Rob87(2., 2.2)
9.52380952380953

julia> (x2, y2) = inv(Rob87)(s, 2.)
(2.0, 2.1999999999999997)

julia> s ≈ Rob87(x2, y2)
true

julia> (x1, y1) = inv(Rob87)(s1, 1.)
(1.0, 1.0999999999999999)

julia> s ≈ Rob87(x1, y1)
true
```
"""
Base.inv(::Type{<:SymmetryFunc})

"""
    convert(T <: SymmetryFunc, U <: SymmetryFunc) -> anonymous function

Return a function to convert from symmetry metric `U` to `T`

# Examples
```jldoctest
julia> s = Zif08(1., 1.1)
3.029234437673622

julia> convert(Rob87, Zif08)(s)
9.52380952380947

julia> convert(Rob87, Zif08)(s) ≈ Rob87(1,1.1)
true
```
"""
function Base.convert(::Type{T}, ::Type{U}) where {T <: SymmetryFunc, U <: SymmetryFunc}
    return (x) -> T(inv(U)(x)...)
end

function limits(S, T::Type{<:Number}=Float64)
    small_x = S(floatmin(T), floatmax(T))
    small_y = S(floatmax(T), floatmin(T))

    small_x = round(small_x; digits=2)
    small_y = round(small_y; digits=2)

    return minmax(small_x, small_y)
end

"""
    Sel86

Calculate the symmetry of `x` and `y` as the ratio ``y/x``

[1] R. Seliktar and J. Mizrahi, “Some Gait Characteristics of Below-Knee Amputees and Their
Reflection on the Ground Reaction Forces,” Engineering in Medicine, vol. 15, no. 1, pp.
27–34, Jan. 1986.
"""
struct Sel86 <: SymmetryFunc; end
const Ratio = Sel86

function (::Type{Sel86})(x::T, y::T) where T <: Number
    return y/x
end

Base.inv(::Type{Sel86}) = _inv_SymmetryRatio

function _inv_SymmetryRatio(s::T, baseline::T=one(T)) where T <: Number
    return (baseline, s / baseline)
end

"""
    Rob87

Calculate the symmetry of `x` and `y` as defined in [1]:
``
(y - x)/(0.5*(x+y))*100%
``

[1] R. O. Robinson, W. Herzog, and B. M. Nigg, “Use of force platform variables to quantify
the effects of chiropractic manipulation on gait symmetry,” J Manipulative Physiol Ther,
vol. 10, no. 4, pp. 172–176, Aug. 1987.
"""
struct Rob87 <: SymmetryFunc; end

function (::Type{Rob87})(x::T, y::T) where T <: Number
    (y - x)/(x+y)*200
end

Base.inv(::Type{Rob87}) = _inv_Rob87

function _inv_Rob87(s::T, baseline::T=one(T)) where T <: Number
    y = -baseline*(s+200)/(s-200)
    return (baseline, y)
end

"""
    Vag92

``
(y-x)/max(x,y)*100
``

[1] G. Vagenas and B. Hoshizaki, “A Multivariable Analysis of Lower Extremity Kinematic
Asymmetry in Running,” Journal of Applied Biomechanics, vol. 8, no. 1, pp.
11–29, Feb. 1992.

See also: [`Rob87`](@ref)
"""
struct Vag92 <: SymmetryFunc; end

function (::Type{Vag92})(x::T, y::T) where T <: Number
    return (y-x)/max(x,y)*100
end

Base.inv(::Type{Vag92}) = _inv_Vag92

function _inv_Vag92(s::T, baseline::T=one(T)) where T <: Number
    y_xmax = baseline*(1 + s/100)
    y_ymax = baseline/(1 - s/100)

    if Vag92(baseline, y_xmax) ≈ s
        y = y_xmax
    elseif Vag92(baseline, y_ymax) ≈ s
        y = y_ymax
    else
        # Logically should never happen
        throw(DomainError(s, "no possible x/y pair found"))
    end

    return (baseline, y)
end

"""
    Plo05

Calculate the symmetry of `x` and `y`, multiplied by 100% by default [1], as defined in [2]:
``
ln(y/x)
``

[1] M. Plotnik, N. Giladi, Y. Balash, C. Peretz, and J. M. Hausdorff, “Is freezing of gait
in Parkinson’s disease related to asymmetric motor function?,” Ann Neurol., vol. 57, no. 5,
pp. 656–663, May 2005, doi:10.1002/ana.20452.
[2] M. Plotnik, N. Giladi, and J. M. Hausdorff, “A new measure for quantifying the
bilateral coordination of human gait: effects of aging and Parkinson’s disease,” Exp Brain
Res, vol. 181, no. 4, pp. 561–570, Aug. 2007, doi:10.1007/s00221-006-0676-3.
"""
struct Plo05 <: SymmetryFunc; end

function (::Type{Plo05})(x::T, y::T; scaled::Bool=true) where T <: Number
    # !(signbit(x) ⊻ signbit(y)) || throw(DomainError((x,y), "Sign of `x` and `y` must match"))

    ga = log(y/x)
    if scaled
        ga *= 100
    end
    return ga
end

Base.abs(::Type{Plo05}) = _abs_Plo05
_abs_Plo05(x,y) = abs(Plo05(x,y))

Base.inv(::Type{Plo05}) = _inv_Plo05

function _inv_Plo05(s::T, baseline::T=one(T); scaled::Bool=true) where T <: Number
    !signbit(s) || throw(DomainError(s, "Plo05 only produces positive values"))

    if scaled
        y = baseline*exp(s/100)
    else
        y = baseline*exp(s)
    end
    return (baseline, y)
end

"""
    Zif08

Calculate the symmetry of `x` and `y`, originally representing the left and right sides,
respectively, as defined in [1]:
``
((45° - arctan(x/y))/90°)×100%
``

If `(45° - arctan(x/y)) > 90°`, the following equation is used:
``
((45° - arctan(x/y) - 180°)/90°)×100%
``

[1] R. A. Zifchock, I. Davis, J. Higginson, and T. Royer, “The symmetry angle: A novel,
robust method of quantifying asymmetry,” Gait & Posture, vol. 27, no. 4, pp. 622–627, May
2008, doi:10.1016/j.gaitpost.2007.08.006.
"""
struct Zif08 <: SymmetryFunc; end

function (::Type{Zif08})(x::T, y::T) where T <: Number
    a = 45 - atand(x, y)
    a -= (a > 90) ? 180 : 0
    return 100*a/90
end

Base.inv(::Type{Zif08}) = _inv_Zif08

function _inv_Zif08(s::T, baseline::T=one(T)) where T <: Number
    a′ = 90*s/100
    a′ += (a′ + 180 > 90) ? 180 : 0
    a′ = tand(45 - a′)
    return (baseline, baseline / a′)
end

"""
    Roc14(x, y)

``
y - x
``

[1] L. Rochester, B. Galna, S. Lord, and D. Burn, “The nature of dual-task interference
during gait in incident Parkinson’s disease,” Neuroscience, vol. 265, pp. 83–94, 2014, doi:
10.1016/j.neuroscience.2014.01.041.
"""
struct Roc14 <: SymmetryFunc; end

function (::Type{Roc14})(x::T, y::T) where T
    return y - x
end

Base.abs(::Type{Roc14}) = _abs_Roc14
_abs_Roc14(x,y) = abs(Roc14(x,y))

Base.inv(::Type{Roc14}) = _inv_Roc14

function _inv_Roc14(s::T, baseline::T=one(T)) where T
    return (baseline, s + baseline)
end

"""
    Que20

``
(x - y)/(max(0,x,y) - min(0,x,y))
``


[1] R. Queen, L. Dickerson, S. Ranganathan, and D. Schmitt, “A novel method for measuring
asymmetry in kinematic and kinetic variables: The normalized symmetry index,” J Biomech,
vol. 99, p. 109531, 2020, doi:10.1016/j.jbiomech.2019.109531.
"""
struct Que20 <: SymmetryFunc; end

function (::Type{Que20})(x::T, y::T) where T
    return (y - x)/(max(0,x,y) - min(0,x,y))
end

Base.inv(::Type{Que20}) = _inv_Que20

function _inv_Que20(s::T, baseline::T=one(T)) where T
    if s > 0
        y = -baseline/(s-1)
    else
        y = s*baseline + baseline
    end

    return (baseline, y)
end

"""
    Alv20(x, y)

``
(y-x)/sqrt(2*(x^2 + y^2))
``

[1] S. A. Alves, R. M. Ehrig, P. C. Raffalt, A. Bender, G. N. Duda, and A. N. Agres,

“Quantifying Asymmetry in Gait: The Weighted Universal Symmetry Index to Evaluate 3D Ground
Reaction Forces,” Frontiers Bioeng Biotechnology, vol. 8, p. 579511, 2020,
doi:10.3389/fbioe.2020.579511.
"""
struct Alv20 <: SymmetryFunc; end

function (::Type{Alv20})(x::T, y::T, σ::Union{Nothing,T}=nothing) where T
    if σ === nothing
        return (y - x)/hypot(√2*x, √2*y)
    else
        return (y - x)/hypot(√2*x, √2*y)*(1-√2*σ/√(2σ^2 + x^2 + y^2))
    end
end

Base.inv(::Type{Alv20}) = _inv_Alv20

function _inv_Alv20(s::T, baseline::T=one(T)) where T
    y = -baseline*(2s*√(1-s^2)+1)/(2s^2 -1)
    return (baseline, y)
end

struct Alv20b{T} <: SymmetryFunc
    σ::T
end

(alv::Alv20b{T})(x::T, y::T) where T = Alv20(x, y, alv.σ)
(::Type{Alv20b})(x::T, y::T) where T = Alv20(x, y)

end # module
