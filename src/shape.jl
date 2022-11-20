abstract type PlasmaShape{T} end

Base.broadcastable(S::PlasmaShape) = (S,)

function Base.copy(S::T) where T <: PlasmaShape
    SS = T((getfield(S,s) for s in fieldnames(T))...)
    return SS
end

function (S::PlasmaShape)(ρ::T,θ::T,ζ::T) where T<:Number
    ρ > minor_radius(S) && throw(DomainError("ρ is larger than the minor radius"))
    r, z = S(ρ,θ)
    x = r*cos(ζ)
    y = r*sin(ζ)
    return x,y,z
end

function (S::PlasmaShape)(x::T) where T <:AbstractVector
    N = length(x)
    if N == 2
        return collect(S(x[1],x[2]))
    elseif N == 3
        return collect(S(x[1],x[2],x[3]))
    else
        throw(ArgumentError("Length of array can only be 2 or 3"))
    end
end

"""
3 element u: Volume element
2 element u: Internal Area element
"""
function dVA(S::PlasmaShape,u::T) where T<:AbstractVector
    @assert length(u) in (2,3)
    abs(det(ForwardDiff.jacobian(S,u)))
end

"""
Surface area element of a plasma shape at a given position
"""
function dS(S::PlasmaShape,u::T) where T<:AbstractVector{R} where R
    N = length(u)
    @assert N == 2
    d1 = zeros(R,3)
    d2 = zeros(R,3)
    r = minor_radius(S)
    for i=1:3
        d1[i] = ForwardDiff.derivative(t->S([r,t,u[2]])[i], u[1])
        d2[i] = ForwardDiff.derivative(t->S([r,u[1],t])[i], u[2])
    end
    return norm(cross(d1,d2))
end

"""
Line element of a a plasma shape at a given position
"""
function dL(S::PlasmaShape,u::T) where T<:Number
    d = zeros(2)
    r = minor_radius(S)
    for i=1:2
        d[i] = ForwardDiff.derivative(t->S([r,t])[i], u)
    end
    return norm(d)
end

"""
Integrate a function F which takes a scalar or vector over the given plasma shape

Example:

```
I = integrate(x->1,S,:line,(0,2pi)) # circumference
I = integrate(x->1,S,:area,(0,minor_radius(S)),(0,2pi)) # area
I = integrate(x->1,S,:surface,(0,2pi),(0,2pi)) # surface area
I = integrate(x->1,S,:volume,(0,minor_radius(S)),(0,2pi),(0,2pi))
```
"""
function integrate(f::Function, S::PlasmaShape, type::Symbol, bds::Vararg{NTuple{2},N}; kwargs...) where N
    kws = pairs((rtol=1e-3, atol=1e-3, kwargs...))
    if N == 1 && type == :line
        θ_min, θ_max = bds[1]
        (F,err) = hquadrature(x->f(x)*dL(S,x),θ_min,θ_max; kws...)
    elseif N == 2 && type == :surface
        θ_min, θ_max = bds[1]
        ζ_min, ζ_max = bds[2]
        (F,err) = hcubature(x->f(x)*dS(S,x),(θ_min,ζ_min),(θ_max,ζ_max); kws...)
    elseif N == 2 && type == :area
        r_min, r_max = bds[1]
        θ_min, θ_max = bds[2]
        (F,err) = hcubature(x->f(x)*dVA(S,x),(r_min, θ_min),(r_max, θ_max); kws...)
    elseif N == 3 && type == :volume
        r_min, r_max = bds[1]
        θ_min, θ_max = bds[2]
        ζ_min, ζ_max = bds[3]
        (F,err) = hcubature(x->f(x)*dVA(S,x),(r_min,θ_min,ζ_min),(r_max,θ_max,ζ_max); kws...)
    else
        throw(ArgumentError("Unsupported Integral Type: $type : Supported types: :line ∫dθ, :area ∬drdθ, :surface ∬dθdζ, :volume ∭drdθdζ"))
    end

    return F
end

@memoize LRU(maxsize=cache_size) function circumference(S::PlasmaShape; kwargs...)
    integrate(x->1, S, :line, (0.0,2pi); kwargs...)
end

@memoize LRU(maxsize=cache_size) function area(S::PlasmaShape; kwargs...)
    integrate(x->1, S, :area, (0.0,minor_radius(S)), (0.0,2pi); kwargs...)
end

@memoize LRU(maxsize=cache_size) function surface_area(S::PlasmaShape; kwargs...)
    integrate(x->1,S,:surface, (0.0,2pi), (0.0,2pi); kwargs...)
end

@memoize LRU(maxsize=cache_size) function volume(S::PlasmaShape; kwargs...)
    integrate(x->1, S, :volume, (0.0,minor_radius(S)), (0.0,2pi), (0.0,2pi); kwargs...)
end

@memoize LRU(maxsize=cache_size) function average(f::Function, S::PlasmaShape, type::Symbol; kwargs...)

    if type == :line
        F_bar = integrate(f, S, :line, (0.0,2pi); kwargs...)/circumference(S; kwargs...)
    elseif type == :surface
        F_bar = integrate(f, S, :surface, (0.0,2pi),(0.0,2pi); kwargs...)/surface_area(S; kwargs...)
    elseif type == :area
        F_bar = integrate(f, S, :area, (0.0, minor_radius(S)), (0.0, 2pi); kwargs...)/area(S; kwargs...)
    elseif type == :volume
        F_bar = integrate(f, S, :volume, (0.0, minor_radius(S)), (0.0, 2pi), (0, 2pi); kwargs...)/volume(S; kwargs...)
    else
        throw(ArgumentError("Unsupported Average Type: $type : Supported types: :line ∫dθ, :area ∬drdθ, :surface ∬dθdζ, :volume ∭drdθdζ"))
    end
end

function limits(s::PlasmaShape, x_point=nothing)
    xlims = (0.8*s(pi)[1], 1.2*s(0.0)[1])
    ylims = (1.2*s(3pi/2)[2], 1.2*s(pi/2)[2])
    if x_point !== nothing
        xlims = (min(xlims[1],x_point[1]), max(xlims[2],x_point[1]))
        ylims = (min(ylims[1],x_point[2]), max(ylims[2],x_point[2]))
    end
    return xlims, ylims
end

function plasma_boundary(S::PlasmaShape; kwargs...)
    x,y = shape(S; kwargs...)
    return PlasmaBoundary(collect(zip(x[1:end-1],y[1:end-1])))
end

function scale_aspect(S::PlasmaShape,s)
    SS = copy(S)
    SS.ϵ = s*S.ϵ
    return SS
end

function Base.getproperty(S::PlasmaShape,s::Symbol)
    if s == :_get_x
        return x -> S(x)[1]
    elseif s == :_get_y
        return x -> S(x)[2]
    end
    return getfield(S,s)
end

"""
MillerShape Structure

Defines the Miller Plasma Shape Parameterization

Fields:\\
`R0` - Major Radius [m]\\
`Z0` - Elevation [m]\\
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`δ`  - Triangularity
"""
struct MillerShape{T} <: PlasmaShape{T}
    R0::T  # Major Radius [m]
    Z0::T  # Elevation
    ϵ::T   # Inverse Aspect Ratio a/R0 (a = minor radius)
    κ::T   # Elongation
    δ::T   # Triangularity
end

const MShape = MillerShape

MillerShape() = MillerShape(0.0, 0.0, 0.0, 0.0, 0.0)

function MillerShape(R0,Z0,ϵ,κ,δ)
    MillerShape(promote(R0,Z0,ϵ,κ,δ)...)
end

# Miller Shape API
aspect_ratio(S::MillerShape) = inv(S.ϵ)
elongation(S::MillerShape) = S.κ
major_radius(S::MillerShape) = S.R0
minor_radius(S::MillerShape) = S.R0*S.ϵ
elevation(S::MillerShape) = S.Z0
triangularity(S::MillerShape) = S.δ

function Base.show(io::IO, S::MillerShape)
    print(io, "$(typeof(S))\n")
    print(io, "  R0 = $(round(major_radius(S),digits=3)) [m]\n")
    print(io, "  Z0 = $(round(elevation(S),digits=3)) [m]\n")
    print(io, "  ϵ  = $(round(inv(aspect_ratio(S)),digits=3))\n")
    print(io, "  κ  = $(round(elongation(S),digits=3))\n")
    print(io, "  δ  = $(round(triangularity(S),digits=3))")
end

function m_rz(r, θ, R0, Z0, κ, δ)
    δ₀ = asin(δ)
    x = R0 + r * cos(θ + δ₀ * sin(θ))
    y = Z0 + r * κ * sin(θ)
    return x, y
end

function shape(S::MillerShape; N=100, normed=false)
    r = S.ϵ*S.R0

    x = zeros(N)
    y = zeros(N)
    τ = range(0,2pi,length=N)
    @inbounds for i=1:N
        x[i], y[i] = m_rz(r, τ[i], S.R0, S.Z0, S.κ, S.δ)
    end
    if !normed
        return x, y
    else
        return x/S.R0, (y .- S.Z0)/S.R0
    end
end

function (S::MillerShape)(θ::T) where T<:Number
    r = S.ϵ*S.R0
    x, y = m_rz(r, θ, S.R0, S.Z0, S.κ, S.δ)
    return x,y
end

function (S::MillerShape)(r::T,θ::T) where T<:Number
    x, y = m_rz(r, θ, S.R0, S.Z0, S.κ, S.δ)
    return x,y
end

"""
AsymmetricMillerShape Structure

Defines the Asymmetric Miller Plasma Shape Parameterization

Fields:\\
`R0` - Major Radius [m]\\
`Z0` - Elevation [m]\\
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`δl` - Lower Triangularity\\
`δu` - Upper Triangularity
"""
struct AsymmetricMillerShape{T} <: PlasmaShape{T}
    R0::T  # Major Radius [m]
    Z0::T  # Elevation
    ϵ::T   # Inverse Aspect Ratio a/R0 (a = minor radius)
    κ::T   # Elongation
    δl::T  # Lower Triangularity
    δu::T  # Lower Triangularity
end

const AMShape = AsymmetricMillerShape

AsymmetricMillerShape() = AsymmetricMillerShape(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function AsymmetricMillerShape(R0,Z0,ϵ,κ,δl,δu)
    MillerShape(promote(R0,Z0,ϵ,κ,δl,δu)...)
end

aspect_ratio(S::AMShape) = inv(S.ϵ)
elongation(S::AMShape) = S.κ
major_radius(S::AMShape) = S.R0
minor_radius(S::AMShape) = S.R0*S.ϵ
elevation(S::AMShape) = S.Z0
triangularity(S::AMShape) = (S.δl,S.δu)

function Base.show(io::IO, S::AMShape)
    print(io, "$(typeof(S))\n")
    print(io, "  R0 = $(round(major_radius(S),digits=3)) [m]\n")
    print(io, "  Z0 = $(round(elevation(S),digits=3)) [m]\n")
    print(io, "  ϵ  = $(round(inv(aspect_ratio(S)),digits=3))\n")
    print(io, "  κ  = $(round(elongation(S),digits=3))\n")
    print(io, "  δ  = $(round.(triangularity(S),digits=3))")
end

function am_rz(r, θ, R0, Z0, κ, δl, δu)
    δ₀l = asin(δl)
    δ₀u = asin(δu)
    y = Z0 + r * κ * sin(θ)
    δ₀ = y < Z0 ? δ₀l : δ₀u
    x = R0 + r * cos(θ + δ₀ * sin(θ))
    return x, y
end

function shape(S::AMShape; N=100, normed=false)
    r = S.ϵ*S.R0

    x = zeros(N)
    y = zeros(N)
    τ = range(0,2pi,length=N)
    @inbounds for i=1:N
        x[i], y[i] = am_rz(r, τ[i], S.R0, S.Z0, S.κ, S.δl, S.δu)
    end
    if !normed
        return x, y
    else
        return x/S.R0, (y .- S.Z0)/S.R0
    end
end

function (S::AMShape)(θ::T) where T<:Number
    r = S.ϵ*S.R0
    x, y = am_rz(r, θ, S.R0, S.Z0, S.κ, S.δl, S.δu)
    return x,y
end

function (S::AMShape)(r::T,θ::T) where T<:Number
    x, y = am_rz(r, θ, S.R0, S.Z0, S.κ, S.δl, S.δu)
    return x,y
end

"""
TurnbullMillerShape Structure

Defines the Turnbull-Miller Plasma Shape Parameterization\\
> Turnbull, A. D., et al. "Improved magnetohydrodynamic stability through optimization of higher order moments in cross-section shape of tokamaks." Physics of Plasmas 6.4 (1999): 1113-1116.

Fields:\\
`R0` - Major Radius [m]\\
`Z0` - Elevation [m]\\
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`δ`  - Triangularity\\
`ζ`  - Squareness
"""
struct TurnbullMillerShape{T} <: PlasmaShape{T}
    R0::T  # Major Radius [m]
    Z0::T  # Elevation
    ϵ::T   # Inverse Aspect Ratio a/R0 (a = minor radius)
    κ::T   # Elongation
    δ::T   # Triangularity
    ζ::T   # Squareness
end

const TMShape = TurnbullMillerShape

TurnbullMillerShape() = TurnbullMillerShape(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function TurnbullMillerShape(R0,Z0,ϵ,κ,δ,ζ)
    TurnbullMillerShape(promote(R0,Z0,ϵ,κ,δ,ζ)...)
end

aspect_ratio(S::TMShape) = inv(S.ϵ)
elongation(S::TMShape) = S.κ
major_radius(S::TMShape) = S.R0
minor_radius(S::TMShape) = S.R0*S.ϵ
elevation(S::TMShape) = S.Z0
triangularity(S::TMShape) = S.δ
squareness(S::TMShape) = S.ζ

function Base.show(io::IO, S::TMShape)
    print(io, "$(typeof(S))\n")
    print(io, "  R0 = $(round(major_radius(S),digits=3)) [m]\n")
    print(io, "  Z0 = $(round(elevation(S),digits=3)) [m]\n")
    print(io, "  ϵ  = $(round(inv(aspect_ratio(S)),digits=3))\n")
    print(io, "  κ  = $(round(elongation(S),digits=3))\n")
    print(io, "  δ  = $(round(triangularity(S),digits=3))\n")
    print(io, "  ζ  = $(round(squareness(S),digits=3))")
end

function tm_rz(r, θ, R0, Z0, κ, δ, ζ)
    δ₀ = asin(δ)
    x = R0 + r * cos(θ + δ₀*sin(θ))
    y = Z0 + r * κ * sin(θ + ζ*sin(2*θ))
    return x, y
end

function shape(S::TMShape; N=100, normed=false)
    r = S.ϵ*S.R0

    x = zeros(N)
    y = zeros(N)
    τ = range(0,2pi,length=N)
    @inbounds for i=1:N
        x[i], y[i] = tm_rz(r, τ[i], S.R0, S.Z0, S.κ, S.δ, S.ζ)
    end
    if !normed
        return x, y
    else
        return x/S.R0, (y .- S.Z0)/S.R0
    end
end

function (S::TMShape)(θ::T) where T<:Number
    r = S.ϵ*S.R0
    x, y = tm_rz(r, θ, S.R0, S.Z0, S.κ, S.δ, S.ζ)
    return x,y
end

function (S::TMShape)(r::T,θ::T) where T<:Number
    x, y = tm_rz(r, θ, S.R0, S.Z0, S.κ, S.δ, S.ζ)
    return x,y
end

function TurnbullMillerShape(S::MillerShape)
    TurnbullMillerShape(S.R0, S.Z0, S.ϵ, S.κ, S.δ, zero(S.δ))
end

"""
MillerExtendedHarmonicShape Structure

Defines the Miller Extended Harmonic Plasma Shape Parameterization\\
> Arbon, Ryan, Jeff Candy, and Emily A. Belli. "Rapidly-convergent flux-surface shape parameterization." Plasma Physics and Controlled Fusion 63.1 (2020): 012001.

Fields:\\
`R0` - Major Radius [m]\\
`Z0` - Elevation [m]\\
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`c0` - Tilt\\
`c`  - Cosine coefficients i.e. [ovality,...]\\
`s`  - Sine coefficients i.e. [asin(triangularity), squareness,...])
"""
struct MillerExtendedHarmonicShape{N, T} <: PlasmaShape{T}
    R0::T           # Major Radius
    Z0::T           # Elevation
    ϵ::T            # inverse aspect ratio a/R0
    κ::T            # Elongation
    c0::T           # Tilt
    c::SVector{N,T} # Cosine coefficients acos.([ovality,...])
    s::SVector{N,T} # Sine coefficients asin.([triangularity,-squareness,...]
end

const MXHShape = MillerExtendedHarmonicShape

function MillerExtendedHarmonicShape(N)
    if N == -2
        return MillerShape()
    elseif N == -1
        return AsymmetricMillerShape()
    end
    return MXHShape(0.0, 0.0, 0.0, 0.0, 0.0, zeros(SVector{N}), zeros(SVector{N}))
end

function MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, c0, c::Vector, s::Vector)
    @assert length(c) == length(s) 
    
    R0, Z0, ϵ, κ, c0 = promote(R0,Z0,ϵ,κ,c0, one(eltype(s)), one(eltype(c)))
    N = length(c)
    if N == 0 && c0 == zero(R0)
        return MillerShape(R0,Z0,ϵ,κ,c0)
    end
    c = convert(SVector{N, typeof(R0)}, c)
    s = convert(SVector{N, typeof(R0)}, s)

    MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, c0, c, s)
end

function MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, δ; tilt=zero(κ), c0=tilt, ovality=one(R0), squareness=zero(R0))

    R0, Z0, ϵ, κ, c0, ovality, squareness = promote(R0,Z0,ϵ,κ,c0,ovality,squareness)
    Z = zero(R0)
    if δ == Z && c0 == Z && ovality == Z && squareness == Z
        return MillerShape(R0,Z0,ϵ,κ,δ)
    end

    c = SVector(ovality, zero(ovality))
    s = SVector(asin(δ), -squareness)
    MillerExtendedHarmonicShape(R0,Z0,ϵ,κ,c0,c,s)
end

function MillerExtendedHarmonicShape(S::MillerShape)
    MXHShape(S.R0,S.Z0,S.ϵ,S.κ,zero(S.δ), SVector(zero(S.δ)), SVector(asin.(S.δ)))
end

function MillerExtendedHarmonicShape(S::TurnbullMillerShape)
    MXHShape(S.R0, S.Z0, S.ϵ, S.κ, zero(S.δ), SVector(zero(S.δ),zero(S.δ)), SVector(asin.(S.δ),-S.ζ))
end

aspect_ratio(S::MXHShape) = inv(S.ϵ)
elongation(S::MXHShape) = S.κ
major_radius(S::MXHShape) = S.R0
minor_radius(S::MXHShape) = S.R0*S.ϵ
elevation(S::MXHShape) = S.Z0
tilt(S::MXHShape) = S.c0
ovality(S::MXHShape) = S.c[1]
triangularity(S::MXHShape) = sin(S.s[1])
# MXH squareness differs from TurnbullMiller, using MXH paper's definition
squareness(S::MXHShape) = S.s[2]

function Base.show(io::IO, S::MXHShape)
    print(io, "$(typeof(S))\n")
    print(io, "  R0 = $(round(major_radius(S),digits=3)) [m]\n")
    print(io, "  Z0 = $(round(elevation(S),digits=3)) [m]\n")
    print(io, "  ϵ  = $(round(inv(aspect_ratio(S)),digits=3))\n")
    print(io, "  κ  = $(round(elongation(S),digits=3))\n")
    print(io, "  c₀ = $(round(tilt(S),digits=3))\n")
    print(io, "  c  = $(round.(S.c,digits=3))\n")
    print(io, "  s  = $(round.(S.s,digits=3))")
end

function mxh_rz(r, θ, R0, Z0, κ, c0, c::SVector{N}, s::SVector{N}) where N

    c_sum = 0.0
    @inbounds for n=1:N
        c_sum += c[n]*cos(n*θ)
    end

    s_sum = 0.0
    @inbounds for n=1:N
        s_sum += s[n]*sin(n*θ)
    end

    θ_R = θ + c0 + c_sum + s_sum
    x = R0 + r*cos(θ_R)
    y = Z0 + κ*r*sin(θ)

    return x, y
end

function shape(S::MXHShape; N=100, normed=false)
    r = S.ϵ*S.R0

    x = zeros(N)
    y = zeros(N)
    θ = range(0,2pi,length=N)
    @inbounds for i=1:N
        x[i], y[i] = mxh_rz(r, θ[i], S.R0, S.Z0, S.κ, S.c0, S.c, S.s)
    end
    if !normed
        return x, y
    else
        return x/S.R0, (y .- S.Z0)/S.R0
    end
end

function (S::MXHShape)(θ::T) where T<:Number
    r = S.ϵ*S.R0
    return mxh_rz(r, θ, S.R0, S.Z0, S.κ, S.c0, S.c, S.s)
end

function (S::MXHShape)(r::T,θ::T) where T<:Number
    return mxh_rz(r, θ, S.R0, S.Z0, S.κ, S.c0, S.c, S.s)
end

"""
Luce Plasma Shape as described in:

> "An analytic functional form for characterization and generation of axisymmetric plasma boundaries",\\
TC Luce, Plasma Phys. Control. Fusion 55 (2013) http://dx.doi.org/10.1088/0741-3335/55/9/095009

Fields:\\
  `R0` - Major Radius [m]\\
  `Z0` - Elevation [m]\\
  `r`  - Minor Radius [m]\\
  `Zᵣₘ` - Z(Rₘₐₓ) [m]\\
  `κ`  - Lower and Upper Elongation\\
  `δ`  - Lower and Upper Triangulation\\
  `ζ`  - Squareness for the I,II,III,IV quadrants
"""
struct LuceShape{T} <: PlasmaShape{T}
    R0::T          # Major Radius
    Z0::T          # Elevation
    r::T           # Minor Radius
    Zᵣₘ::T         # Zoff = Z(Rₘₐₓ)
    κ::NTuple{2,T} # Lower and Upper Elongation
    δ::NTuple{2,T} # Lower and Upper Triangularity
    ζ::NTuple{4,T} # Squareness for the 4 quadrants ζ_(uo,ui,li,lo)
end

const LShape = LuceShape

LuceShape() = LuceShape(0.0, 0.0, 0.0, 0.0, (0.0,0.0), (0.0,0.0), (0.0,0.0,0.0,0.0))

function LuceShape(R0,Z0,r,Zrm,κ::NTuple{2},δ::NTuple{2},ζ::NTuple{4})
    R0,Z0,r,Zrm = promote(R0,Z0,r,Zrm)
    κ = convert.(typeof(R0), κ)
    δ = convert.(typeof(R0), δ)
    ζ = convert.(typeof(R0), ζ)
    LuceShape(R0,Z0,r,Zrm,κ,δ,ζ)
end

function LuceShape(G::PlasmaGeometricParameters)
    LuceShape(getfield.(G,fieldnames(typeof(G)))...)
end

function Base.show(io::IO, G::LuceShape)
    print(io, "$(typeof(G))\n")
    print(io, "  R0 = $(round(G.R0,digits=3)) [m]\n")
    print(io, "  Z0 = $(round(G.Z0,digits=3)) [m]\n")
    print(io, "  r  = $(round(G.r,digits=3))  [m]\n")
    print(io, "  Zᵣₘ= $(round.(G.Zᵣₘ,digits=3)) [m]\n")
    print(io, "  κ  = $(round.(G.κ,digits=3))\n")
    print(io, "  δ  = $(round.(G.δ,digits=3))\n")
    print(io, "  ζ  = $(round.(G.ζ,digits=3))")
end

aspect_ratio(S::LShape) = S.r/S.R0
elongation(S::LShape) = S.κ
major_radius(S::LShape) = S.R0
minor_radius(S::LShape) = S.r
elevation(S::LShape) = S.Z0
triangularity(S::LShape) = S.δ
squareness(S::LShape) = S.ζ

function superellipse(t,A,B,n)
    st, ct = sincos(t)
    x = abs(ct)^(2/n) * A*sign(ct)
    y = abs(st)^(2/n) * B*sign(st)
    return x, y
end

function luce_rz(r, θ, R0, Z0, Zrm, κ::NTuple{2}, δ::NTuple{2}, ζ::NTuple{4})

    θ = mod2pi(θ)
    if 0 <= θ < pi/2
        t = θ
        A = r*(1 + δ[2])
        B = κ[2]*r
        n = -log(2)/log(inv(sqrt(2)) + ζ[1]*(1 - inv(sqrt(2))))

        x, y = superellipse(t,A,B,n)
        R = x + r*(inv(r/R0) - δ[2])
        Z = y + Zrm
    elseif pi/2 <= θ < pi
        t = pi - θ
        A = r*(1 - δ[2])
        B = κ[2]*r
        n = -log(2)/log(inv(sqrt(2)) + ζ[2]*(1 - inv(sqrt(2))))

        x, y = superellipse(t,A,B,n)
        R = r*(inv(r/R0) - δ[2]) - x
        Z = y + Zrm
    elseif pi <= θ < 3pi/2
        t = θ - pi
        A = r*(1 - δ[1])
        B = κ[1]*r
        n = -log(2)/log(inv(sqrt(2)) + ζ[3]*(1 - inv(sqrt(2))))

        x, y = superellipse(t,A,B,n)
        R = r*(inv(r/R0) - δ[1]) - x
        Z = Zrm - y
    else
        t = 2pi - θ
        A = r*(1 + δ[1])
        B = κ[1]*r
        n = -log(2)/log(inv(sqrt(2)) + ζ[4]*(1 - inv(sqrt(2))))

        x, y = superellipse(t,A,B,n)
        R = x + r*(inv(r/R0) - δ[1])
        Z = Zrm - y
    end
    return R, Z
end

function shape(S::LuceShape; N=100, normed=false)
    x = zeros(N)
    y = zeros(N)
    θ = range(0,2pi,length=N)
    @inbounds for i=1:N
        x[i], y[i] = luce_rz(S.r, θ[i], S.R0, S.Z0, S.Zᵣₘ, S.κ, S.δ, S.ζ)
    end
    if !normed
        return x, y
    else
        return x/S.R0, (y .- S.Z0)/S.R0
    end
end

function (S::LuceShape)(θ::T) where T<:Number
    return luce_rz(S.r, θ, S.R0, S.Z0, S.Zᵣₘ, S.κ, S.δ, S.ζ)
end

function (S::LuceShape)(r::T,θ::T) where T<:Number
    return luce_rz(r, θ, S.R0, S.Z0, S.Zᵣₘ, S.κ, S.δ, S.ζ)
end

# --- special cases ----
function plasma_geometry(S::PlasmaShape)
    return plasma_geometry(plasma_boundary(S))
end

function plasma_geometry(S::Union{MShape,TMShape})
    Z = zero(S.δ)
    return PlasmaGeometricParameters(S.R0,S.Z0,S.R0*S.ϵ,S.Z0,(S.κ,S.κ), (S.δ, S.δ),(Z,Z,Z,Z))
end

function plasma_geometry(S::AMShape)
    Z = zero(S.δ)
    return PlasmaGeometricParameters(S.R0,S.Z0,S.R0*S.ϵ, S.Z0, (S.κ,S.κ), (S.δl, S.δu),(Z,Z,Z,Z))
end

#--- Curvature calculation via AutoDiff ---
_d1x(S) = (S._get_x)'
_d1y(S) = (S._get_y)'
_d2x(S) = (S._get_x)''
_d2y(S) = (S._get_y)''

function curvature(S::PlasmaShape,θ)
    xp  = _d1x(S)(θ)
    yp  = _d1y(S)(θ)
    xpp = _d2x(S)(θ)
    ypp = _d2y(S)(θ)

    κ = abs(yp*xpp - ypp*xp)/(xp^2 + yp^2)^1.5
    return κ
end

