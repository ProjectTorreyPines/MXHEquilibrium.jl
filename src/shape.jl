abstract type PlasmaShape{T} end

Base.broadcastable(S::PlasmaShape) = (S,)

function Base.copy(S::T) where T <: PlasmaShape
    SS = T((getfield(S,s) for s in fieldnames(T))...)
    return SS
end

function Base.show(io::IO, S::T) where T<:PlasmaShape
    print(io, "$(T)\n")
    print(io, "  R0 = $(major_radius(S)) [m]\n")
    print(io, "  Z0 = $(elevation(S)) [m]\n")
    print(io, "  ϵ  = $(inv(aspect_ratio(S)))\n")
    print(io, "  κ  = $(elongation(S))\n")
    print(io, "  δ  = $(triangularity(S))\n")
    print(io, "  ζ  = $(squareness(S))\n")
    print(io, "  ξ  = $(ovality(S))\n")
    print(io, "  τ  = $(tilt(S))\n")
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

function boundary(S::PlasmaShape; N=100)
    x,y = shape(S,N=N)
    return Boundary(collect(zip(x[1:end-1],y[1:end-1])))
end

function circumference(S::PlasmaShape; N=100)
    bdry = boundary(S, N=N)
    return circumference(bdry)
end

function scale_aspect(S::PlasmaShape,s)
    SS = copy(S)
    SS.ϵ = s*S.ϵ
    return SS
end

function plasma_geometry(bdry::Boundary)
    rmin,rmax = extrema(bdry.r)
    zmin,zmax = extrema(bdry.z)

    r = (rmax - rmin)/2
    κ = (zmax - zmin)/(2*r)
    R0 = (rmax + rmin)/2
    Z0 = (zmax + zmin)/2

    Ru = bdry.r[argmax(bdry.z)]
    Rl = bdry.r[argmin(bdry.z)]
    δu = (R0 - Ru)/r
    δl = (R0 - Rl)/r

    return R0, Z0, r, κ, δl, δu
end

function plasma_geometry(S::PlasmaShape)
    return plasma_geometry(boundary(S))
end

function plasma_geometry(M::AbstractEquilibrium)
    return plasma_geometry(shape(M))
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
mutable struct MillerShape{T} <: PlasmaShape{T}
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
tilt(S::MillerShape) = zero(S.δ)
ovality(S::MillerShape) = zero(S.δ)
squareness(S::MillerShape) = zero(S.δ)

function shape(S::MillerShape; N=100)
    τ = range(0,2pi,length=N)
    δ₀ = asin(S.δ)
    x = S.R0*(1 .+ S.ϵ .* cos.(τ .+ δ₀*sin.(τ)))
    y = S.Z0 .+ S.R0*(S.ϵ*S.κ*sin.(τ))
    return x,y
end

function (S::MillerShape)(θ)
    δ₀ = asin(S.δ)
    x = S.R0*(1 + S.ϵ * cos(θ + δ₀*sin(θ)))
    y = S.Z0 + S.R0*(S.ϵ*S.κ*sin(θ))
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
mutable struct AsymmetricMillerShape{T} <: PlasmaShape{T}
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
tilt(S::AMShape) = zero(S.δ)
ovality(S::AMShape) = zero(S.δ)
squareness(S::AMShape) = zero(S.δ)

function shape(S::AMShape; N=100)
    τ = range(0,2pi,length=N)
    δ₀l = asin(S.δl)
    δ₀u = asin(S.δu)
    y = S.Z0 .+ S.R0*(S.ϵ*S.κ*sin.(τ))
    δ₀ = [t < S.Z0 ? δ₀l : δ₀u for t in y]
    x = S.R0*(1 .+ S.ϵ .* cos.(τ .+ δ₀ .* sin.(τ)))
    return x,y
end

function (S::AMShape)(θ)
    δ₀l = asin(S.δl)
    δ₀u = asin(S.δu)
    y = S.Z0 + S.R0*(S.ϵ*S.κ*sin(θ))
    x = S.R0*(1 + S.ϵ * cos(θ + (y < S.Z0 ? δ₀l : δ₀u)*sin(θ)))
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
mutable struct TurnbullMillerShape{T} <: PlasmaShape{T}
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
tilt(S::TMShape) = zero(S.δ)
ovality(S::TMShape) = zero(S.δ)
squareness(S::TMShape) = S.ζ

function shape(S::TurnbullMillerShape; N=100)
    τ = range(0,2pi,length=N)
    δ₀ = asin(S.δ)
    x = S.R0*(1 .+ S.ϵ .* cos.(τ .+ δ₀*sin.(τ)))
    y = S.Z0 .+ S.R0*(S.ϵ*S.κ*sin.(τ .+ S.ζ*sin.(2*τ)))
    return x,y
end

function (S::TurnbullMillerShape)(θ)
    δ₀ = asin(S.δ)
    x = S.R0*(1 + S.ϵ * cos(θ + δ₀*sin(θ)))
    y = S.Z0 + S.R0*(S.ϵ*S.κ*sin(θ + S.ζ*sin(2*θ)))
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
mutable struct MillerExtendedHarmonicShape{N, T} <: PlasmaShape{T}
    R0::T           # Major Radius
    Z0::T           # Elevation
    ϵ::T            # inverse aspect ratio a/R0
    κ::T            # Elongation
    c0::T           # Tilt
    c::SVector{N,T} # Cosine coefficients acos.([ovality,...])
    s::SVector{N,T} # Sine coefficients asin.([triangularity,-squareness,...]
end

const MXHShape = MillerExtendedHarmonicShape

MillerExtendedHarmonicShape(N) = MillerExtendedHarmonicShape(0.0, 0.0, 0.0, 0.0, 0.0, zeros(SVector{N}), zeros(SVector{N}))

function MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, c0, c::Vector, s::Vector)
    @assert length(c) == length(s) 
    
    R0, Z0, ϵ, κ, c0 = promote(R0,Z0,ϵ,κ,c0, one(eltype(s)), one(eltype(c)))
    N = length(c)
    c = convert(SVector{N, typeof(R0)}, c)
    s = convert(SVector{N, typeof(R0)}, s)

    MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, c0, c, s)
end

function MillerExtendedHarmonicShape(R0, Z0, ϵ, κ, δ; tilt=zero(κ), c0=tilt, ovality=one(R0), squareness=zero(R0))

    R0, Z0, ϵ, κ, c0, ovality, squareness = promote(R0,Z0,ϵ,κ,c0,ovality,squareness)
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

function shape(S::MillerExtendedHarmonicShape{L,T}; N=100) where {L,T}
    x = zeros(N)
    y = zeros(N)
    r = S.ϵ*S.R0
    θ = range(0,2pi,length=N)
    
    @inbounds for i=1:N
        c_sum = 0.0
        for n=1:L
            c_sum += S.c[n]*cos(n*θ[i])
        end

        s_sum = 0.0
        for n=1:L
            s_sum += S.s[n]*sin(n*θ[i])
        end

        θ_R = θ[i] + S.c0 + c_sum + s_sum
        x[i] = S.R0 + r*cos(θ_R)
        y[i] = S.Z0 + S.κ*r*sin(θ[i])
    end

    return x, y
end

function (S::MillerExtendedHarmonicShape{N,T})(θ) where {N,T}
    r = S.ϵ*S.R0

    c_sum = 0.0
    @inbounds for n=1:N
        c_sum += S.c[n]*cos(n*θ)
    end

    s_sum = 0.0
    @inbounds for n=1:N
        s_sum += S.s[n]*sin(n*θ)
    end

    θ_R = θ + S.c0 + c_sum + s_sum
    x = S.R0 + r*cos(θ_R)
    y = S.Z0 + S.κ*r*sin(θ)

    return x, y
end

function plasma_geometry(S::Union{MillerShape,TurnbullMillerShape})
    return S.R0, S.Z0, S.R0*S.ϵ, S.κ, S.δ, S.δ
end

function plasma_geometry(S::AsymmetricMillerShape)
    return S.R0, S.Z0, S.R0*S.ϵ, S.κ, S.δl, S.δu
end

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

