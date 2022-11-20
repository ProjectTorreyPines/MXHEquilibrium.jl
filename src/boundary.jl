struct Boundary{T}
    points::Vector{NTuple{2,T}}

    function Boundary(points::Vector{NTuple{2,T}}) where T
        if first(points) != last(points)
            push!(points,first(points))
        end
        return new{T}(points)
    end
end

const PlasmaBoundary = Boundary
const FluxSurface = Boundary
const Wall = Boundary

Base.length(b::Boundary) = length(b.points)

function Base.show(io::IO, b::Boundary)
    print(io,"Boundary: npoints = $(length(b.points))")
end

Base.broadcastable(b::Boundary) = (b,)

function Boundary(x::Vector,y::Vector)
    Boundary(collect(zip(x,y)))
end

function Base.getproperty(b::Boundary,s::Symbol)
    if s == :points
        return getfield(b,s)
    end
    if s in (:r,:x)
        v = getfield(b,:points)
        return getindex.(v,1)
    end
    if s in (:z,:y)
        v = getfield(b,:points)
        return getindex.(v,2)
    end
    error("type $(typeof(b)) has no field $s")
end

in_boundary(b::Boundary, p) = PolygonOps.inpolygon(p,b.points,in=true,on=true,out=false)
in_boundary(b::Boundary,x,y) = in_boundary(b,(x,y))

const in_vessel = in_boundary
const in_plasma = in_boundary

function flux_surface(M::AbstractEquilibrium, psi::Float64, dx::Float64=0.01, dy::Float64=0.01; n_interp=0, raise_error=true)
    maxis = magnetic_axis(M)
    psimag, psibdry = psi_limits(M)
    if abs(psi) >= abs(psimag)
        return Boundary([maxis[1]],[maxis[2]])
    end

    xlims, ylims = limits(M)
    x = range(xlims...,step=dx)
    y = range(ylims...,step=dy)
    Psi = [M(xx,yy) for xx in x, yy in y]
    bnd = flux_surface(x, y, Psi, psi; maxis = maxis)
    if isnothing(bnd) && raise_error
        throw("Could not trace closed flux surface at ψ=$psi. Note that ψ limits are $(psi_limits(M))")
    end

    if n_interp > 0
        N = length(bnd.r)
        t = 0:(N-1)
        itp_r = linear_interpolation(t, bnd.r, extrapolation_bc = Periodic())
        itp_z = linear_interpolation(t, bnd.z, extrapolation_bc = Periodic())
        τ = range(0,N-1,length=n_interp)
        ri = itp_r.(τ)
        zi = itp_z.(τ)
        return Boundary(ri,zi)
    end

    return bnd
end

function flux_surface(x::AbstractRange, y::AbstractRange, Psi::Matrix, psi::Float64; maxis = ())
    cl = Contour.contour(x, y, Psi, psi)
    # pick a flux surface that closes and isn't far away from Z=0
    for cc in Contour.lines(cl)
        xc, yc = Contour.coordinates(cc)
        if !((xc[1] == xc[end]) & (yc[1] == yc[end]))
            continue
        end
        if length(maxis) > 0
            N = length(xc)
            ym = sum(yc)/N
            xm = sum(xc)/N
            R0, Z0 = maxis
            r = hypot(xm/R0 - 1,ym/R0 - Z0/R0)
            if r > 0.25
                continue
            end
        end
        return Boundary(xc,yc)
    end
    return nothing
end

plasma_boundary(M::AbstractEquilibrium) = flux_surface(M, plasma_boundary_psi(M))

limits(b::Boundary) = extrema(b.r), extrema(b.z)

@memoize LRU(maxsize = 5) function circumference(b::Boundary)
    p = b.points
    return sum(@inbounds norm(p[i+1] .- p[i]) for i=1:(length(p)-1))
end

@memoize LRU(maxsize = 5) function average(b::Boundary, F)
    p = b.points
    return sum(@inbounds norm(p[i+1] .- p[i])*(F(p[i][1],p[i][2])) for i=1:(length(p)-1))/circumference(b)
end

area(s::Boundary) = PolygonOps.area(s.points)

@memoize LRU(maxsize=cache_size) function surface_area(b::Boundary)
    p = b.points
    return 2pi*sum(@inbounds norm(p[i+1] .- p[i])*(p[i+1][1] + p[i][1])/2 for i=1:(length(p)-1))
end

@memoize LRU(maxsize=cache_size) function surface_area_average(b::Boundary, F)
    A = surface_area(b)
    p = b.points
    return (2pi/A)*sum(@inbounds norm(p[i+1] .- p[i])*((p[i+1][1] + p[i][1])/2)*
                   (F(p[i+1][1],p[i+1][2]) + F(p[i][1],p[i][2]))/2 for i=1:(length(p)-1))
end

@memoize LRU(maxsize=cache_size) function surface_area_average(b::Boundary, F::Vector)
    A = surface_area(b)
    p = b.points
    return (2pi/A)*sum(@inbounds norm(p[i+1] .- p[i])*((p[i+1][1] + p[i][1])/2)*
                       (F[i+1] + F[i])/2 for i=1:(length(p)-1))
end

@memoize LRU(maxsize=cache_size) function area_average(b::Boundary, F; dx=0.01, dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    FId = [F(xi,yi)*in_boundary(b,xi,yi) for xi=x,yi=y]
    A = trapz((x,y),FId)

    return A/area(b)
end

@memoize LRU(maxsize=cache_size) function area_average(b::Boundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    Id = [in_boundary(b,xi,yi) for xi=x,yi=y]
    A = trapz((x,y),F.*Id)
    return A/area(b)
end

@memoize LRU(maxsize=cache_size) function volume(b::Boundary; dx=0.01,dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    JId = [xx*in_boundary(b,(xx,yy)) for xx=x,yy=y]

    V = trapz((x,y),JId)

    return 2pi*V
end

@memoize LRU(maxsize=cache_size) function volume_average(b::Boundary, F; dx=0.01, dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    FJId = [F(xx,yy)*xx*in_boundary(b,(xx,yy)) for xx=x,yy=y]

    A = trapz((x,y),FJId)

    return 2pi*A/volume(b)
end

@memoize LRU(maxsize=cache_size) function volume_average(b::Boundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    FJId = [F[ix,iy]*xx[ix]*in_boundary(b,(xx[ix],yy[iy])) for ix in eachindex(x),yy in eachindex(y)]

    A = trapz((x,y),FJId)

    return 2pi*A/volume(b)
end
