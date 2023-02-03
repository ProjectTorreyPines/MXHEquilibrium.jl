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

"""
    outer, top, inner, bottom = boundary_extrema(bdry; interp=true)

Returns the outermost, top, innermost, and bottom points on the boundary. Defaults to using interpolation to find extrema
"""
function boundary_extrema(bdry::Boundary; interp=true)

    if interp
        t = range(0,1,length=length(bdry.points))
        rz = hcat(bdry.r,bdry.z)
        itp = scale(Interpolations.interpolate(rz, (BSpline(Cubic(Periodic(OnGrid()))), NoInterp())), t, 1:2)

        t0 = t[argmax(bdry.r)]
        res = Optim.optimize(x -> -itp(x[1],1), [0.0], [1.0], [t0], Optim.Fminbox(Optim.GradientDescent()))
        rmax = -Optim.minimum(res)
        ir_max = Optim.minimizer(res)[1]
        z_rmax = itp(ir_max,2)

        t0 = t[argmin(bdry.r)]
        res = Optim.optimize(x -> itp(x[1],1), [0.0], [1.0], [t0], Optim.Fminbox(Optim.GradientDescent()))
        rmin = Optim.minimum(res)
        ir_min = Optim.minimizer(res)[1]
        z_rmin = itp(ir_min,2)

        t0 = t[argmax(bdry.z)]
        res = Optim.optimize(x -> -itp(x[1],2), [0.0], [1.0], [t0], Optim.Fminbox(Optim.GradientDescent()))
        zmax = -Optim.minimum(res)
        iz_max = Optim.minimizer(res)[1]
        r_zmax = itp(iz_max,1)

        t0 = t[argmin(bdry.z)]
        res = Optim.optimize(x -> itp(x[1],2), [0.0], [1.0], [t0], Optim.Fminbox(Optim.GradientDescent()))
        zmin = Optim.minimum(res)
        iz_min = Optim.minimizer(res)[1]
        r_zmin = itp(iz_min,1)
    else
        rmin, ir_min = findmin(bdry.r)
        rmax, ir_max = findmax(bdry.r)
        zmin, iz_min = findmin(bdry.z)
        zmax, iz_max = findmax(bdry.z)

        r_zmax = bdry.r[iz_max]
        r_zmin = bdry.r[iz_min]
        z_rmax = bdry.z[ir_max]
        z_rmin = bdry.z[ir_min]
    end

    outer = (rmax,z_rmax)
    inner = (rmin,z_rmin)
    top = (r_zmax, zmax)
    bottom = (r_zmin, zmin)

    return outer, top, inner, bottom
end

function flux_surface(M::AbstractEquilibrium, psi::Float64, dx::Float64=0.01, dy::Float64=0.01; n_interp=0, raise_error=true, kwargs...)
    maxis = magnetic_axis(M)
    psimag, psibdry = psi_limits(M)
    if abs(psi) >= abs(psimag)
        return Boundary([maxis[1]],[maxis[2]])
    end

    xlims, ylims = limits(M; kwargs...)
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

plasma_boundary(M::AbstractEquilibrium; kwargs...) = flux_surface(M, psi_boundary(M); kwargs...)

function limits(b::Boundary; pad=0.2)
    xlims, ylims = extrema(b.r), extrema(b.z)
    xpad = xlims[2] - xlims[1]
    ypad = ylims[2] - ylims[1]
    xlims = (max(xlims[1] - xpad * pad, 0.0), xlims[2] + xpad * pad)
    ylims = (ylims[1] - ypad * pad, ylims[2] + ypad * pad)
    return xlims, ylims
end

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
