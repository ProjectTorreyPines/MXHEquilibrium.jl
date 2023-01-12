"""
Plasma Geometry Parameters as described in:

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
struct PlasmaGeometricParameters{T}
    R0::T          # Major Radius
    Z0::T          # Elevation
    r::T           # Minor Radius
    Zᵣₘ::T         # Zoff = Z(Rₘₐₓ)
    κ::NTuple{2,T} # Lower and Upper Elongation
    δ::NTuple{2,T} # Lower and Upper Triangularity
    ζ::NTuple{4,T} # Squareness for the 4 quadrants ζ_(uo,ui,li,lo)
end

Base.broadcastable(G::PlasmaGeometricParameters) = (G,)

function Base.show(io::IO, G::PlasmaGeometricParameters)
    print(io, "$(typeof(G))\n")
    print(io, "  R0 = $(round(G.R0,digits=3)) [m]\n")
    print(io, "  Z0 = $(round(G.Z0,digits=3)) [m]\n")
    print(io, "  r  = $(round(G.r,digits=3))  [m]\n")
    print(io, "  Zᵣₘ= $(round.(G.ζ,digits=3)) [m]\n")
    print(io, "  κ  = $(round.(G.κ,digits=3))\n")
    print(io, "  δ  = $(round.(G.δ,digits=3))\n")
    print(io, "  ζ  = $(round.(G.ζ,digits=3))")
end

function intersection(a, b, origin, p1, p2)
    m = (p2[2] - p1[2])/(p2[1] - p1[1])
    c = p1[2] - m*p1[1]

    x0, y0 = origin

    # intersection between a line and an ellipse
    x1 = (-sqrt(a^2*b^2*(a^2*m^2 + b^2 - c^2 - 2*c*m*x0 +
          2*c*y0 - m^2 * x0^2 + 2*m*x0*y0 - y0^2)) +
          a^2*m*(y0 - c) + b^2 * x0)/(a^2 * m^2 + b^2)
    x2 = (sqrt(a^2*b^2*(a^2*m^2 + b^2 - c^2 - 2*c*m*x0 +
          2*c*y0 - m^2*x0^2 + 2*m*x0*y0 - y0^2)) +
          a^2*m*(y0 - c) + b^2*x0)/(a^2*m^2 + b^2)

    y1 = m*x1 + c
    y2 = m*x2 + c

    return (x1,y1), (x2,y2)
end

function intersection(B::Boundary, p1, p2; tol=1e-6)

    p1_in = in_boundary(B,p1)
    p2_in = in_boundary(B,p2)

    if p2_in && !p1_in
        p1, p2 = p2, p1
    elseif p1_in && p2_in
        error("Initial points inside boundary")
    elseif !p1_in && !p2_in
        error("Initial points outside boundary")
    end

    L = Inf
    while tol < L
        pm = (p1 .+ p2)./2
        pm_in = in_boundary(B,pm)
        if pm_in
            p1 = pm
        else
            p2 = pm
        end
        L = norm(p2 .- p1)
    end
    return (p1 .+ p2) ./ 2
end

function quadrant(p)
    p[1] > 0 && p[2] > 0 && return 1
    p[1] < 0 && p[2] > 0 && return 2
    p[1] < 0 && p[2] < 0 && return 3
    p[1] > 0 && p[2] < 0 && return 4
end

function quadrant_squareness(bdry,A,B,p1,p2,quad)

    ## intersection with ellipse
    pi1, pi2 = intersection(A,B,p1,p1,p2)
    if quadrant(pi1 .- p1) == quad
        L_OC = norm(pi1 .- p1)
        L_CE = norm(p2 .- pi1)
        p_e = pi1
    else
        L_OC = norm(pi2 .- p1)
        L_CE = norm(p2 .- pi2)
        p_e = pi2
    end

    ## intersection with boundary
    if in_boundary(bdry,p_e)
        p_bi = intersection(bdry,p_e,p2)
    else
        p_bi = intersection(bdry,p1,p_e)
    end
    L_OD = norm(p_bi .- p1)

    return (L_OD - L_OC)/L_CE
end

function plasma_geometry(bdry::Boundary)

    outer,top,inner,bottom = boundary_extrema(bdry)

    rmax, z_rmax = outer
    r_zmax, zmax = top
    rmin, z_rmin = inner
    r_zmin, zmin = bottom

    r = (rmax - rmin)/2
    κ = (zmax - zmin)/(2*r)
    R0 = (rmax + rmin)/2
    Z0 = (zmax + zmin)/2

    δᵤ = (R0 - r_zmax)/r
    δₗ = (R0 - r_zmin)/r

    κᵤ = (zmax - z_rmax)/r
    κₗ = (z_rmax - zmin)/r

    # Quadrant I
    p1 = r_zmax, z_rmax
    p2 = rmax, zmax
    A, B = p2 .- p1
    ζ_I = quadrant_squareness(bdry,A,B,p1,p2,1)

    # Quadrant II
    p1 = r_zmax, z_rmin
    p2 = rmin, zmax
    A, B = p2 .- p1
    ζ_II = quadrant_squareness(bdry,A,B,p1,p2,2)

    # Quadrant III
    p1 = r_zmin, z_rmin
    p2 = rmin, zmin
    A, B = p2 .- p1
    ζ_III = quadrant_squareness(bdry,A,B,p1,p2,3)

    # Quadrant IV
    p1 = r_zmin, z_rmax
    p2 = rmax, zmin
    A, B = p2 .- p1
    ζ_IV = quadrant_squareness(bdry,A,B,p1,p2,4)

    return PlasmaGeometricParameters(R0, Z0, r, z_rmax, (κₗ,κᵤ), (δₗ, δᵤ), (ζ_I,ζ_II,ζ_III,ζ_IV))
end

function plasma_geometry(M::AbstractEquilibrium)
    return plasma_geometry(shape(M))
end

function flux_surface_geometry(bdry::Boundary)
    return plasma_geometry(bdry)
end
