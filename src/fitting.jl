function MXH_parameters(bdry::Boundary; N=7, kwargs...)

    mxh = MXH(bdry.r, bdry.z, N; kwargs...)
    R0 = mxh.R0
    Z0 = mxh.Z0
    κ = mxh.κ
    ϵ = mxh.ϵ
    c0 = -mxh.c0
    c = -mxh.c
    s = mxh.s

    return R0, Z0, ϵ, κ, c0, c, s
end

function fit(bdry::Boundary, S::MillerExtendedHarmonicShape{N,T}; normalize=false, optimize_ncoeffs=false, kwargs...) where {N,T}
    if length(bdry) == 2
        if normalize
            R0 = one(bdry.r[1])
            Z0 = zero(bdry.z[1])
        else
            R0 = bdry.r[1]
            Z0 = bdry.z[1]
        end
        s = MXHShape(R0, Z0, 0.0, 0.0, 0.0, zeros(SVector{N}), zeros(SVector{N}))
        return s
    end
    R0, Z0, ϵ, κ, c0, c, s = MXH_parameters(bdry; N=N, kwargs...)

    mxh = MXHShape(R0,Z0,ϵ,κ,c0,c,s)
    if optimize_ncoeffs
        imax = 0
        omax = 0.0
        for i=1:N
            bdry_mxh = Boundary(shape(mxh, ncoeffs=i)...)
            o = overlap_metric(bdry, bdry_mxh)
            if o > omax
                imax = i
                omax = o
            end
        end
        c[imax+1:end] .= 0.0
        s[imax+1:end] .= 0.0
    end

    if normalize
        R0 = one(R0)
        Z0 = zero(Z0)
    end

    return MXHShape(R0,Z0,ϵ,κ,c0,c,s)
end

fit(bdry::Boundary) = fit(bdry,MillerShape())

function fit(bdry::Boundary, S::MillerShape; normalize = false, kwargs...)
    if length(bdry) == 2
        if normalize
            R0 = one(bdry.r[1])
            Z0 = zero(bdry.z[1])
        else
            R0 = bdry.r[1]
            Z0 = bdry.z[1]
        end
        s = MShape(R0, Z0, 0.0, 0.0, 0.0)
        return s
    end

    G = plasma_geometry(bdry)
    if normalize
        return MillerShape(one(G.R0), zero(G.Z0), G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2)
    else
        return MillerShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2)
    end
end

function fit(bdry::Boundary, S::AsymmetricMillerShape; normalize = false, kwargs...)
    if length(bdry) == 2
        if normalize
            R0 = one(bdry.r[1])
            Z0 = zero(bdry.z[1])
        else
            R0 = bdry.r[1]
            Z0 = bdry.z[1]
        end
        s = AMShape(R0, Z0, 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry)
    if normalize
        return AMShape(zero(G.R0), one(G.Z0), G.r/G.R0, sum(G.κ)/2, G.δ...)
    else
        return AMShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, G.δ...)
    end
end

function fit(bdry::Boundary, S::TurnbullMillerShape; normalize=false,kwargs...)
    if length(bdry) == 2
        if normalize
            R0 = one(bdry.r[1])
            Z0 = zero(bdry.z[1])
        else
            R0 = bdry.r[1]
            Z0 = bdry.z[1]
        end
        s = TMShape(R0, Z0, 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry)
    if normalize
        return TMShape(zero(G.R0), one(G.Z0), G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2, sum(G.ζ)/4)
    else
        return TMShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2, sum(G.ζ)/4)
    end
end

function fit(bdry::Boundary, S::LShape; normalize=false, kwargs...)
    if length(bdry) == 2
        if normalize
            R0 = one(bdry.r[1])
            Z0 = zero(bdry.z[1])
        else
            R0 = bdry.r[1]
            Z0 = bdry.z[1]
        end
        s = TMShape(R0, Z0, 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry; normalize=normalize)
    return LuceShape(G)
end

@memoize LRU(maxsize=cache_size) function shape(M::AbstractEquilibrium; N=10)
    bdry = plasma_boundary(M)
    return fit(bdry,MXHShape(N))
end
