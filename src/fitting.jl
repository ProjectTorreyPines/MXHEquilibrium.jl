function MXH_parameters(bdry::Boundary; N=7, kwargs...)

    mxh = MXH(bdry.r, bdry.z, N; kwargs...)
    R0 = mxh.R0
    Z0 = mxh.Z0
    κ = mxh.κ
    ϵ = mxh.ϵ
    c0 = mxh.c0
    c = mxh.c
    s = mxh.s

    return R0, Z0, ϵ, κ, c0, c, s
end

function fit(bdry::Boundary, S::MillerExtendedHarmonicShape{N,T}; kwargs...) where {N,T}
    if length(bdry) == 2
        s = MXHShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0, zeros(SVector{N}), zeros(SVector{N}))
        return s
    end
    return MXHShape(MXH_parameters(bdry; N=N, kwargs...)...)
end

fit(bdry::Boundary) = fit(bdry,MillerShape())

function fit(bdry::Boundary, S::MillerShape; kwargs...)
    if length(bdry) == 2
        s = MShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0)
        return s
    end

    G = plasma_geometry(bdry)
    return MillerShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2)
end

function fit(bdry::Boundary, S::AsymmetricMillerShape; kwargs...)
    if length(bdry) == 2
        s = AMShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry)
    return AMShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, G.δ...)
end

function fit(bdry::Boundary, S::TurnbullMillerShape; kwargs...)
    if length(bdry) == 2
        s = TMShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry)
    return TMShape(G.R0, G.Z0, G.r/G.R0, sum(G.κ)/2, sum(G.δ)/2, sum(G.ζ)/4)
end

function fit(bdry::Boundary, S::LShape; kwargs...)
    if length(bdry) == 2
        s = TMShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0, 0.0)
        return s
    end
    G = plasma_geometry(bdry)
    return LuceShape(G)
end

@memoize LRU(maxsize=cache_size) function shape(M::AbstractEquilibrium; N=10)
    bdry = plasma_boundary(M)
    return fit(bdry,MXHShape(N))
end
