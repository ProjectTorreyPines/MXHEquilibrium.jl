function _reorder!(r, z, R0, Z0)
    # largest r above Z0 on low field side
    j = argmax((r .- R0) .+ -Inf*(z .< Z0))
    r .= circshift(r,1-max(j))
    z .= circshift(z,1-max(j))
    if z[2] - z[1] < 0
        reverse!(z)
        reverse!(r)
    end
end

function MXH_angles(bdry::Boundary; n_interp=0, debug=false)

    R0, Z0, r, κ = plasma_geometry(bdry)
    
    r_s = copy(bdry.r)
    z_s = copy(bdry.z)

    _reorder!(r_s,z_s,R0,Z0) # put start in first quadrant
    
    θ = zero(r_s)
    Δθᵣ = zero(r_s)
    th = 0.0
    thr = 0.0
    quad = 1
    @inbounds for j in eachindex(θ)
        th_old = th
        thr_old = thr

        aa = clamp((z_s[j] - Z0)/(κ*r) ,-1,1)
        th = asin(aa) 

        bb = clamp((r_s[j] - R0) / r,-1,1)
        thr = acos(bb)

        # asin(-1) = -pi/2, asin(1) = pi/2
        # acos(-1) = pi, acos(1) = 0
        if (j == 1) || ((th > th_old) && (thr > thr_old)) #first quadrant
            debug && println("index = $j, quad = 1: sin($(round(aa,digits=3))) = \
                             $(round(th,digits=3)), cos($(round(bb,digits=3))) = \
                             $(round(thr,digits=3)), r-R0 = $(r_s[j] - R0), z-Z0= $(z_s[j] - Z0)")
            θ[j] = th
            Δθᵣ[j] = thr
        elseif (th < th_old) && (thr > thr_old) # second quadrant
            debug && println("index = $j, quad = 2: sin($(round(aa,digits=3))) = \
                             $(round(th,digits=3)), cos($(round(bb,digits=3))) = \
                             $(round(thr,digits=3)), r-R0 = $(r_s[j] - R0), z-Z0= $(z_s[j] - Z0)")
            θ[j] = π - th
            Δθᵣ[j] = thr
        elseif (th < th_old) && (thr < thr_old) # third quadrant
            debug && println("index = $j, quad = 3: sin($(round(aa,digits=3))) = \
                             $(round(th,digits=3)), cos($(round(bb,digits=3))) = \
                             $(round(thr,digits=3)), r-R0 = $(r_s[j] - R0), z-Z0= $(z_s[j] - Z0)")
            θ[j] = π - th
            Δθᵣ[j] = 2pi - thr
        elseif (th > th_old) && (thr < thr_old) # fourth quadrant
            debug && println("index = $j, quad = 4: sin($(round(aa,digits=3))) = \
                             $(round(th,digits=3)), cos($(round(bb,digits=3))) = \
                             $(round(thr,digits=3)), r-R0 = $(r_s[j] - R0), z-Z0= $(z_s[j] - Z0)")
            θ[j] = 2π + th
            Δθᵣ[j] = 2π - thr
        end
        Δθᵣ[j] -= θ[j]
    end
    w = sortperm(θ)
    θ .= θ[w]
    Δθᵣ .= Δθᵣ[w]

    if n_interp > 0
        Interpolations.deduplicate_knots!(θ)
        itp = linear_interpolation((θ,),Δθᵣ,extrapolation_bc=Periodic())
        τ = range(0,2pi,length=n_interp)
        return τ, itp.(τ)
    end
    return θ, Δθᵣ
end

function MXH_moment(x, F, w)
    return trapz(x,F.*w)/trapz(x,w.^2)
end

function MXH_coeffs(θ::AbstractVector{T},θ_r::AbstractVector{T};N=7) where T
    cs = zeros(T,N)
    ss = zeros(T,N)

    c0 = MXH_moment(θ, θ_r, ones(length(θ_r)))
    for n = 1:N
        Fn = cos.(n .* θ)
        cs[n] = MXH_moment(θ, θ_r, Fn)
    end
    
    for n=1:N
        Fn = sin.(n .* θ)
        ss[n] = MXH_moment(θ, θ_r, Fn)
    end
    return c0, cs, ss
end

function MXH_parameters(bdry::Boundary; N=7, kwargs...)
    
    R0, Z0, r, κ, = plasma_geometry(bdry)
    
    θ, θr = MXH_angles(bdry; kwargs...)
    
    c0, c, s = MXH_coeffs(θ, θr; N=N)
    
    return R0, Z0, r/R0, κ, c0, c, s
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
    R0, Z0, r, κ, δl, δu = plasma_geometry(bdry)
    return MillerShape(R0, Z0, r/R0, κ, (δl + δu)/2)
end

function fit(bdry::Boundary, S::AsymmetricMillerShape; kwargs...)
    if length(bdry) == 2
        s = AMShape(bdry.r[1], bdry.z[1], 0.0, 0.0, 0.0, 0.0)
        return s
    end
    R0, Z0, r, κ, δl, δu = plasma_geometry(bdry)
    return AsymmetricMillerShape(R0, Z0, r/R0, κ, δl, δu)
end
