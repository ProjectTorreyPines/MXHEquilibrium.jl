function _reorder!(r, z, R0, Z0)
    # first z above Z0 on low field side
    i = argmin(abs.(z .- Z0) .+ (r .< R0) .+ (z .< Z0))
    circshift!(r,1-i)
    circshift!(z,1-i)
end

function MXH_angles(bdry::Boundary)

    R0, Z0, r, κ = plasma_geometry(bdry)
    
    r_s = copy(bdry.r)
    z_s = copy(bdry.z)

    _reorder!(r_s,z_s,R0,Z0)
    
    θ = zero(r_s)
    Δθᵣ = zero(r_s)
    th = 0.0
    thr = 0.0
    @inbounds for j in eachindex(θ)
        th_old = th
        thr_old = thr

        aa = clamp((z_s[j] - Z0)/(κ*r) ,-1,1)
        th = asin(aa)
        
        bb = clamp((r_s[j] - R0) / r,-1,1)
        thr = acos(bb)
        
        if (j == 1) || ((th > th_old) && (thr > thr_old))
            θ[j] = th
            Δθᵣ[j] = thr
        elseif (th < th_old) && (thr > thr_old)
            θ[j] = π - th
            Δθᵣ[j] = thr
        elseif (th < th_old) && (thr < thr_old)
            θ[j] = π - th
            Δθᵣ[j] = 2π - thr
        elseif (th > th_old) && (thr < thr_old)
            θ[j] = 2π + th
            Δθᵣ[j] = 2π - thr
        end
        Δθᵣ[j] -= θ[j]
    end
    w = sortperm(θ)
    return θ[w], Δθᵣ[w]
end

function MXH_moment(x, F, w)
    return trapz(x,F.*w)/trapz(x,w.^2)
end

function MXH_coeffs(θ,θ_r;N=7)
    cs = zeros(N)
    ss = zeros(N)
    
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

function MXH_parameters(bdry::Boundary; N=7)
    
    R0, Z0, r, κ, = plasma_geometry(bdry)
    
    θ, θr = MXH_angles(bdry)
    
    c0, c, s = MXH_coeffs(θ, θr; N=N)
    
    return R0, Z0, r/R0, κ, c0, c, s
end

function fit(bdry::Boundary, S::MillerExtendedHarmonicShape{N,T}) where {N,T}
    return MXHShape(MXH_parameters(bdry,N=N)...)
end

fit(bdry::Boundary) = fit(bdry,MXHShape(10))

function fit(bdry::Boundary, S::MillerShape)
    R0, Z0, ϵ, κ, δl, δu = plasma_geometry(bdry)
    return MillerShape(R0, Z0, ϵ, κ, (δl + δu)/2)
end

function fit(bdry::Boundary, S::AsymmetricMillerShape)
    R0, Z0, ϵ, κ, δl, δu = plasma_geometry(bdry)
    return AsymmetricMillerShape(R0, Z0, ϵ, κ, δl, δu)
end
