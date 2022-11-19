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

    G = plasma_geometry(bdry)
    R0 = G.R0
    Z0 = G.Z0
    r = G.r
    κ = sum(G.κ...)/2

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
    
    G = plasma_geometry(bdry)
    R0 = G.R0
    Z0 = G.Z0
    r = G.r
    κ = sum(G.κ...)/2
    
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

function x_point_constraints(G::PlasmaGeometricParameters, quad::Int, x::Vector)
    α, β, n = x
    a = G.r
    κl, κu = G.κ
    δl, δu = G.δ
    ζ_I,ζ_II,ζ_III,ζ_IV = G.ζ
    ir2 = inv(sqrt(2))
    #x-point constraint
    if quad == 1
        C1 = abs(1 - (a/α)*(1 + δu))^n + abs(κu*a/β)^n - 1
    
        C2 = abs(1 + a*(1 + δu)*(1 - ir2)*(ζ_I - 1)/α)^n + 
             abs(κu*a*(ir2 + ζ_I*(1 - ir2))/β)^n - 1

        C3 = β*abs((α - a*(1 + δu))/(κu*a))^(n-1) * abs(β/α)^(n-1) -
             (δu >= 0.0 ? (1 - δu) : inv(1 + δu))

    elseif quad == 2
        C1 = abs(1 - a/α*(1 - δu))^n + abs(κu*a/β)^n - 1

        C2 = abs(1 + a*(1 - δu)*(1 - ir2)*(ζ_II - 1)/α)^n + 
             abs(κu*a*(ir2 + ζ_II*(1 - ir2))/β)^n - 1

        C3 = β*abs((α - a*(1 - δu))/κu*a)^(n-1) * abs(β/α)^(n-1) -
             (δu >= 0.0 ? inv(1 + δu) : (1 - δu))

    elseif quad == 3
        C1 = abs((1 - a*(1 - δl))/α)^n + abs(κl*a/β)^n - 1

        C2 = abs((1 + a*(1 - δl)*(1 - ir2)*(ζ_III - 1))/α)^n + 
             abs(κl*a*(ir2 + ζ_III*(1 - ir2))/β)^n - 1

        C3 = β*abs((α - α*(1 - δl))/(κl*a))^(n-1) * abs(β/a)^(n-1) -
             (δl >= 0.0 ? inv(1 + δl) : (1 - δl))
    else
        C1 = abs((1 - a*(1 + δl))/α)^n + abs(κl*a/β)^n - 1

        C2 = abs((1 + a*(1 + δl)*(1 - ir2)*(ζ_IV - 1))/α)^n + 
             abs(κl*a*(ir2 + ζ_IV*(1 - ir2))/β)^n - 1

        C3 = β*abs((α - α*(1 - δl))/(κl*a))^(n-1) * abs(β/α)^(n-1) -
             (δl >= 0.0 ? (1 + δl) : inv(1 - δl))


    end

    return (C1 + C2 + C3)^2
end

function fit(bdry::Boundary, S::LuceShape; x_point = 0, kwargs...)
    if length(bdry) == 2
        Z = zero(bdry.r[1])
        s = LuceShape(bdry.r[1], bdry.z[1], Z, Z,(Z,Z),(Z,Z),
                      (Z,Z,Z,Z),(Z,Z,Z,Z))
        return s
    end

    G = plasma_geometry(bdry)

    if x_point == 0
        return LuceShape(G)
    else
        quad = x_point
        alphas = zeros(4)
        res = Optim.optimize(x->x_point_constraints(G,quad,x), [G.r,G.r,2.0], Optim.NelderMead())
        α,β,n = res.minimizer
        println("A = ",α)
        println("B = ",β)
        println("n = ",n)
        println("C = ", x_point_constraints(G,quad,[α,β,n]))
        alphas[quad] = α
        alphas = ntuple(i->alphas[i],4)
        return LuceShape(G,alphas)
    end
end
