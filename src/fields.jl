function fields(M::T, r, z) where {T<:AbstractEquilibrium}
    B = Bfield(M, r, z)
    E = Efield(M, r, z)
    return B, E
end

function fields(M::AbstractEquilibrium, x, y, z)
    r = hypot(x, y)
    phi = atan(y, x)
    B, E = fields(M, r, z)
    sp, cp = sincos(phi)
    B = SVector{3}(B[1] * cp - B[2] * sp, B[1] * sp + B[2] * cp, B[3])
    E = SVector{3}(E[1] * cp - E[2] * sp, E[1] * sp + E[2] * cp, E[3])
    return B, E
end

function Bfield(M::T, r, z) where {T<:AbstractEquilibrium}
    psi = M(r, z)
    gval = poloidal_current(M, psi)
    grad_psi = psi_gradient(M, r, z)

    # Calculate B-Field
    cc = cocos(M)
    invr = inv(r)
    #invr not part of cocos factor but added here for performance
    cocos_factor = invr * cc.sigma_RpZ * cc.sigma_Bp / ((2pi)^cc.exp_Bp)

    BR = cocos_factor * grad_psi[2]
    Bz = -cocos_factor * grad_psi[1]
    Bt = gval * invr

    B = SVector{3}(cylindrical_cocos(cc, BR, Bt, Bz))
    return B
end

function Bfield(M::T, x, y, z) where {T<:AbstractEquilibrium}
    r = hypot(x, y)
    phi = atan(y, x)
    B = Bfield(M, r, z)
    sp, cp = sincos(phi)
    B_xyz = SVector{3}(B[1] * cp - B[2] * sp, B[1] * sp + B[2] * cp, B[3])
    return B_xyz
end

function poloidal_Bfield(M::T, r, z) where {T<:AbstractEquilibrium}
    rmaxis, zmaxis = magnetic_axis(M)
    cc = cocos(M)
    B = Bfield(M, r, z)

    ir, iphi, iz = cylindrical_cocos_indices(cc)
    sign_theta = cc.sigma_RpZ * cc.sigma_rhotp # + CW, - CCW
    sign_Bp = sign_theta * sign((z - zmaxis) * B[ir] - (r - rmaxis) * B[iz]) # sign(theta)*sign(r x B)
    Bpol = sign_Bp * sqrt(B[ir]^2 + B[iz]^2)
    return Bpol
end

function Jfield(M::T, r, z; curl=false) where {T<:AbstractEquilibrium}
    if curl
        return curlB(M, r, z) / mu0
    end

    psi = M(r, z)
    gval = poloidal_current(M, psi)
    grad_psi = psi_gradient(M, r, z)

    gp = poloidal_current_gradient(M, psi)
    pp = pressure_gradient(M, psi)

    cc = cocos(M)

    invrmu0 = inv(r * mu0)
    cocos_factor = cc.sigma_RpZ * invrmu0
    jr = -cocos_factor * gp * grad_psi[2]
    jz = cocos_factor * gp * grad_psi[1]

    cocos_factor = -cc.sigma_Bp * ((2pi)^cc.exp_Bp)
    jt = cocos_factor * (r * pp + gval * gp * invrmu0)

    return SVector{3}(cylindrical_cocos(cc, jr, jt, jz))
end

function Jfield(M::T, x, y, z; kwargs...) where {T<:AbstractEquilibrium}
    r = hypot(x, y)
    phi = atan(y, x)
    J = Jfield(M, r, z; kwargs...)
    sp, cp = sincos(phi)
    J_xyz = SVector{3}(J[1] * cp - J[2] * sp, J[1] * sp + J[2] * cp, J[3])
    return J_xyz
end

function poloidal_Jfield(M::T, r, z; kwargs...) where {T<:AbstractEquilibrium}
    rmaxis, zmaxis = magnetic_axis(M)
    cc = cocos(M)
    J = Jfield(M, r, z; kwargs...)

    ir, iphi, iz = cylindrical_cocos_indices(cc)
    sign_theta = cc.sigma_RpZ * cc.sigma_rhotp # + CW, - CCW
    sign_Jp = sign_theta * sign((z - zmaxis) * J[ir] - (r - rmaxis) * J[iz]) # sign(theta)*sign(r x B)
    Jpol = sign_Jp * sqrt(J[ir]^2 + J[iz]^2)
    return Jpol
end

function Efield(M::T, r, z) where {T<:AbstractEquilibrium}

    psi = M(r, z)
    phi_grad = electric_potential_gradient(M, psi)
    if phi_grad == 0
        return SVector{3}(zero(r), zero(r), zero(r))
    end

    Bpol = poloidal_Bfield(M, r, z) #signed

    grad_psi = psi_gradient(M, r, z)
    grad_psi = grad_psi / norm(grad_psi)

    # is Bpol signed in this equation?
    Er = -r * Bpol * phi_grad

    ER = Er * grad_psi[1] # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = Er * grad_psi[2] # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = zero(Ez)

    return SVector{3}(cylindrical_cocos(cocos(M), ER, Et, Ez))
end

function Efield(M::T, x, y, z) where {T<:AbstractEquilibrium}
    r = hypot(x, y)
    phi = atan(y, x)
    E = Efield(M, r, z)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1] * cp - E[2] * sp, E[1] * sp + E[2] * cp, E[3])
    return E_xyz
end

function Efield(M::T, r, z, vrot::AbstractVector) where {T<:AbstractEquilibrium}
    cc = cocos(M)
    psi = M(r, z)
    gval = poloidal_current(M, psi)
    B = Bfield(M, r, z)

    grad_psi = psi_gradient(M, r, z)
    grad_psi = grad_psi / norm(grad_psi)

    qval = safety_factor(M, psi)
    dp = pressure_gradient(M, psi)

    gradp = SVector{3}(cylindrical_cocos(cc, dp * grad_psi[1], zero(grad_psi[1]), dp * grad_psi[2]))

    E = cross(vrot, B) .+ gradp / qval

    return E
end

function Efield(M::T, x, y, z, vrot::AbstractVector) where {T<:AbstractEquilibrium}
    r = hypot(x, y)
    phi = atan(y, x)
    E = Efield(M, r, z, vrot)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1] * cp - E[2] * sp, E[1] * sp + E[2] * cp, E[3])
    return E_xyz
end

function psi_range(M::T; N=100, offset=0.01) where {T<:AbstractEquilibrium}
    psimag, psibdry = psi_limits(M)

    if psibdry > 0 && psimag < 0
        psibdry *= (1 - offset)
    elseif psibdry > 0 && psimag > 0
        psibdry *= (1 + offset)
    elseif psibdry < 0 && psimag < 0
        psibdry *= (1 + offset)
    elseif psibdry < 0 && psimag > 0
        psibdry *= (1 - offset)
    end
    psimag *= (1 - offset)

    return range(psimag, psibdry; length=N)
end

function rho_p(M::T, psi) where {T<:AbstractEquilibrium}
    psimag, psibry = psi_limits(M)
    flux = abs(psibry - psimag)
    return sqrt((psi - psimag) / flux)
end

function rho_p(M::T, r, z) where {T<:AbstractEquilibrium}
    return rho_p(M, M(r, z))
end

function rho_p(M::T, x, y, z) where {T<:AbstractEquilibrium}
    r = hypot(x, y)
    return rho_p(M, r, z)
end

function toroidal_flux(M::T, psi; N=100) where {T<:AbstractEquilibrium}
    cc = cocos(M)
    cocos_factor = cc.sigma_Bp * cc.sigma_rhotp * (2pi)^(1 - cc.exp_Bp)

    psimag, psibry = psi_limits(M)
    psi_arr = range(psimag * 0.99, psi; length=N)
    q = safety_factor.(M, psi_arr)
    phi = cocos_factor * trapz(psi_arr, q)
    return phi
end

function safety_factor(M::AbstractEquilibrium, psi)
    g = poloidal_current(M, psi)
    fs = flux_surface(M, psi)
    q = (g / (2pi)) * average(fs, (x, y) -> inv(poloidal_Bfield(M, x, y) * x^2)) * circumference(fs)
    return q
end

function plasma_current(M::AbstractEquilibrium)
    cc = cocos(M)
    S = shape(M)
    Ip = cc.sigma_rhotp * integrate(x -> poloidal_Bfield(M, S(x)...), S, :line, (0.0, 2pi))
    return Ip / mu0
end

function beta_generic(M::T) where {T<:AbstractEquilibrium}
    S = shape(M)
    p_bar = average(x -> pressure(M, M(S(x[1], x[2])...)), S, :volume)
    B_bar = average(x -> norm(Bfield(M, S(x[1], x[2])...)), S, :volume)
    return 100 * p_bar / (B_bar^2 / (2 * mu0))
end
beta(M::AbstractEquilibrium) = beta_generic(M)

function beta_p_generic(M::T) where {T<:AbstractEquilibrium}
    S = shape(M)
    p_bar = average(x -> pressure(M, M(S(x[1], x[2])...)), S, :volume)
    Bp_bar = average(x -> poloidal_Bfield(M, S(x[1], x[2])...), S, :volume)
    return 100 * p_bar / (Bp_bar^2 / (2 * mu0))
end
beta_p(M::AbstractEquilibrium) = beta_p_generic(M)

function beta_t_generic(M::T) where {T<:AbstractEquilibrium}
    S = shape(M)
    p_bar = average(x -> pressure(M, M(S(x[1], x[2])...)), S, :volume)
    Bt_bar = average(x -> Bfield(M, S(x[1], x[2])...)[2], S, :volume)
    return 100 * p_bar / (Bt_bar^2 / (2 * mu0))
end
beta_t(M::AbstractEquilibrium) = beta_t_generic(M)

function beta_n_generic(M::T) where {T<:AbstractEquilibrium}
    Ip = plasma_current(M) * 1e-6
    β = beta(M)
    a = minor_radius(shape(M))
    Bt = Bfield(M, magnetic_axis(M)...)[2]
    βn = β * (a * Bt) / Ip
    return βn
end
beta_n(M::AbstractEquilibrium) = beta_n_generic(M)

function gradB_autodiff(M::AbstractEquilibrium, r, z)
    gB_rz = ForwardDiff.gradient(x -> norm(Bfield(M, x[1], x[2])), SVector{2}(r, z))
    return SVector{3}(cylindrical_cocos(cocos(M), gB_rz[1], 0.0, gB_rz[2]))
end

function gradB_autodiff(M::AbstractEquilibrium, x, y, z)
    gB_xyz = ForwardDiff.gradient(x -> norm(Bfield(M, x[1], x[2], x[3])), SVector{3}(x, y, z))
    return SVector{3}(gB_xyz[1], gB_xyz[2], gB_xyz[3])
end

gradB(M::AbstractEquilibrium, r, z) = gradB_autodiff(M, r, z)

function gradB(M::AbstractEquilibrium, x, y, z)
    r = hypot(x, y)
    phi = atan(y, x)
    gB = gradB(M, r, z)
    sp, cp = sincos(phi)
    return gB_xyz = SVector{3}(gB[1] * cp - gB[2] * sp, gB[1] * sp + gB[2] * cp, gB[3])
end

function curlB_autodiff(M::AbstractEquilibrium, r, z)
    cc = cocos(M)
    B = Bfield(M, r, z)

    # Complete version
    #J = ForwardDiff.jacobian(x->Bfield(M,x[1],x[3]),SVector{3}(r,zero(r),z))
    #return cc.sigma_RpZ*SVector{3}(cylindrical_cocos(cc, J[3,2]/r - J[2,3], J[1,3] - J[3,1], (B[2]/r + J[2,1]) - J[1,2]/r))

    # Exploit axisymmetry to simplify
    J = ForwardDiff.jacobian(x -> Bfield(M, x[1], x[2]), SVector{2}(r, z))
    return cc.sigma_RpZ * SVector{3}(cylindrical_cocos(cc, -J[2, 2], J[1, 2] - J[3, 1], (B[2] / r + J[2, 1])))
end

function curlB_autodiff(M::AbstractEquilibrium, x, y, z)
    J = ForwardDiff.jacobian(x -> Bfield(M, x[1], x[2], x[3]), SVector{3}(x, y, z))
    return SVector{3}(J[3, 2] - J[2, 3], J[1, 3] - J[3, 1], J[2, 1] - J[1, 2])
end

curlB(M::AbstractEquilibrium, r, z) = curlB_autodiff(M, r, z)

function curlB(M::AbstractEquilibrium, x, y, z)
    r = hypot(x, y)
    phi = atan(y, x)
    cB = curlB(M, r, z)
    sp, cp = sincos(phi)
    cB_xyz = SVector{3}(cB[1] * cp - cB[2] * sp, cB[1] * sp + cB[2] * cp, cB[3])
    return cB_xyz
end
