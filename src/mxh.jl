struct MXHEquilibrium{S<:AbstractVector,Q<:AbstractVector,N} <: AbstractEquilibrium
    r::S                      # Normalized radius
    ψ::Q                      # Poloidal Flux
    surf::Vector{MHXShape{N}} # Vector of plasma shapes
    sigma_B0::Int
    sigma_Ip::Int
    x_point::Union{Nothing,NTuple{2}}
#    cocos::COCOS
end

function MXHEquilibrium(surfs,ψ; B0_dir=1, Ip_dir=1, x_point=nothing)
    ϵ_s = getfield.(surfs,:ϵ)
    R0s = getfield.(surfs,:R0)
    r = ϵ_s .* R0s
    r *= inv(R0[1])
    return MXHEquilibrium(r, ψ, surfs, B0_dir, Ip_dir, x_point)
end

function (M::MXHEquilibrium)(r,θ)
end

function shape(M::MXHEquilibrium)
    return S.surf[end]
end

function Base.show(io::IO, M::MXHEquilibrium{S,Q,N}) where {S,Q,N}
    print(io, "MXHEquilibrium{$(N)\n")
    print(io, "  B0 = $(M.B0) [T]\n")
    print(io, "  βp = $(M.beta_p)\n")
    print(io, "  βt = $(M.beta_t)\n")
    print(io, "  σ_B0 = $(M.sigma_B0)\n")
    print(io, "  σ_Ip = $(M.sigma_Ip)\n")
    print(io, "  x_point = $(M.x_point)\n")
end

Base.broadcastable(M::MXHEquilibrium) = (M,)

cocos(S::MXHEquilibrium) = M.cocos

function B0Ip_sign(M::MXHEquilibrium)
    return M.sigma_B0*M.sigma_Ip
end

function limits(M::MXHEquilibrium)
    xlims, ylims = limits(shape(M),M.x_point)
end

#function psi_gradient(S::MXHEquilibrium,r,θ)
#end

function magnetic_axis(M::MXHEquilibrium)
    return major_radius(M.surf[1]), elevation(M.surf[1])
end

function psi_limits(M::MXHEquilibrium)
    psimag = M.psi[1]
    psibry = M.psi[end]
    return (psimag, psibry)
end

#function pressure_gradient(S::MXHEquilibrium)
#end

function pressure_gradient(M::MXHEquilibrium, psi)
    return pressure_gradient(M) + psi*0
end

#function pressure(M::MXHEquilibrium, psi; p0=zero(psi))
#end

#function poloidal_current(M::MXHEquilibrium,psi)
#end

#function poloidal_current_gradient(M::MXHEquilibrium,psi)
#end

electric_potential(M::MXHEquilibrium,psi) = zero(psi)

electric_potential_gradient(M::MXHEquilibrium,psi) = zero(psi)

#beta_p(M::MXHEquilibrium) = M.beta_p
#beta_t(M::MXHEquilibrium) = M.beta_t
