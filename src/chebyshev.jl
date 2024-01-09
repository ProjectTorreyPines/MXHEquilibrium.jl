using Polynomials

struct ChebyshevEquilibrium{T,F} <: AbstractEquilibrium
    cocos::COCOS
    psi_coeffs::Vector{T}
    nx::Int
    ny::Int
    g::F
    rlims::NTuple{2,T}
    zlims::NTuple{2,T}
    psi_lims::NTuple{2,T}
end

chebT(x,n) = abs(x) <= 1.0 ? cos(n*acos(x)) : error("x out of range")
function chebT(x,n::Union{Range,Vector{Int}})
    abs(x) > 1.0 && error("x out of range")
    acx = acos(x)
    return [cos(nn*acx) for nn in n]
end

function _chebyshev_coeffs(x, y, psi, nx, ny)
    psi_coeffs = transpose(hcat((vec([chebT(xx,nn)*chebT(yy,mm)
                                      for mm=0:(ny-1), nn=0:(nx-1)]) 
                                 for yy in y, xx in x)...)) \ vec(psi);
    return psi_coeffs
end

function chebyshev(r, z, psi, g::Function, psi_lims; nx=10, ny=20, cc=3)
    rmin, rmax = extrema(r)
    x = (2/(rmax - rmin))*(r .- rmin) .- 1
    zmin, zmax = extrema(z)
    y = (2/(zmax - zmin))*(z .- zmin) .- 1
    psi_coeffs = _chebyshev_coeffs(x,y,psi,nx,ny)

    ChebyshevEquilibrium(cocos(cc),psi_coeffs,nx,ny,g,(rmin,rmax),(zmin,zmax),psi_lims)
end

function (CE::ChebyshevEquilibrium)(r,z)
    rmin,rmax = CE.rlims
    zmin,zmax = CE.zlims
    x = (2/(rmax - rmin))*(r - rmin) - 1
    y = (2/(zmax - zmin))*(z - zmin) - 1
    return vec([chebT(x,nn)*chebT(y,mm) for mm=0:(CE.ny-1),nn=0:(CE.nx-1)])'*CE.psi_coeffs
end

cocos(CE::ChebyshevEquilibrium) = CE.cocos

function limits(CE::ChebyshevEquilibrium)
    return CE.rlims, CE.zlims
end

function psi_limits(CE::ChebyshevEquilibrium)
    return CE.psi_lims
end

function poloidal_current(CE::ChebyshevEquilibrium, psi)
    return CE.g(psi)
end
