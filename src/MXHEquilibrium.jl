module MXHEquilibrium

using EFIT
using LinearAlgebra
using Interpolations
using ForwardDiff
using Zygote
using StaticArrays
using HCubature
using Trapz

using Memoize
using LRUCache

using MeshTools

import PolygonOps
import Contour
import Optim

using CoordinateConventions
import CoordinateConventions: cocos, identify_cocos, check_cocos, transform_cocos
import MillerExtendedHarmonic: MXH

const mu0 = 4 * pi * 1e-7

cache_size = 5

function set_cache_size!(N::Int)
    m = @__MODULE__
    cache = LRU(; maxsize=cache_size)
    global cache_size = N
    for f in names(m)
        cachename = Memoize.cache_name(f)
        try
            cache = getproperty(m, cachename)
        catch e
        finally
            resize!(cache; maxsize=cache_size)
        end
    end
end

abstract type AbstractEquilibrium end
export AbstractEquilibrium

# Equilibrium Fallbacks
_not_implemented(M) = error("$(typeof(M)) has not implemented this functionality")
(M::AbstractEquilibrium)(x, y) = _not_implemented(M)
magnetic_axis(M::AbstractEquilibrium, r, z) = _not_implemented(M)
limits(M::AbstractEquilibrium) = _not_implemented(M)
psi_limits(M::AbstractEquilibrium) = _not_implemented(M)
psi_gradient(M::AbstractEquilibrium) = _not_implemented(M)
electric_potential(M::AbstractEquilibrium, psi) = _not_implemented(M)
electric_potential_gradient(M::AbstractEquilibrium) = _not_implemented(M)
pressure(M::AbstractEquilibrium, psi) = _not_implemented(M)
poloidal_current(M::AbstractEquilibrium, psi) = _not_implemented(M)
pressure_gradient(M::AbstractEquilibrium, psi) = _not_implemented(M)
poloidal_current_gradient(M::AbstractEquilibrium, psi) = _not_implemented(M)
cocos(M::AbstractEquilibrium) = _not_implemented(M)
B0Ip_sign(M::AbstractEquilibrium) = _not_implemented(M)
psi_boundary(M::AbstractEquilibrium) = _not_implemented(M)
plasma_boundary_psi(M::AbstractEquilibrium) = (plasma_boundary(M), psi_boundary(M))

# Equilibrium API
export magnetic_axis, limits, psi_limits, psi_boundary, plasma_boundary_psi, psi_gradient, electric_potential, electric_potential_gradient
export pressure, poloidal_current, pressure_gradient, poloidal_current_gradient, safety_factor, B0Ip_sign
export plasma_current, beta, beta_n, beta_p, beta_t, psi_range, rho_p, toroidal_flux

include("cocos.jl")
export COCOS, cocos, check_cocos, identify_cocos, transform_cocos
export cylindrical_cocos, poloidal_cocos, cylindrical_cocos_indices, poloidal_cocos_indices

include("boundary.jl")
export Boundary, PlasmaBoundary, FluxSurface, Wall, in_boundary, in_plasma, in_vessel
export boundary, flux_surface, plasma_boundary, boundary_extrema
export circumference, average, area, area_average, volume, volume_average, surface_area

include("geometry.jl")
export PlasmaGeometricParameters, plasma_geometry, flux_surface_geometry

include("shape.jl")
export PlasmaShape, MillerShape, TurnbullMillerShape, MillerExtendedHarmonicShape, LuceShape
export AsymmetricMillerShape, AMShape, MShape, TMShape, MXHShape, LShape, shape
export curvature, triangularity, squareness, tilt, elevation, ovality, twist
export scale_aspect, elongation, aspect_ratio, major_radius, minor_radius
export convert_eltype

# Shape Fallbacks
(S::PlasmaShape)(x) = _not_implemented(S)
shape(S::PlasmaShape) = _not_implemented(S)
aspect_ratio(S::PlasmaShape) = _not_implemented(S)
elongation(S::PlasmaShape) = _not_implemented(S)
major_radius(S::PlasmaShape) = _not_implemented(S)
minor_radius(S::PlasmaShape) = _not_implemented(S)
elevation(S::PlasmaShape) = _not_implemented(S)
triangularity(S::PlasmaShape) = _not_implemented(S)
tilt(S::PlasmaShape) = _not_implemented(S)
ovality(S::PlasmaShape) = _not_implemented(S)
twist(S::PlasmaShape) = _not_implemented(S)
squareness(S::PlasmaShape) = _not_implemented(S)

include("fitting.jl")
export fit

include("solovev.jl")
export SolovevEquilibrium, solovev, clear_cache

include("fields.jl")
export Bfield, Efield, Jfield, EMFields, fields, gradB, curlB, poloidal_Bfield, poloidal_Jfield

include("efit.jl")
export EFITEquilibrium, efit

include("efit_io.jl")
export read_geqdsk

include("transp_io.jl")
export transp_potential!

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module
