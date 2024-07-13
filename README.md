# MXHEquilibrium.jl

MXHEquilibrium.jl is a fork of Equilibrium.jl and provides functionality for Miller Extended Harmonic fitting.

## AbstractEquilibrium API

```julia
using MXHEquilibrium
typeof(S) <: AbstractEquilibrium

psi = S(r, z) # Poloidal flux at r,z
gradpsi = psi_gradient(S, r, z)

B = Bfield(S, r, z)
Bp = poloidal_Bfield(S, r, z)

J = Jfield(S, r, z)
Jp = poloidal_Jfield(S, r, z)

F = poloidal_current(S, psi)
Fprime = poloidal_current_gradient(S, psi)

p = pressure(S, psi)
pprime = pressure_gradient(S, psi)

V = electric_potential(S, psi)
gradV = electric_potential_gradient(S, psi)

q = safety_factor(S, psi)

maxis = magnetic_axis(S)

btip = B0Ip_sign(S)

rlims, zlims = limits(S)

psi_lims = psi_limits(S)

cc = cocos(S) # Return COCOS structure

fs = flux_surface(S, psi) # returns a boundary object

```

## Solov'ev Equilibrium
Solov'ev Equilibrium are analytic solutions to the Grad-Shafranov equation where the p' and FF' are constant.
The resulting Grad-Shafranov equation takes the form `Δ⋆ψ = α + (1-α)x²` where `α` is some constant.
The boundary conditions are found using a plasma shape parameterization.

```julia

# ITER parameters
δ = 0.33        # Triangularity
ϵ = 0.32        # Inverse aspect ratio a/R0
κ = 1.7         # Elongation
B0 = 5.3        # Magnitude of Toroidal field at R0 [T]
R0 = 6.2        # Major Radius [m]
Z0 = 0.0        # Elevation
qstar = 1.57    # Kink safety factor
alpha = -0.155  # constant

S = solovev(B0, MillerShape(R0, Z0, ϵ, κ, δ), alpha, qstar, B0_dir=1, Ip_dir=1)


SolovevEquilibrium
  B0 = 2.0 [T]
  S  = MillerShape{Float64}(6.2, 0.0, 0.32, 1.7, 0.33)
  α  = -0.155
  q⋆ = 1.57
  βp = 1.1837605469381924
  βt = 0.049177281028224634
  σ  = 1
  diverted  = false
  symmetric = true
```

## EFIT Equilibrium
EFIT geqdsk files are a commonly used file format.
Here we provide routines for converting the GEQDSK files into an Equilibrium object.

```julia
using EFIT

g = readg("g000001.01000")
M = efit(g, clockwise_phi=false) # direction of phi needed to determine COCOS ID
wall = Wall(g)

in_vessel(wall, r, z)

# or
# M, wall = read_geqdsk("g000001.01000",clockwise_phi=false)

```

## COCOS: Tokamak Coordinate Conventions
We provide routines for working determining, transforming, and checking COCOS's.
```julia
julia> cocos(3)
COCOS = 3
 e_Bp  = 0
 σ_Bp  = -1
 σ_RΦZ = (R,Φ,Z): 1
 σ_ρθΦ = (ρ,Φ,θ): -1
 Φ from top: CCW
 θ from front: CCW
 ψ_ref: Decreasing assuming +Ip, +B0
 sign(q) = -1 assuming +Ip, +B0
 sign(p') = 1 assuming +Ip, +B0

julia> transform_cocos(3,1)
Dict{Any, Any} with 14 entries:
  "Z"        => 1.0
  "Q"        => -1
  "P"        => 1.0
  "B"        => 1.0
  "F_FPRIME" => -1.0
  "ψ"        => -1.0
  "TOR"      => 1.0
  "Φ"        => 1.0
  "PSI"      => -1.0
  "I"        => 1.0
  "J"        => 1.0
  "R"        => 1.0
  "F"        => 1.0
  "PPRIME"   => -1.0
```

## Boundaries
MXHEquilibrium.jl also provides routines for working with boundries such as walls or flux surfaces. Internally boundaries are stored as a list of points forming a polygon.

```julia

fs = flux_surface(S, psi)

in_plasma(fs, r, z) # or in_vessel(fs, r, z), in_boundary(fs, r, z)

cicumference(fs)

area(fs) # Area enclosed by the boundary

volume(fs) # assuming toroidal symmetry. F can be a vector with the same length as fs or a function of (r,z)

average(fs, F) # Average F over the boundary

area_average(fs, F) # average F over the area

volume_average(fs, F) # average F over the volume
```

## Parameterized Plasma Shapes
MXHEquilibrium.jl provides the commonly used plasma shape parameterizations.

```julia
help?> MillerShape # alias = MShape

  MillerShape Structure

  Defines the Miller Plasma Shape Parameterization

  Fields:
  R0 - Major Radius [m]
  Z0 - Elevation [m]
  ϵ - Inverse Aspect Ratio a/R0 where a = minor radius
  κ - Elongation
  δ - Triangularity

help?> TurnbullMillerShape # alias = TMShape

  TurnbullMillerShape Structure

  Defines the Turnbull-Miller Plasma Shape Parameterization
  > Turnbull, A. D., et al. "Improved magnetohydrodynamic stability through optimization of higher order moments in cross-section shape of tokamaks." Physics of Plasmas 6.4 (1999): 1113-1116.
 
  Fields:
  R0 - Major Radius [m]
  Z0 - Elevation [m]
  ϵ - Inverse Aspect Ratio a/R0 where a = minor radius
  κ - Elongation
  δ - Triangularity
  ζ - Squareness

help?> AsymmetricMillerShape # alias = AMShape

  AsymmetricMillerShape Structure

  Defines the Asymmetric Miller Plasma Shape Parameterization
 
  Fields:
  R0 - Major Radius [m]
  Z0 - Elevation [m]
  ϵ - Inverse Aspect Ratio a/R0 where a = minor radius
  κ - Elongation
  δl - Lower Triangularity
  δu - Upper Triangularity

help?> MillerExtendedHarmonicShape # alias = MXHShape 

  MillerExtendedHarmonicShape Structure

  Defines the Miller Extended Harmonic Plasma Shape Parameterization
  > Arbon, Ryan, Jeff Candy, and Emily A. Belli. "Rapidly-convergent flux-surface shape parameterization." Plasma Physics and Controlled Fusion 63.1 (2020): 012001.

  Fields:
  R0 - Major Radius [m]
  Z0 - Elevation [m]
  ϵ - Inverse Aspect Ratio a/R0 where a = minor radius
  κ - Elongation
  c0 - Tilt
  c - Cosine coefficients i.e. [ovality,...]
  s - Sine coefficients i.e. [asin(triangularity), squareness,...]
```

## Flux Surface Fitting
MXHEquilibrium.jl provides routines for fitting flux surfaces to a PlasmaShape

```julia
julia> g = readg(@__DIR__*"/test/g150219.03200");

julia> M = efit(g,clockwise_phi=false);

julia> bdry = boundary(M);

julia> MXH = fit(bdry,MXHShape(10)) # use 10 fourier components
MillerExtendedHarmonicShape{10, Float64}
  R0 = 1.6667475890195798 [m]
  Z0 = -0.14918543837718545 [m]
  ϵ  = 0.3278211449811159
  κ  = 1.8166930389369045
  δ  = 0.4745349500498265
  ζ  = 0.0640865222160445
  ξ  = -0.12517240241115193
  τ  = -0.0847654065659985

julia> bdry = boundary(MXH; N=100)

julia> fMXH = fit(bdry,MXHShape(10))
MillerExtendedHarmonicShape{10, Float64}
  R0 = 1.6667278385191975 [m]
  Z0 = -0.14918543837718545 [m]
  ϵ  = 0.3277787473904293
  κ  = 1.8167208515952422
  δ  = 0.4745418627448689
  ζ  = 0.06313553088494958
  ξ  = -0.11595260388302828
  τ  = -0.08051086241624195
```

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/MXHEquilibrium.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/MXHEquilibrium.jl/actions/workflows/make_docs.yml/badge.svg)
