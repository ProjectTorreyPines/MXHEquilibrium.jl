using Test

using ForwardDiff
using MXHEquilibrium
using CoordinateConventions
using EFIT


mu0 = 4pi*10^-7

@testset "Equilibrium" begin

include("shape.jl")

include("solovev.jl")

include("efit.jl")

end
