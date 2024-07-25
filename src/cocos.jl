# ----- GEQDSK COCOS Interface -----

"""
    check_cocos(g::GEQDSKFIle, cc::Union{Int,COCOS}) -> Bool
Checks if `GEQDSKFile` is consistant with the given cocos, `cc`
"""
function check_cocos(g::GEQDSKFile, cc::Union{Int,COCOS}; kwargs...)
    return check_cocos(g.bcentr, g.current, g.fpol, g.pprime, g.qpsi, g.psi, cocos(cc); kwargs...)
end

"""
    identify_cocos(g::GEQDSKFIle; clockwise_phi=nothing) -> List of possible COCOS IDs
Identifies possible `GEQDSKFile` COCOS. A unique identification requires setting the `clockwise_phi` keyword.
"""
function identify_cocos(g::GEQDSKFile; clockwise_phi = nothing)
    return filter(x -> x < 10, identify_cocos(g.bcentr, g.current, g.qpsi, g.psi, clockwise_phi, nothing))
end

"""
    cocos(g::GEQDSKFile; kwargs...) -> COCOS
Identifies and returns `GEQDSKFile` `COCOS`. A unique identification requires setting the `clockwise_phi` keyword.
"""
function cocos(g::GEQDSKFile; kwargs...)
    cc = identify_cocos(g; kwargs...)
    if length(cc) > 1
        @error "Unable to determine unique COCOS. Try providing clockwise_phi::Bool keyword. Possibilities are $cc"
    end
    return cocos(cc[1])
end

"""
    transform_cocos(g::GEQDSKFile, cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS}; kwargs...) -> GEQDSKFile
Transforms the given `GEQDSKFile` with `COCOS=cc_in`, and returns a `GEQDSKFile` with `COCOS=cc_out`.
"""
function transform_cocos(g::GEQDSKFile, cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS}; kwargs...)
    T = transform_cocos(cc_in, cc_out; kwargs...)

    g_new = GEQDSKFile(g.file*" w/ cocos = $(cocos(cc_out).cocos)", g.time, g.nw, g.nh,
                       g.r*T["R"], g.z*T["Z"], g.rdim*T["R"], g.zdim*T["Z"],
                       g.rleft*T["R"], g.zmid*T["Z"], g.nbbbs, g.rbbbs*T["R"], g.zbbbs*T["Z"],
                       g.limitr, g.rlim*T["R"], g.zlim*T["Z"], g.rcentr*T["R"], g.bcentr*T["B"],
                       g.rmaxis*T["R"], g.zmaxis*T["Z"], g.simag*T["PSI"], g.sibry*T["PSI"],
                       g.psi*T["PSI"], g.current*T["I"], g.fpol*T["F"],
                       g.pres*T["P"], g.ffprim*T["F_FPRIME"], g.pprime*T["PPRIME"],
                       g.qpsi*T["Q"], g.psirz*T["PSI"], g.rhovn)

    return g_new
end
