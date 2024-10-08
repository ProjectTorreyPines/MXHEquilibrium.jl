using NetCDF

function transp_potential!(M, transp_file, time)
    tt = Float64.(ncread(transp_file, "TIME"))
    ind = argmin(abs.(tt .- 1e-3 * time))
    psi = Float64.(ncread(transp_file, "PLFLX")[:, ind])
    epot = Float64.(ncread(transp_file, "EPOTNC")[:, ind])
    phi_itp = linear_interpolation(psi, epot; extrapolation_bc=Flat())
    phi = phi_itp.(M.psi)
    phi_itp = cubic_spline_interpolation(M.psi, phi; extrapolation_bc=Flat())
    M.phi = phi_itp
    return nothing
end
