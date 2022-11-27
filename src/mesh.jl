import Gmsh
using FileIO, MeshIO
import Meshes

"""
    create_mesh(ext::Boundary; interior=Boundary[], side_length=circumference(ext)/50, interior_side_length=fill(side_length,length(interior)))

Create a Triangular mesh within `ext` using Gmsh. Interior boundaries define holes. Returns Meshes.SimpleMesh
"""
function create_mesh(ext::Boundary; interior=Boundary[], side_length=circumference(ext)/50, interior_side_length=fill(side_length,length(interior)))
    Gmsh.initialize(["-v","0"])

    name = tempname()
    Gmsh.gmsh.model.add(name)

    factory = Gmsh.gmsh.model.geo

    p_t = 0  # point tag
    l_t = 0  # line/spline/curve tag
    c_t = 0  # curve loop tag
    s_t = 0  # surface tag

    # Add interior curves
    Ntot = 0
    for (j,h) in enumerate(interior)
        hN = length(h)
        for i=1:(hN-1)
            p = h.points[i]
            p_t += 1
            factory.addPoint(p[1], p[2], 0.0, interior_side_length[j], p_t)
        end
        l_t += 1
        factory.addSpline(append!(collect((Ntot+1):p_t),(Ntot+1)),l_t)
        c_t += 1
        factory.addCurveLoop([l_t],c_t)
        Ntot += hN-1
    end

    # Add Exterior Surface
    N = length(ext.points)
    for i = 1:(N-1)
        p = ext.points[i]
        p_t += 1
        factory.addPoint(p[1], p[2], 0.0, side_length, p_t)
    end
    l_t += 1
    factory.addSpline(append!(collect((Ntot+1):p_t),(Ntot+1)),l_t)
    c_t += 1
    factory.addCurveLoop([l_t],c_t)
    s_t += 1 
    factory.addPlaneSurface(collect(c_t:-1:1),s_t)

    factory.synchronize()

    Gmsh.gmsh.model.mesh.generate(2)
    
    filename = name*".msh"
    Gmsh.gmsh.write(filename)
    Gmsh.finalize()
    
    m = load(filename)

    # Convert to Meshes.jl format
    points = [Tuple([p[1],p[2]]) for p in Set(m.position)]
    indices = Dict(p=>i for (i,p) in enumerate(points))
    connect = map(m) do el
        Meshes.connect(Tuple(indices[Tuple([p[1],p[2]])] for p in el))
    end
    return Meshes.SimpleMesh(points,connect)
end

"""
    create_mesh(ext::PlasmaShape; kwargs...)

Create a mesh from a plasma shape.
"""
function create_mesh(ext::PlasmaShape; kwargs...)
    b = plasma_boundary(ext; N=50)
    return create_mesh(b; kwargs...)
end
