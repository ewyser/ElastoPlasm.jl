"""
    get_thermal(mesh, mat, solver::S; ni=2) where {S<:AbstractSolver}

Initialize geometry and material point fields for a thermal test problem.

# Arguments
- `mesh`: Mesh object with geometry and boundary self.
- `mat`: Material parameters (Dict or NamedTuple).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`).
- `ni`: Number of intervals per element (default: 2).

# Returns
- `ni`: Number of intervals per element.
- `nmp`: Number of material points.
- `fields`: NamedTuple with coordinates and material properties.
"""
function get_thermal(mesh,mat,solver; ni = 2)
    props = mesh.prprt
    out = mpts_populate(props,mat,solver; ni=ni)

    xp = copy(out.x)

    nmp  = size(xp,2)

    return (; xp, ni, nmp, material_fields(mat, nmp)...)
end