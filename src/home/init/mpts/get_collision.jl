"""
    get_collision(mesh::Mesh{T1,T2,D}, mat, solver::S; ni=2, r=0.1, v=20.5) where {D,S<:AbstractSolver}

Initialize geometry and material point fields for a disk collision test problem.

Two circular disks of particles with opposing velocities collide with each other.
This is a classic SPH benchmark problem for testing kernel functions, density
evolution, and momentum conservation.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object with geometry and boundary self.
- `mat`: Material parameters (Dict or NamedTuple).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`).
- `ni`: Number of intervals per element (default: 2).
- `r`: Relative disk radius as fraction of domain height (default: 0.1).
- `v`: Velocity magnitude (default: 20.5).

# Returns
- NamedTuple with:
  - `xp`: Material point coordinates
  - `coh0`: Initial cohesion
  - `cohr`: Residual cohesion
  - `phi`: Friction angle
  - `T`: Temperature
  - `c`: Specific heat capacity
  - `k`: Thermal conductivity
  - `ni`: Number of intervals per element
  - `nmp`: Number of material points
  - `vp`: Velocities (optional for disk collision)

# Example
```julia
geom = get_collision(mesh, mat, solver; ni=3, r=0.1, v=20.5)
mpts = setup_mpts(mesh, solver, mat; geom=geom)
```
"""
function get_collision(mesh::Mesh{T1,T2,D}, mat, solver::S; ni=2, r=0.1, v=10.5) where {T1,T2,D,S<:AbstractSolver}
    props = mesh.prprt
    lz = props.L[end]  # Domain height

    # Populate candidate material point positions
    out = mpts_populate(props, mat, solver; ni=ni)

    # Get candidate positions
    id = findall(x -> x ≤ lz - (0.5*props.h[end]/ni), out.x[end,:])

    if D == 2
        xp, zp = out.x[1, id], out.x[2, id]
        c = out.c0[id]

        # Disk 1 center (left disk)
        xc0 = 0.25 * (props.L[1])
        zc0 = 0.25 * lz

        # Disk 2 center (right disk)
        xc1 = 0.75 * (props.L[1])
        zc1 = 0.75 * lz

        # Velocity direction (from disk 1 to disk 2)
        dx = xc1 - xc0
        dz = zc1 - zc0
        norm = sqrt(dx * dx + dz * dz)
        vx = dx / norm
        vz = dz / norm

        # Disk radius
        radius = r * lz

        # Lists to store material points
        xlt, zlt = Float64[], Float64[]
        vxlt, vzlt = Float64[], Float64[]
        clt = Float64[]

        # Generate particles for disk 1
        for mpts in eachindex(xp)
            Δx = xp[mpts] - xc0
            Δz = zp[mpts] - zc0
            d = sqrt(Δx * Δx + Δz * Δz)

            if d <= radius
                push!(xlt, xp[mpts])
                push!(zlt, zp[mpts])
                push!(vxlt, v * vx)
                push!(vzlt, v * vz)
                push!(clt, c[mpts])
            end
        end

        # Generate particles for disk 2
        for mpts in eachindex(xp)
            Δx = xp[mpts] - xc1
            Δz = zp[mpts] - zc1
            d = sqrt(Δx * Δx + Δz * Δz)

            if d <= radius
                push!(xlt, xp[mpts])
                push!(zlt, zp[mpts])
                push!(vxlt, -v * vx)
                push!(vzlt, -v * vz)
                push!(clt, c[mpts])
            end
        end

        # Combine into position array
        xp = vcat(xlt', zlt')
        vp = vcat(vxlt', vzlt')

    elseif D == 3
        xp, yp, zp = out.x[1, id], out.x[2, id], out.x[3, id]
        c = out.c0[id]

        # Disk 1 center (left disk)
        xc0 = 0.4 * (props.L[1])
        yc0 = 0.5 * (props.L[2])
        zc0 = 0.4 * lz

        # Disk 2 center (right disk)
        xc1 = 0.6 * (props.L[1])
        yc1 = 0.5 * (props.L[2])
        zc1 = 0.6 * lz

        # Velocity direction (from disk 1 to disk 2)
        dx = xc1 - xc0
        dy = yc1 - yc0
        dz = zc1 - zc0
        norm = sqrt(dx * dx + dy * dy + dz * dz)
        vx = dx / norm
        vy = dy / norm
        vz = dz / norm

        # Disk radius
        radius = r * lz

        # Lists to store material points
        xlt, ylt, zlt = Float64[], Float64[], Float64[]
        vxlt, vylt, vzlt = Float64[], Float64[], Float64[]
        clt = Float64[]

        # Generate particles for disk 1 (sphere in 3D)
        for mpts in eachindex(xp)
            Δx = xp[mpts] - xc0
            Δy = yp[mpts] - yc0
            Δz = zp[mpts] - zc0
            d = sqrt(Δx * Δx + Δy * Δy + Δz * Δz)

            if d <= radius
                push!(xlt, xp[mpts])
                push!(ylt, yp[mpts])
                push!(zlt, zp[mpts])
                push!(vxlt, v * vx)
                push!(vylt, v * vy)
                push!(vzlt, v * vz)
                push!(clt, c[mpts])
            end
        end

        # Generate particles for disk 2 (sphere in 3D)
        for mpts in eachindex(xp)
            Δx = xp[mpts] - xc1
            Δy = yp[mpts] - yc1
            Δz = zp[mpts] - zc1
            d = sqrt(Δx * Δx + Δy * Δy + Δz * Δz)

            if d <= radius
                push!(xlt, xp[mpts])
                push!(ylt, yp[mpts])
                push!(zlt, zp[mpts])
                push!(vxlt, -v * vx)
                push!(vylt, -v * vy)
                push!(vzlt, -v * vz)
                push!(clt, c[mpts])
            end
        end

        # Combine into position array
        xp = vcat(xlt', ylt', zlt')
        vp = vcat(vxlt', vylt', vzlt')
    end

    # Number of material points
    nmp = size(xp, 2)

    # Material properties
    coh0 = clt
    cohr = ones(nmp) .* mat[:cr]
    phi = ones(nmp) .* mat[:ϕ0]

    # Thermal properties
    c = ones(nmp) .* mat[:specific_heat_capacity]
    k = ones(nmp) .* mat[:thermal_conductivity]
    T = ones(nmp) .* mat[:initial_temperature]

    return (;xp=xp, vp=vp, coh0=coh0, cohr=cohr, phi=phi, T=T, c=c, k=k, ni=ni, nmp=nmp)
end
