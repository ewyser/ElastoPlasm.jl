
"""
    build_solid_mesh_phase(T1,T2,D,prop,bcs,nn_total,SVD) -> MeshSolidPhase

Build the nodal solid-phase state (`MeshSolidPhase`) — mass, out-of-balance force, momentum,
velocity, acceleration — all zero-initialized.
"""
function build_solid_mesh_phase(T1,T2,D,prop,bcs,nn_total,SVD)
    return MeshSolidPhase{T1,T2,D}(
        prop,
        bcs,
        T2.(zeros(nn_total             )), # m
        T2.(zeros(D,nn_total          )), # oobf  (Matrix: atomic writes)
        [zero(SVD) for _ in 1:nn_total]  , # a
        T2.(zeros(D,nn_total          )), # mv    (Matrix: atomic writes)
        [zero(SVD) for _ in 1:nn_total]  , # v
    )
end

"""
    build_thermal_mesh_phase(T1,T2,D,prop,bcs,nn_total; thermal::Bool=false) -> Union{Nothing,MeshThermalPhase}

Build the nodal thermal-phase state (`MeshThermalPhase`), zero-initialized and reusing
the same Dirichlet boundary layout as the solid phase, when `thermal=true` — otherwise
`nothing`.
"""
function build_thermal_mesh_phase(T1,T2,D,prop,bcs,nn_total; thermal::Bool=false)
    return thermal ? MeshThermalPhase{T1,T2,D}(
        prop,
        bcs,
        T2.(zeros(nn_total)), # cᵢ
        T2.(zeros(nn_total)), # oobq
        T2.(zeros(nn_total)), # dT
        T2.(zeros(nn_total)), # mcT
        T2.(zeros(nn_total)), # T
    ) : nothing
end

"""
    mesh_coords(L, h, nel) -> (x, nel, nno)

The nodes of the grid `[0, L]` with spacing `h`: coordinates (D×nno, x fastest, then y, then z) and
the element and node counts per axis followed by their product. A floating-point range that misses a
node (count ≠ `nel + 1`) is an error.
"""
function mesh_coords(L::AbstractVector{T2}, h::AbstractVector{T2}, nel) where {T2}
    axs = [collect(zero(T2):h[d]:L[d]) for d ∈ eachindex(L)]
    n   = length.(axs)
    n == nel .+ 1 || error("the node ranges have $(n) nodes for $(nel) elements (floating-point range of L/nel)")
    x   = reshape(stack(Iterators.product(axs...)), length(L), :)
    return x, [n .- 1; prod(n .- 1)], [n; prod(n)]
end

"""
    setup_mesh(L, nel, solver; thermal=false) -> Mesh

The background grid `[0, L]` with `nel` elements per axis: node coordinates, element size,
boundary conditions (`solver.bcs.dirichlet`, `get_bc`) and
zero nodal fields of the solid phase and, with `thermal = true`, the thermal phase. Connectivity is
built by `setup_basis`, once the material points exist.
"""
function setup_mesh(L::AbstractVector, nel::AbstractVector, solver::S; thermal::Bool=false) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    h          = T2.(L ./ nel)
    xn,nel,nno = mesh_coords(T2.(L), h, nel)
    status,xB  = get_bc(xn,solver)
    SVD = SVector{D,T2}                         # node coordinates, velocity, acceleration
    nn_total = nno[end]

    prop = MeshProperties{T1,T2,D}(
        T1.(nel                         ), # nel
        T1.(nno                         ), # nno
        T2.(L                           ), # L
        SVD(h                           ), # h
        T2.(xB                          ), # xB
    )
    bcs = MeshBoundary(
        status
    )
    s = build_solid_mesh_phase(T1,T2,D,prop,bcs,nn_total,SVD)
    t = build_thermal_mesh_phase(T1,T2,D,prop,bcs,nn_total; thermal=thermal)
    mesh = Mesh{T1,T2,D}(
        prop,
        # nodal quantities
        SVD(T2.(vec(minimum(xn,dims=2))) )     , # x₀
        [SVD(T2.(xn[:,i])) for i in 1:nn_total], # x
        T2.(zeros(nn_total                  )), # ΔJ
        # solid phase
        s                                , # solid phase
        # thermal phase
        t                                , # thermal phase
    )
    return mesh
end
