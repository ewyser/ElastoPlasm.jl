"""
    init_shpfun(dim::Number, instr::NamedTuple; what::String="nothing")

Initialize shape function and topology kernels for the MPM algorithm.

# Arguments
- `dim::Number`: Spatial dimension (1, 2, or 3).
- `instr::NamedTuple`: Instruction/configuration dictionary.
- `what::String`: (Optional) Additional selector for kernel type.

# Returns
- Named tuple of kernel functions for topology, shape function, and delta function.
"""
function init_ignite(instr::NamedTuple; ignite::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}(),)
    backend = exec_backend(instr)   # the solver's device (backend.exec), CPU() for the host
    # topology function
    ignite[:tplgy!]  = p2e2n(backend)
    # cache shape values/gradients for all basis kinds, once per particle per timestep
    ignite[:shpfun!] = shpfun!(backend)
    # nodal solid fraction, for excluding ill-conditioned nodes (`sfrac.jl`)
    ignite[:sfrac!]  = solid_fraction(backend)
    if instr[:basis][:trsfr] == "apic"
        ignite[:Dᵢⱼ!] = Dij_nd(backend) 
    end    
    return (;ignite...)
end
"""
    shpfun(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, instr::NamedTuple) where {T1,T2,D}

Initialize and compute shape functions and topological relations for material points.

# Arguments
- `mpts::Point{T1,T2,D}`: Material point data structure.
- `mesh::Mesh{T1,T2,D}`: Mesh data structure.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function ignite(mpts::Point,mesh::Mesh,basis::Basis,solver::AbstractSolver)
    bk = get_backend(mpts.x)   # the backend the arrays live on
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    solver.cairn.ignite.tplgy!(mpts,mesh,basis; ndrange=(mpts.nmp));sync(bk)
    # cache shape values/gradients for all NN neighbors of every particle
    solver.cairn.ignite.shpfun!(mpts,mesh,basis; ndrange=(mpts.nmp));sync(bk)
    # calculate identity shape functions
    if solver.basis.trsfr == "apic"
        solver.cairn.ignite.Dᵢⱼ!(mpts,mesh,basis; ndrange=(mpts.nmp));sync(bk)
    end
    return nothing
end