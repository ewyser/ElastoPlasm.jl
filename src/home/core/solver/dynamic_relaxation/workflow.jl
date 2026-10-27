export elastoquasistatic!

"""
    elastoquasistatic!(mpts, mesh, basis, time, solver; lstps=10, tol=1e-9, plastic=false, G=9.81)

Quasi-static self-weight loading by implicit dynamic relaxation (Sutcliffe et al. 2026,
Algorithm 1): gravity is applied in `lstps` uniform load steps, `g = (k/lstps)·G` at step
`k`, each solved to `‖R‖/‖fᵉˣᵗ‖ < tol` by `relax!`. Returns the DR iteration count of each
load step. `plastic=true` adds the `material.plastic` return mapping inside every DR
iteration (non-local regularization is not applied). Finite strain only. To change `lstps` or
`tol` from `execute!`, pass a closure, e.g.
`workflows = [(a...) -> elastoquasistatic!(a...; lstps=40)]`.
"""
function elastoquasistatic!(mpts::Point{T1,T2,D,CM,ST},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},time::Time{T1,T2},instr::S; lstps::Int=10, tol::Real=1e-9, plastic::Bool=false, G::Real=9.81) where {T1,T2,D,CM,ST,S<:AbstractSolver{T1,T2,D}}
    ST <: LogarithmicStrain || throw(ArgumentError("elastoquasistatic! is finite-strain only: material.elastic=\"hypoelastic\" is not supported"))
    its  = zeros(Int, lstps)
    sf   = node_exclusion(mpts,mesh,instr)
    prog = Progress(lstps;dt=0.5,desc="Solving elastoquasistatic...",barlen=10)
    for k ∈ 1:lstps
        g  = get_g(mesh.prprt; G=T2(G*k/lstps))
        ignite(mpts,mesh,basis,instr)
        its[k] = it = relax!(mpts,mesh,basis,g,instr; tol=T2(tol), plastic=plastic, sf=sf)
        instr.cairn.implicit.update!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        time.t[1] += T2(1.0)
        bake(mpts,mesh,time.t[1],instr; title="load step $k/$lstps")
        record!(mpts,mesh,time,instr; step=k, loadstep=true)
        next!(prog;showvalues = [("load step",k),("DR iterations",it)])
    end
    finish!(prog)
    return its
end

# the workflow on a whole problem (`execute!(sim)`)
elastoquasistatic!(p::Problem; kwargs...) = elastoquasistatic!(p.mpts, p.mesh, p.basis, p.time, p.solver; kwargs...)
