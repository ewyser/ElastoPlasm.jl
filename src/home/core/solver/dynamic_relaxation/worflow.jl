export elastoquasistatic!  

function elastoquasistatic!(mpts::Point{T1,T2,D,CM},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},time::Time{T1,T2},instr::S) where {T1,T2,D,CM,S<:AbstractSolver{T1,T2,D}}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoquasistatic...",barlen=10)
    lstps = collect(range(0, 9.8; length=10))
    for (it,lstp) ∈ enumerate(lstps)

        tic      = time_ns()


        dt = T2(1.0)
        g  = T2[T2(0), -lstp]

        ignite(mpts,mesh,basis,instr)
        relax(mpts,mesh,basis,g,dt,instr)
        instr.cairn.implicit.update!(mpts,mesh,basis,dt; ndrange=mpts.nmp);sync(CPU())

        time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)

        bake(mpts,mesh,time.t[1],instr)
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end

# how to operate end-to-end simulation with elastoquasistatic! workflow
# L,nel  = [64.1584,64.1584/4.0],[40,10];
# jld2   = slump_problem(L,nel;cli()...);
# out    = elastoplasm(jld2; workflows = [elastoquasistatic!]);