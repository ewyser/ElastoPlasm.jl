export elastoquasistatic!  

function elastoquasistatic!(mpts::Point{T1,T2,D,NN,<:AbstractBasis,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::S) where {T1,T2,D,NN,E,R,S<:AbstractSolver{T1,T2,D}}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoquasistatic...",barlen=10)
    lstps = collect(range(0, 9.8; length=10))
    for (it,lstp) ∈ enumerate(lstps)

        tic      = time_ns()


        dt = T2(1.0)
        g  = T2[T2(0), -lstp]

        ignite(mpts,mesh,instr)
        relax(mpts,mesh,cmpr,g,dt,instr)
        instr.cairn.implicit.update!(mpts,mesh,dt; ndrange=mpts.nmp);sync(CPU())

        time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)

        savlot(mpts,mesh,time.t[1],instr)
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end