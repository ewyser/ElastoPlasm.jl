I will explain in greater details in this file what the agent should do.

I designed the solver and implemented the following logic for memory layout

in .\ElastoPlasm.jl\src\boot\needs\types\concrete\lagrangian.jl, the agent will find all the struct related to mps whereas in .\ElastoPlasm.jl\src\boot\needs\types\concrete\eulerian.jl, the agent will find all the sruct related to mesh

My problem is the following: in the function ElastoPlasm.jl\src\home\core\workflow\implicit\fint.jl, there is a lot of small-size matrices. The current memory layout seems not optimal to me, for instance:

    Fn = SMatrix{2,2,T2}(mpts.s.Fn[:, :, mp])
    ϵn = SMatrix{2,2,T2}(mpts.s.ϵn[:, :, mp])

would be better as

    Fn = SMatrix{2,2,T2}(mpts.s.Fn[mp])
    ϵn = SMatrix{2,2,T2}(mpts.s.ϵn[mp])

where mpts.s.ϵn[mp] = SMatrix

In general, i want to use a memory layout as such: 

mpts.s.Fn = SVector{nmp,SMatrix} with SMatrix{2,2,Float32} if TwoDimension

I want to migrate toward static vectors as containers 

The agent should modify accordingly all the necessary file in /src following this new logic

The same is true for mesh-realted quantities. Take mesh.x for example. mesh.x = SVector{nmp, SVector}

As of inspiration, FEMTools package .\FEMTools.jl should be carefully used and mimic for that specific concern of memory layout.