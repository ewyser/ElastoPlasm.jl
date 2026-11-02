abstract type AbstractSolver{T1<:Integer,T2<:Real,D} end
# how a case is solved: `Solver(case)` builds one (an `ExplicitSolver` or an `ImplicitSolver`), see its docstring
const Solver = AbstractSolver

export Solver,ExplicitSolver,ImplicitSolver,DynamicRelaxationSolver

struct ExplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,NL,PT,OT,BK,RN,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    nonloc   ::NL
    plot     ::PT
    output   ::OT
    backend  ::BK
    run      ::RN       # the run settings (`solver_run`): te, tg, tep, workflow, lstps
    cairn    ::CN
end
function ExplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,cairn::CN) where {T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    return ExplicitSolver{T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}(solution,dtype,basis,stab,nonloc,plot,output,backend,run,cairn)
end

struct DynamicRelaxationSolver{T1<:Integer,T2<:Real,D,DT,BA,FW,BC,GR,MT,NL,PT,BK,CN} <: AbstractSolver{T1,T2,D}
    dtype  ::DT
    basis  ::BA
    fwrk   ::FW
    bcs    ::BC
    grf    ::GR
    material::MT
    nonloc ::NL
    plot   ::PT
    backend::BK
    cairn  ::CN
end
function DynamicRelaxationSolver{T1,T2,D}(dtype::DT,basis::BA,fwrk::FW,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,FW,BC,GR,MT,NL,PT,BK,CN}
    return DynamicRelaxationSolver{T1,T2,D,DT,BA,FW,BC,GR,MT,NL,PT,BK,CN}(dtype,basis,fwrk,bcs,grf,material,nonloc,plot,backend,cairn)
end

struct ImplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,NL,PT,OT,BK,RN,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    nonloc   ::NL
    plot     ::PT
    output   ::OT
    backend  ::BK
    run      ::RN       # the run settings (`solver_run`): te, tg, tep, workflow, lstps
    cairn    ::CN
end
function ImplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,cairn::CN) where {T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}
    return ImplicitSolver{T1,T2,D,DT,BA,SB,NL,PT,OT,BK,RN,CN}(solution,dtype,basis,stab,nonloc,plot,output,backend,run,cairn)
end
