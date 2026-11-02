abstract type AbstractSolver{T1<:Integer,T2<:Real,D} end
# how a case is solved: `Solver(case)` builds one (an `ExplicitSolver` or an `ImplicitSolver`), see its docstring
const Solver = AbstractSolver

export Solver,ExplicitSolver,ImplicitSolver,DynamicRelaxationSolver

struct ExplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    material ::MT
    nonloc   ::NL
    plot     ::PT
    output   ::OT
    backend  ::BK
    run      ::RN       # the run settings (`case_run`): te, tg, tep, workflow, lstps, folder
    case     ::CS       # the case solved (`Solver(case)`); `nothing` for `get_solver`'s and a stored solver
    cairn    ::CN
end
function ExplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,case::CS,cairn::CN) where {T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN}
    return ExplicitSolver{T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN}(solution,dtype,basis,stab,bcs,grf,material,nonloc,plot,output,backend,run,case,cairn)
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

struct ImplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    material ::MT
    nonloc   ::NL
    plot     ::PT
    output   ::OT
    backend  ::BK
    run      ::RN       # the run settings (`case_run`): te, tg, tep, workflow, lstps, folder
    case     ::CS       # the case solved (`Solver(case)`); `nothing` for `get_solver`'s and a stored solver
    cairn    ::CN
end
function ImplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,output::OT,backend::BK,run::RN,case::CS,cairn::CN) where {T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN}
    return ImplicitSolver{T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,OT,BK,RN,CS,CN}(solution,dtype,basis,stab,bcs,grf,material,nonloc,plot,output,backend,run,case,cairn)
end
