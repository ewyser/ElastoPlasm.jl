abstract type AbstractSolver{T1<:Integer,T2<:Real,D} end

export ExplicitSolver,ImplicitSolver,DynamicRelaxationSolver

struct ExplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    material ::MT
    nonloc   ::NL
    plot     ::PT
    perf     ::PF
    backend  ::BK
    cairn    ::CN
end
function ExplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN}
    return ExplicitSolver{T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN}(solution,dtype,basis,stab,bcs,grf,material,nonloc,plot,perf,backend,cairn)
end
@adapt_struct ExplicitSolver

struct DynamicRelaxationSolver{T1<:Integer,T2<:Real,D,DT,BA,FW,BC,GR,MT,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    dtype  ::DT
    basis  ::BA
    fwrk   ::FW
    bcs    ::BC
    grf    ::GR
    material::MT
    nonloc ::NL
    plot   ::PT
    perf   ::PF
    backend::BK
    cairn  ::CN
end
function DynamicRelaxationSolver{T1,T2,D}(dtype::DT,basis::BA,fwrk::FW,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,FW,BC,GR,MT,NL,PT,PF,BK,CN}
    return DynamicRelaxationSolver{T1,T2,D,DT,BA,FW,BC,GR,MT,NL,PT,PF,BK,CN}(dtype,basis,fwrk,bcs,grf,material,nonloc,plot,perf,backend,cairn)
end
@adapt_struct DynamicRelaxationSolver

struct ImplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    material ::MT
    nonloc   ::NL
    plot     ::PT
    perf     ::PF
    backend  ::BK
    cairn    ::CN
end
function ImplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,stab::SB,bcs::BC,grf::GR,material::MT,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN}
    return ImplicitSolver{T1,T2,D,DT,BA,SB,BC,GR,MT,NL,PT,PF,BK,CN}(solution,dtype,basis,stab,bcs,grf,material,nonloc,plot,perf,backend,cairn)
end
@adapt_struct ImplicitSolver
