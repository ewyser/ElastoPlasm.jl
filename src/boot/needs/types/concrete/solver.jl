export ExplicitSolver,ImplicitSolver,DynamicRelaxationSolver, Time

struct ExplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    strain   ::ST
    transfer ::TR
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    plast    ::PL
    nonloc   ::NL
    plot     ::PT
    perf     ::PF
    backend  ::BK
    cairn    ::CN
end
function ExplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,strain::ST,transfer::TR,stab::SB,bcs::BC,grf::GR,plast::PL,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN}
    return ExplicitSolver{T1,T2,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN}(solution,dtype,basis,strain,transfer,stab,bcs,grf,plast,nonloc,plot,perf,backend,cairn)
end
@adapt_struct ExplicitSolver

struct DynamicRelaxationSolver{T1<:Integer,T2<:Real,D,DT,BA,FW,BC,GR,PL,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    dtype  ::DT
    basis  ::BA
    fwrk   ::FW
    bcs    ::BC
    grf    ::GR
    plast  ::PL
    nonloc ::NL
    plot   ::PT
    perf   ::PF
    backend::BK
    cairn  ::CN
end
function DynamicRelaxationSolver{T1,T2,D}(dtype::DT,basis::BA,fwrk::FW,bcs::BC,grf::GR,plast::PL,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,FW,BC,GR,PL,NL,PT,PF,BK,CN}
    return DynamicRelaxationSolver{T1,T2,D,DT,BA,FW,BC,GR,PL,NL,PT,PF,BK,CN}(dtype,basis,fwrk,bcs,grf,plast,nonloc,plot,perf,backend,cairn)
end
@adapt_struct DynamicRelaxationSolver

struct ImplicitSolver{T1<:Integer,T2<:Real,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN} <: AbstractSolver{T1,T2,D}
    solution ::String
    dtype    ::DT
    basis    ::BA
    strain   ::ST
    transfer ::TR
    stab     ::SB
    bcs      ::BC
    grf      ::GR
    plast    ::PL
    nonloc   ::NL
    plot     ::PT
    perf     ::PF
    backend  ::BK
    cairn    ::CN
end
function ImplicitSolver{T1,T2,D}(solution::String,dtype::DT,basis::BA,strain::ST,transfer::TR,stab::SB,bcs::BC,grf::GR,plast::PL,nonloc::NL,plot::PT,perf::PF,backend::BK,cairn::CN) where {T1,T2,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN}
    return ImplicitSolver{T1,T2,D,DT,BA,ST,TR,SB,BC,GR,PL,NL,PT,PF,BK,CN}(solution,dtype,basis,strain,transfer,stab,bcs,grf,plast,nonloc,plot,perf,backend,cairn)
end
@adapt_struct ImplicitSolver



struct Time{T1<:Integer,T2<:Real}
    t  ::Vector{T2}
    te ::T2
    tg ::T2
    tep::T2
end
@adapt_struct Time
