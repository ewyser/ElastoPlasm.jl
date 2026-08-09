export ExplicitSolver,ImplicitSolver,DynamicRelaxationSolver,Instruction, Time

struct ExplicitSolver{T1<:Integer,T2<:Real,D<:AbstractDimension} <: AbstractSolver{T1,T2,D}
    dtype  ::NamedTuple
    basis  ::NamedTuple
    fwrk   ::NamedTuple
    bcs    ::NamedTuple
    grf    ::NamedTuple
    plast  ::NamedTuple
    nonloc ::NamedTuple
    plot   ::NamedTuple
    perf   ::NamedTuple
    backend::NamedTuple
    cairn  ::NamedTuple
end
@adapt_struct ExplicitSolver

struct DynamicRelaxationSolver{T1<:Integer,T2<:Real,D<:AbstractDimension} <: AbstractSolver{T1,T2,D} 
    dtype  ::NamedTuple
    basis  ::NamedTuple
    fwrk   ::NamedTuple
    bcs    ::NamedTuple
    grf    ::NamedTuple
    plast  ::NamedTuple
    nonloc ::NamedTuple
    plot   ::NamedTuple
    perf   ::NamedTuple
    backend::NamedTuple
    cairn  ::NamedTuple
end
@adapt_struct DynamicRelaxationSolver

struct ImplicitSolver{T1<:Integer,T2<:Real,D<:AbstractDimension} <: AbstractSolver{T1,T2,D}
    dtype  ::NamedTuple
    basis  ::NamedTuple
    fwrk   ::NamedTuple
    bcs    ::NamedTuple
    grf    ::NamedTuple
    plast  ::NamedTuple
    nonloc ::NamedTuple
    plot   ::NamedTuple
    perf   ::NamedTuple
    backend::NamedTuple
    cairn  ::NamedTuple
end
@adapt_struct ImplicitSolver






export Instruction, Time

struct Instruction{T1<:Integer,T2<:Real,D<:AbstractDimension}
    dtype  ::NamedTuple
    basis  ::NamedTuple
    fwrk   ::NamedTuple
    bcs    ::NamedTuple
    grf    ::NamedTuple
    plast  ::NamedTuple
    nonloc ::NamedTuple
    plot   ::NamedTuple
    perf   ::NamedTuple
    backend::NamedTuple
    cairn  ::NamedTuple
end
@adapt_struct Instruction



struct Time{T1<:Integer,T2<:Real}
    t  ::Vector{T2}
    te ::T2
    tg ::T2
    tep::T2
end
@adapt_struct Time
