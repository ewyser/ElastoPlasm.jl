# Include dependencies
using Revise,Pkg,Test
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,StaticArrays,SparseArrays,Random
using JLD2,HDF5
using KernelAbstractions,Adapt,Base.Threads
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.Kernel as Cairn
import KernelAbstractions.synchronize as sync
import Adapt.adapt as user_adapt
import Adapt.@adapt_structure as @adapt_struct

# Include types 
include(joinpath(SRC,"boot/include.jl"))
for f ∈ [
    "self.jl",
    "solver.jl",
    "problem/constitutive.jl",
    "problem/geometry.jl",
    "problem/eulerian.jl",
    "problem/tensor.jl",
    "problem/lagrangian.jl",
    "problem/time.jl",
    "problem/problem.jl",
    "basis/basis.jl",
    "basis/bsmpm.jl",
    "basis/gimpm.jl",
    "basis/mlsmpm.jl",
    "basis/smpm.jl",
    "basis/transfer.jl",
]
    include(joinpath(SRC,"boot/needs/types",f))
end

# Include utility and backend files
include(joinpath(SRC,"boot/needs/utils.jl"))
include(joinpath(SRC,"boot/needs/backend.jl"))
include(joinpath(SRC,"boot/needs/distributed.jl"))

# Create self object
self = Self(
    sys = Path(
        root = SRC,
	    dump  = joinpath(dirname(SRC),"dump"),
	    test = joinpath(dirname(SRC),"test"),
    ), 
    ui = UI(), 
    bckd = Execution(), 
    mpi = Distributed()
)  

# Flushing
rootflush(self.sys.dump)

# Find & printout active backend(s)
add_backend!(self.bckd, Val(:x86_64))

# Automatically include .jl files
lists = ["home/api","home/init","home/plot","home/utils","home/script","home/core"]
@info join(superInc(lists; root=SRC, lib=self.sys.lib),"\n")