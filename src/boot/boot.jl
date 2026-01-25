# include dependencies
using Revise,Pkg,Test
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,SparseArrays,Random
using JLD2,HDF5
using KernelAbstractions,Adapt,Base.Threads
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.Kernel as Cairn
import KernelAbstractions.synchronize as sync
import Adapt.adapt as user_adapt
import Adapt.@adapt_structure as @adapt_struct

# include types &
include(joinpath(SRC,"boot/include.jl"))
sucess = superInc(["boot/needs/types"]; root=SRC)

# create primitive structs
info = Self(
    sys = Path(
        root = SRC,
	    out  = joinpath(dirname(SRC),"dump"),
	    test = joinpath(dirname(SRC),"test"),
    ), 
    ui = UI(), 
    bckd = Execution(), 
    mpi = Distributed()
)  

# include
include(joinpath(SRC,"boot/needs/utils.jl"))
include(joinpath(SRC,"boot/needs/backend.jl"))

# flushing
@info join(rootflush(info),"\n")

# find & printout active backend(s)
add_backend!(Val(:x86_64),info)

# include .jl files
lists = ["home/api","home/program","home/script"]
@info join(superInc(lists; root=SRC, lib=info.sys.lib),"\n")