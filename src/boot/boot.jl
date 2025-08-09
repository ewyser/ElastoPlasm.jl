# include dependencies
using Revise,Pkg,Test
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,SparseArrays,Random
using KernelAbstractions,Adapt,Base.Threads
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
import Adapt.adapt as user_adapt
import Adapt.@adapt_structure as @adapt_struct

# arithmetic precision & relative path for figs & data
const typeD = Float64  

# create primitive structs
include(joinpath(ROOT,"boot/needs/types.jl"))
info = Self(sys = Path(), ui = UI(), bckd = Execution(), mpi = Distributed())  

# include
include(joinpath(ROOT,"boot/needs/utils.jl"))
include(joinpath(ROOT,"boot/needs/backend.jl"))

lists = ["home/api","home/program","home/script"]

# flushing
@info join(rootflush(info),"\n")

# find & printout active backend(s)
add_backend!(Val(:x86_64),info)

# include .jl files
@info join(superInc(lists),"\n")