# include dependencies
using Revise,Pkg,Test,Suppressor
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,SparseArrays,Random
using KernelAbstractions,Adapt,Base.Threads
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
import Adapt.adapt as user_adapt
import Adapt.@adapt_structure as @user_struct

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
@info join(superList(lists),"\n")


#=
sys = moduleCore()
out = ["ElastoPlasm.jl location:\n\t- "*sys.root,]
if !isdir(sys.out)
	mkpath(sys.out) 
	push!(out,"\n\tcreating directory at:\n- "*sys.out)
else
	push!(out,"\n\talready existing directory at:\n- "*sys.out)
end
=#