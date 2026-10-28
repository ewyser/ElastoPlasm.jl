# Include dependencies
using Revise,Pkg,Test
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,StaticArrays,SparseArrays,Random
using JLD2,HDF5
using TOML
using EasyDraw: EasyDraw, AbstractShape, Polygon, Block, Ball, Region, Zone, inside, inpolygon, polygon_rings, shape_dict, parse_shape, outside,
                zone_at, zone_index, property, points, partition, write_geo, save_layout, write_catalog
using KernelAbstractions,Adapt,Base.Threads
import UnifiedBackend    # the execution backends; qualified: KernelAbstractions also exports Backend, get_backend
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
    "problem/eulerian.jl",
    "problem/tensor.jl",
    "problem/elasticity.jl",
    "problem/lagrangian.jl",
    "problem/time.jl",
    "problem/case.jl",
    "basis/basis.jl",
    "basis/bsmpm.jl",
    "basis/gimpm.jl",
    "basis/mlsmpm.jl",
    "basis/smpm.jl",
    "basis/transfer.jl",
    "problem/problem.jl",          # needs Point, Mesh, Basis, Time, Case
    "problem/simulation.jl",       # needs Problem
]
    include(joinpath(SRC,"boot/needs/types",f))
end

# Include utility files
include(joinpath(SRC,"boot/needs/utils.jl"))
include(joinpath(SRC,"boot/needs/distributed.jl"))

# Create self object
self = Self(
    sys = Path(
        root = SRC,
	    dump  = joinpath(dirname(SRC),"dump"),
	    test = joinpath(dirname(SRC),"test"),
    ), 
    ui = UI(), 
    mpi = Distributed()
)  

# the default output folder (not tracked by git), emptied only on request: clean_dump()
mkpath(self.sys.dump)

# Automatically include .jl files
lists = ["home/api","home/init","home/plot","home/utils","home/script","home/core"]
@info join(superInc(lists; root=SRC, lib=self.sys.lib),"\n")