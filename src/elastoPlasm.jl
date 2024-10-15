module elastoPlasm
# define module location as const
const ROOT = dirname(@__FILE__)
# include startup file
include(joinpath(ROOT,"init/startup.jl"))
# include .jl files
sucess = ["welcome to elastoPlasm:\nsucessful superInclude()"]
for (k,child) ∈ enumerate(sys.lib)
	list = superInc(joinpath(sys.init,child))
	if isempty(list)
		push!(sucess,"\n✗ "*child)
	else	
		push!(sys.method,("$(child)"=>list))
		push!(sucess    ,"\n✓ "*child      )
	end
	if haskey(ENV,"TREE") && ENV["TREE"]=="true"
		push!(sucess,join(tree(list)))
	end
end
@info join(sucess)
@info """new comer ?
- copy-paste the followings:
  julia> L,nel = [64.1584,12.80],40
  julia> instr = slump(L,nel);
- wait for the simulation to end
"""
end # module elastoPlasm

