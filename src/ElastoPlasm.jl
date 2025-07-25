module ElastoPlasm
# define module location as const
const ROOT = @__DIR__

# include boot file
include(joinpath(ROOT,"boot/boot.jl"))

@info """New comer ?
- copy-paste the following:
  L,nel = [64.1584,64.1584],[40,40];
  instr = slump(L,nel);
- wait for the simulation to end
"""

end # module elastoPlasm