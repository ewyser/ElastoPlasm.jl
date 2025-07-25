module ElastoPlasm
# define module location as const
const ROOT = @__DIR__

# include boot file
include(joinpath(ROOT,"boot/boot.jl"))

@info """New comer ?
- copy-paste the following:
  instr = slump([64.1584,64.1584],40);
- wait for the simulation to end
"""

end # module elastoPlasm