module ElastoPlasm
# define module location as const
const ROOT = @__DIR__

# include boot file
include(joinpath(ROOT,"boot/boot.jl"))

welcome_log() # TODO(!pending): find a way to force this to be called when `using ElastoPlasm`

end # module elastoPlasm