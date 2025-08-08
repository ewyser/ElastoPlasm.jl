module ElastoPlasm
# define module location as const
const ROOT = @__DIR__

# include boot file
include(joinpath(ROOT,"boot/boot.jl"))

welcome_log()

end # module elastoPlasm