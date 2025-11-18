module ElastoPlasm
# define module location as const
const ROOT = @__DIR__

# include boot file
include(joinpath(ROOT,"boot/boot.jl"))

function __init__()
    welcome_log() 
end

end # module ElastoPlasm