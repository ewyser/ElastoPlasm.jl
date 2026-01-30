module ElastoPlasm
# define module location as const
const SRC = @__DIR__

# include boot file
include(joinpath(SRC,"boot/boot.jl"))

function __init__()
    welcome_log() 
end

end # module ElastoPlasm