__precompile__(false)

module CUDAExt

@info "📦 Including CUDAExt.jl extension module"

using ElastoPlasm

try
    @info "🔧 Using CUDA backend"
    using CUDA
    @info "🧠 CUDA 🔁 overloading stub functions..."
    #include(joinpath(@__DIR__, "CUDAExt", "CUDA_backend.jl"))
    #add_backend!(Val(:CUDA), info)
catch
    @info "🧠 CUDA loaded, but no CUDA backend found..."
end

end