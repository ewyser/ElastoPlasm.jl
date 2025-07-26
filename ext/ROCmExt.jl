__precompile__(false)

module ROCmExt

@info "📦 Including ROCmExt.jl extension module"

using ElastoPlasm

try
    @info "🔧 Using ROCm backend"
    using AMDGPU
    @info "🧠 ROCm 🔁 overloading stub functions..."
    #include(joinpath(@__DIR__, "ROCmExt", "ROCm_backend.jl"))
    #add_backend!(Val(:AMDGPU), info)
catch
    @info "🧠 ROCm loaded, but no ROCm backend found..."
end

end
