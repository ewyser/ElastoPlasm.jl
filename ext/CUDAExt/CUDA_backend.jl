"""
    add_backend!(::Val{:CUDA},info::Self)

Description:
---
Return Dicts of gpu backend based on hard-coded supported backends. 
"""
function add_backend!(::Val{:CUDA},info::Self)
    availables::Dict{Symbol, Dict{Symbol, Any}} = Dict( 
        :CUDA => Dict(
            :dev       => "gpu",
            :Backend   => CUDABackend(),
            :brand     => "NVIDIA",
            :wrapper   => CuArray,
            :devices   => CUDA.devices,
            :name      => CUDA.name,
            :handle    => CuDevice,
            :functional => CUDA.functional(),
        )
    )
    dev_id = length(keys(info.bckd.gpu))
    for (platform, backend) ∈ availables
        if backend[:functional]
            for (k, device) ∈ enumerate(backend[:devices]())
                dev_id += 1
                info.bckd.gpu[Symbol("dev$(dev_id)")] = (; 
                    dev      = backend[:dev],
                    platform = platform,
                    brand    = backend[:brand],
                    name     = backend[:name](device),
                    Backend  = backend[:Backend],
                    wrapper  = backend[:wrapper],
                    handle   = device,
                )
            end
            push!(info.bckd.functional, "✓ $(backend[:brand]) $(platform)")
        else
            push!(info.bckd.functional, "✗ $(backend[:brand]) $(platform)")
        end
    end
    @info join(info.bckd.functional,"\n")
    return nothing
end

"""
    device_wakeup!(handle::CuDevice)

Description:
---
Return Dicts of effective cpu and gpu backend based on hard-coded supported backends. 

"""
function device_wakeup!(handle::CuDevice)
    @info "waking $(handle) up"
    CUDA.device!(handle)
    return nothing
end

"""
    device_free(mesh::Mesh,::Val{:CUDA})

Description:
---
Return Dicts of effective cpu and gpu backend based on hard-coded supported backends. 
"""
function device_free!(mesh::Mesh,::Val{:CUDA})
    msg = ["freed GPU memory (VRAM):\n"]
    for i in 1:nfields(mesh)
        field = getfield(mesh,i)
        if typeof(field) <: AbstractArray
            CUDA.unsafe_free!(getfield(mesh,i))
        end
    end
    CUDA.reclaim()
    @info join(msg)
    mesh = nothing
    GC.gc()
    return nothing
end