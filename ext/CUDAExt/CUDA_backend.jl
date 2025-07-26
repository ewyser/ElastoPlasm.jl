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
#=

=#





















#=
function available_backends(::Val{:ROCm}; msg::Vector{String} = ["functional execution platform(s):"], gpus::Dict{Symbol,NamedTuple} = Dict{Symbol,NamedTuple}())
    availables::Dict{Symbol, Dict{Symbol, Any}} = Dict( 
        :ROCm => Dict(
            :dev  => "gpu",
            :Backend => ROCBackend(),
            :brand => "AMD",
            :wrapper => ROCArray,
            :devices => AMDGPU.devices,
            :name => AMDGPU.HIP.name,
            :handle => HIPDevice ,
            :functional => AMDGPU.functional(),
        ), 

    )
    dev_id = length(keys(cORIUm.getinfobckd.gpu))
    gpus = Dict{Symbol, NamedTuple}()
    for (platform, backend) in availables
        if haskey(backend, :dev)
            if backend[:functional]
                for (k, device) in enumerate(backend[:devices]())
                    dev_id += 1
                    gpu = (; 
                        dev      = backend[:dev],
                        platform = platform,
                        brand    = backend[:brand],
                        name     = backend[:name](device),
                        Backend  = backend[:Backend],
                        wrapper  = backend[:wrapper],
                        handle   = device,
                    )
                    gpus[Symbol("dev$(dev_id)")] = gpu
                end
                push!(msg, "✓ $(backend[:brand]) $(platform)")
            else
                push!(msg, "✗ $(backend[:brand]) $(platform)")
            end
        end
    end
    return gpus, msg
end
function available_backends(::Val{:Mtl}; msg::Vector{String} = ["functional execution platform(s):"], gpus::Dict{Symbol,NamedTuple} = Dict{Symbol,NamedTuple}())
    availables::Dict{Symbol, Dict{Symbol, Any}} = Dict( 
        #=:Metal => Dict(
            :dev  => "gpu",
            :Backend =>              ,
            :brand => "Apple"           ,
            :wrapper => MtlArray,
            :devices => Metal.MTL.devices,
            :name =>                ,
            :handle =>           ,
            :functional =>                    ,
        ),=# 

    )
    dev_id = length(keys(cORIUm.getinfobckd.gpu))
    gpus = Dict{Symbol, NamedTuple}()
    for (platform, backend) in availables
        if haskey(backend, :dev)
            if backend[:functional]
                for (k, device) in enumerate(backend[:devices]())
                    dev_id += 1
                    gpu = (; 
                        dev      = backend[:dev],
                        platform = platform,
                        brand    = backend[:brand],
                        name     = backend[:name](device),
                        Backend  = backend[:Backend],
                        wrapper  = backend[:wrapper],
                        handle   = device,
                    )
                    gpus[Symbol("dev$(dev_id)")] = gpu
                end
                push!(msg, "✓ $(backend[:brand]) $(platform)")
            else
                push!(msg, "✗ $(backend[:brand]) $(platform)")
            end
        end
    end
    return gpus, msg
end

function device_free(mesh::Mesh,bckd::NamedTuple)
    if haskey(bckd,:dev) && bckd.platform == :CUDA
        msg = ["freed GPU memory (VRAM):\n"]
        for i in 1:nfields(mesh)
            field = getfield(mesh,i)
            if typeof(field) <: AbstractArray
                CUDA.unsafe_free!(getfield(mesh,i))
            end
        end
        CUDA.reclaim()
        @info join(msg)
    elseif haskey(bckd,:dev) && bckd.platform == :ROCm
        msg = ["freed GPU memory (VRAM):\n"]
        #=map = mapsto(mesh.dim[1],mesh.dim[2])
        for (i,field) ∈ enumerate(collect(keys(map)))
            if map[field][:exec] == :device
                if field == :Dn
                    for (j,key) ∈ enumerate([:i,:j,:Q])
                        AMDGPU.unsafe_free!(getfield(mesh.D,key))
                    end
                    push!(msg,"(✓) $(field), ")
                else
                    AMDGPU.unsafe_free!(getfield(mesh,field))
                    push!(msg,"(✓) $(field), ")
                end
            end
        end
        =#
        @info join(msg)
    end
    mesh = nothing
    GC.gc()
    return nothing
end


function device_wakeup!(handle::HIPDevice)
    @info "waking $(handle) up"
    AMDGPU.device!(handle)
    return nothing
end

=#