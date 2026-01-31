export add_backend!, device_wakeup!, device_free!

function list_host_backend()
	return Dict( 
		:x86_64  => Dict(
            :host => "cpu",
            :Backend => CPU(),
            :brand => ["Intel(R)","AMD"],
            :wrapper => Array,
            :devices => nothing,
            :name => nothing,
            :handle => Val{:Host},
            :functional => Sys.ARCH==:x86_64,
        ),
		:aarch64 => Dict(
            :host => "cpu",
            :Backend => CPU(),
            :brand => ["Apple","AMD"],
            :wrapper => Array,
            :devices => nothing,
            :name => nothing,
            :handle => Val{:Host},
            :functional => Sys.ARCH==:aarch64,
        ),
	)
end

function list_cpu_devices()
    return [String(split(string(cpu), ":")[1]) for cpu in Sys.cpu_info()]
end

"""
    add_backend!(bckd::Execution, ::Val{:x86_64})

Description:
---
Populate Execution struct with effective cpu and gpu backend based on hard-coded supported backends. 
"""
function add_backend!(bckd::Execution,::Val{:x86_64})
	for (k,(platform,backend)) ∈ enumerate(list_host_backend())
		if backend[:functional]
            cpu_info = Sys.cpu_info()
            if !isempty(cpu_info) && !isempty(cpu_info[1].model)
				for brand ∈ backend[:brand]
					if occursin(brand,cpu_info[1].model)
                        bckd.cpu = Dict{Symbol,Any}()
                        for (k,device) ∈ enumerate(list_cpu_devices())
                            bckd.cpu[Symbol("dev$(k)")] = Dict(
                                :host     => "cpu",   
                                :platform => :CPU,        
                                :brand    => brand,            
                                :name     => cpu_info[1].model,
                                :Backend  => backend[:Backend],
                                :wrapper  => backend[:wrapper],
                                :handle   => nothing,
                            )      
                        end
						push!(bckd.functional,"✓ $(brand) $(platform)")
                        break
					end
				end
            else
                throw(ErrorException("Could not retrieve CPU model"))
            end
		end
	end
    @info join(bckd.functional,"\n")
	return nothing
end

"""
    get_host(bckd::Execution; prompt::Bool=false, distributed::Bool=false)

Description:
---
Return a NamedTuple of effective cpu device(s).
"""
function get_host(bckd::Execution; prompt::Bool=false, distributed::Bool=false)
    cpus,devs,names = Dict(),collect(keys(bckd.cpu)),Vector{String}()
    for key ∈ devs
        push!(names,bckd.cpu[key][:name])
    end
    if distributed
        for dev ∈ request("select device(s):",MultiSelectMenu(names))
            cpus[devs[dev]] = bckd.cpu[devs[dev]]
        end 
        return NamedTuple(cpus)
    elseif prompt
        dev = request("select device:",RadioMenu(names))
        return NamedTuple(Dict(devs[dev] => bckd.cpu[devs[dev]]))
    else
        return NamedTuple(Dict(:dev1 => bckd.cpu[:dev1]))
    end
end

"""
    get_device(bckd::Execution; prompt::Bool=false, distributed::Bool=false)

Description:
---
Return a NamedTuple of gpu(s). When `prompt=true`, show interactive menu for device selection.
"""
function get_device(bckd::Execution; prompt::Bool=false, distributed::Bool=false)
    devs,names = collect(keys(bckd.gpu)),Vector{String}()
    for key ∈ devs
        push!(names,bckd.gpu[key][:name])
    end
    gpus = Dict()
    if distributed
        for dev ∈ request("select device(s):",MultiSelectMenu(names))
            gpus[devs[dev]] = bckd.gpu[devs[dev]]
        end 
        return NamedTuple(gpus)
    elseif prompt
        dev = request("select device:",RadioMenu(names))
        return NamedTuple(Dict(devs[dev] => bckd.gpu[devs[dev]]))
    else
        return NamedTuple(Dict(devs[1]=>bckd.gpu[devs[1]]))
    end
end

"""
    device_wakeup!()

Description:
---
Return Dicts of effective cpu and gpu backend based on hard-coded supported backends. 

"""
function device_wakeup!()
    throw(ErrorException("🚧 `device_wakeup!()` is a stub. It must be overloaded in CUDAExt, ROCmExt or MtlExt."))
end

"""
    device_free(mesh::Mesh,::Val{:CPU})

Description:
---
Return Dicts of effective cpu and gpu backend based on hard-coded supported backends. 
"""
function device_free!(mesh::Mesh,::Val{:CPU})
    mesh = nothing
    GC.gc()
    return nothing
end

"""
    select_execution_backend(bckd::Execution, select::String="host"; prompt::Bool=false, distributed::Bool=false)

Description:
---
Select execution backend and return a NamedTuple of device(s).

Arguments:
- `select::String`: "host" for CPU, "device" for GPU
- `prompt::Bool`: Show interactive menu to choose specific device/core
- `distributed::Bool`: Enable distributed multi-device selection

Example:
---
```julia
julia> select_execution_backend(info.bckd, "host")  # Direct CPU selection
julia> select_execution_backend(info.bckd, "host", prompt=true)  # Choose which CPU core
julia> select_execution_backend(info.bckd, "host", distributed=true)  # Choose multiple CPU cores
```
"""
function select_execution_backend(bckd::Execution, select::String="host"; prompt::Bool=false, distributed::Bool=false)
    mode = distributed ? "distributed" : (prompt ? "interactive" : "default")
    if select == "host"
        @info "Using CPU backend ($mode mode)"
        return get_host(bckd; prompt=prompt, distributed=distributed)
    elseif select == "device"
        if isempty(bckd.gpu)
            @info "No GPU available, falling back to CPU backend"
            return get_host(bckd; prompt=prompt, distributed=distributed)
        else
            @info "Using GPU backend ($mode mode)"
            return get_device(bckd; prompt=prompt, distributed=distributed)
        end
    else
        throw(ArgumentError("Invalid backend selection: '$(select)'. Use 'host' or 'device'"))
    end
end