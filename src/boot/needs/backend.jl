# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
export add_backend!, device_wakeup!, device_free!
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
"""
    add_backend!(::Val{:x86_64},info::Self)

Description:
---
Return Dicts of effective cpu and gpu backend based on hard-coded supported backends. 
"""
function add_backend!(::Val{:x86_64},info::Self)
    String.([split(string(Sys.cpu_info()[k]),":")[1] for k ∈ 1:length(Sys.cpu_info())])
	arch = Sys.ARCH
	availables::Dict{Symbol, Dict{Symbol, Any}} = Dict( 
		:x86_64  => Dict(:host => "cpu",:Backend => CPU(),:brand => ["Intel(R)","AMD","Apple"],:wrapper => Array,:devices => nothing,:name => nothing,:handle => Val{:Host},:functional => arch==:x86_64 ,),
		:aarch64 => Dict(:host => "cpu",:Backend => CPU(),:brand => ["AMD","Apple"]           ,:wrapper => Array,:devices => nothing,:name => nothing,:handle => Val{:Host},:functional => arch==:aarch64,),
	)
	for (k,(platform,backend)) ∈ enumerate(availables)
		if backend[:functional]
            cpu_info = Sys.cpu_info()
            if !isempty(cpu_info) && !isempty(cpu_info[1].model)
				for brand ∈ backend[:brand]
					if occursin(brand,cpu_info[1].model)
                        cpu = Dict{Symbol,NamedTuple}()
                        cpu_devices = String.([split(string(Sys.cpu_info()[k]),":")[1] for k ∈ 1:length(Sys.cpu_info())])
                        for (k,device) ∈ enumerate(cpu_devices)
                            cpu[Symbol("dev$(k)")] = (;
                                host     = "cpu",   
                                platform = :CPU,        
                                brand    = brand,            
                                name     = cpu_info[1].model,
                                Backend  = backend[:Backend],
                                wrapper  = backend[:wrapper],
                                handle   = nothing,
                            )      
                        end
                        info.bckd.cpu[:dev0] = NamedTuple(cpu)
						push!(info.bckd.functional,"✓ $(brand) $(platform)")
					end
				end
            else
                throw(ErrorException("Could not retrieve CPU model"))
            end
		end
	end
    @info join(info.bckd.functional,"\n")
	return nothing
end

"""
    get_host(;user_select::Bool=true,mpi::Bool=info.mpi.is_active)

Description:
---
Return a NamedTuple of effective cpu. 
"""
function get_host(;user_select::Bool=true,mpi::Bool=info.mpi.is_active)
    cpus,devs,names = Dict(),collect(keys(info.bckd.cpu[:dev0])),Vector{String}()
    for key ∈ devs
        push!(names,info.bckd.cpu[:dev0][key][:name])
    end
    if mpi
        for dev ∈ request("Select device(s):",MultiSelectMenu(names))
            cpus[devs[dev]] = info.bckd.cpu[:dev0][devs[dev]]
        end 
        return NamedTuple(cpus)
    elseif info.ui.bckd && user_select
        for dev ∈ request("Select device(s):",RadioMenu(names))
            cpus[devs[dev]] = info.bckd.cpu[:dev0][devs[dev]]
        end 
        return NamedTuple(cpus)
    else
        return NamedTuple(Dict(:dev0 => info.bckd.cpu[:dev0][:dev1]))
    end
end

"""
    get_device(;user_select::Bool=true,mpi::Bool=info.mpi.is_active)

Description:
---
Return a NamedTuple of gpu(s) based on an interactive selection (when `user_select=true`) of available device(s) on the system. Otherwise, only return the first device found amongst all.
"""
function get_device(;user_select::Bool=true,mpi::Bool=info.mpi.is_active)
    devs,names = collect(keys(info.bckd.gpu)),Vector{String}()
    for key ∈ devs
        push!(names,info.bckd.gpu[key][:name])
    end
    gpus = Dict()
    if mpi
        for dev ∈ request("Select device(s):",MultiSelectMenu(names))
            gpus[devs[dev]] = info.bckd.gpu[devs[dev]]
        end 
        return NamedTuple(gpus)
    elseif info.ui.bckd && user_select
        for dev ∈ request("Select device(s):",RadioMenu(names))
            gpus[devs[dev]] = info.bckd.gpu[devs[dev]]
        end 
        return NamedTuple(gpus)
    else
        return NamedTuple(Dict(devs[1]=>info.bckd.gpu[devs[1]]))
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
    select_execution_backend(select::{Bool|String})

Description:
---
Call functions `get_host()` and/or `get_device()`, and return a `NamedTuple` of selected backends (if `select=true`) or of predefined backends (if `select={"cpu"|"gpu"}`). 

Example:
---
```julia
julia> cORIUm.select_execution_backend(true)
[ Info: unable to select execution device : defaulting host CPU
(dev0 = (host = "cpu", platform = :x86_64, brand = "Intel(R)", name = "Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz", Backend = KernelAbstractions.CPU(false), wrapper = Array, handle = nothing),)

julia> cORIUm.select_execution_backend("cpu")
[ Info: select CPU as execution device
(dev0 = (host = "cpu", platform = :x86_64, brand = "Intel(R)", name = "Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz", Backend = KernelAbstractions.CPU(false), wrapper = Array, handle = nothing),)

julia> 

```
"""
function select_execution_backend(select::Bool; user_select::Bool=true,mpi::Bool=info.mpi.is_active)
    if isempty(info.bckd.gpu)
        @info "Unable to select execution device : defaulting host CPU"
        return get_host(; user_select=select,mpi=mpi)
    else
        @info "Select execution device : (✓) host CPU and (✓) device GPU:" 
        backends = [:cpu,:gpu]
        select   = request("Please use up/down arrows to choose an option:",RadioMenu(String.(backends)))
        if backends[select] == :cpu
            return get_host(; user_select=user_select,mpi=mpi)
        else
            return get_device(; user_select=user_select,mpi=mpi)
        end
    end  
end
function select_execution_backend(select::String; user_select::Bool=true,mpi::Bool=info.mpi.is_active)
    if select == "host"
        @info "Select CPU as execution device"
        return get_host(; user_select=user_select,mpi=mpi)
    elseif select == "device"
        if isempty(info.bckd.gpu)
            @info "Unable to select execution device : defaulting host CPU"
            return get_host(; user_select=user_select,mpi=mpi)
        else
            @info "Select GPU as execution device"
            return get_device(;user_select=false)
        end
    else
        return throw("$(select) is not a valid backend. Use select = {host|device}")
    end
end