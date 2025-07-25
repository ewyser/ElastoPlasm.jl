export add_backend!, device_wakeup!, device_free!

"""
    add_backend!(::Val{:x86_64}, info::Self)

Detects and registers available CPU backends for x86_64 architecture, populating the `info.bckd` structure with supported devices and brands.

# Arguments
- `info::Self`: The backend info object to populate.

# Behavior
- Scans system CPU info and matches against supported brands.
- Populates `info.bckd.cpu` and updates `info.bckd.functional` with detected devices.
- Throws an error if CPU model cannot be retrieved.

# Returns
- `nothing`
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
    get_host(; user_select::Bool=true, mpi::Bool=info.mpi.is_active)

Returns a NamedTuple of selected CPU devices, optionally allowing interactive selection.

# Keyword Arguments
- `user_select::Bool=true`: If true, prompts user to select devices.
- `mpi::Bool=info.mpi.is_active`: If true, uses MPI selection logic.

# Returns
- `NamedTuple`: Selected CPU devices.
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
    get_device(; user_select::Bool=true, mpi::Bool=info.mpi.is_active)

Returns a NamedTuple of selected GPU devices, optionally allowing interactive selection.

# Keyword Arguments
- `user_select::Bool=true`: If true, prompts user to select devices.
- `mpi::Bool=info.mpi.is_active`: If true, uses MPI selection logic.

# Returns
- `NamedTuple`: Selected GPU devices.
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

Stub for device backend initialization. Must be overloaded in platform-specific extensions.

# Throws
- `ErrorException`: Always, indicating the function is a stub.
"""
function device_wakeup!()
    throw(ErrorException("🚧 `device_wakeup!()` is a stub. It must be overloaded in CUDAExt, ROCmExt or MtlExt."))
end

"""
    device_free!(mesh::Mesh, ::Val{:CPU})

Frees resources associated with a CPU mesh and triggers garbage collection.

# Arguments
- `mesh::Mesh`: The mesh object to free.

# Returns
- `nothing`
"""
function device_free!(mesh::Mesh,::Val{:CPU})
    mesh = nothing
    GC.gc()
    return nothing
end

"""
    select_execution_backend(select::Bool; user_select::Bool=true, mpi::Bool=info.mpi.is_active)

Interactively selects and returns a NamedTuple of CPU or GPU backends based on user input or available devices.

# Arguments
- `select::Bool`: If true, prompts user for device selection.

# Keyword Arguments
- `user_select::Bool=true`: If true, enables interactive selection.
- `mpi::Bool=info.mpi.is_active`: If true, uses MPI selection logic.

# Returns
- `NamedTuple`: Selected backend devices.
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

"""
    select_execution_backend(select::String; user_select::Bool=true, mpi::Bool=info.mpi.is_active)

Selects and returns a NamedTuple of CPU or GPU backends based on the provided string ("host" or "device").

# Arguments
- `select::String`: "host" for CPU, "device" for GPU.

# Keyword Arguments
- `user_select::Bool=true`: If true, enables interactive selection.
- `mpi::Bool=info.mpi.is_active`: If true, uses MPI selection logic.

# Returns
- `NamedTuple`: Selected backend devices.

# Throws
- Error if an invalid backend string is provided.
"""
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