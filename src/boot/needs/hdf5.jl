function superstore(mpi::NamedTuple,startup::String; deflate::Int=3)
    config = load(startup)["config"]
    instr  = config.instr
    checks = config.launch.chck

    Nx,Ny = Int(mpi.globdim[1]),Int(mpi.globdim[2])
    ix,iy = mpi.index[:i],mpi.index[:j]
    nx,ny = length(ix),length(iy)

    # set path to h5file
    path = joinpath(info.sys.out,"$(join(instr[:save][:name][1:end-1])).h5")
    # set args for h5open
    if mpi.status
        h5args = ("cw",mpi.comm,mpi.info)
    else
        h5args = ("cw",)
    end
    # open & write to h5
    h5open(load(startup)["h5file"], h5args...) do file
        # create test structure for mpi
        ranks       = create_group(file,"mpi-rank")
        dims        = (Nx,Ny)
        dset        = create_dataset(ranks,"data",datatype(Int32),dataspace(dims))
        dset[ix,iy] = Int32(mpi.me+1).*ones(Int32,nx,ny)
        # create B-tree set structure
        set  = create_group(file,"set")   
        all  = unique(vcat(["z"],instr[:save][:what]))
        excl = vcat("z",filter(x -> occursin("max", x), instr[:save][:what]))
        for field ∈ all
            GRP = create_group(set,"$(field)")
            if !occursin("max", field)
                for k ∈ 1:length(checks)
                    name = "$(field)$(lpad(k,round(Int,1+log10(length(checks))),"0"))"
                    grp = create_group(GRP,name)
                end
            end
        end
    end
    return load(startup)["h5file"]::String
end
function superstore(ic,instr,checks; deflate::Int=3)
    nx,ny = ic.dim
    Δx,Δy = ic.res
    # set path to h5file
    path = joinpath(info.sys.out,"$(join(instr[:save][:name][1:end-1])).h5")
    # open & write to h5
    h5open(path,"cw") do file
        # create B-tree set structure
        set  = create_group(file,"set")   
        all  = unique(vcat(["z"],instr[:save][:what]))
        excl = vcat("z",filter(x -> occursin("max", x), instr[:save][:what]))
        for field ∈ all
            GRP = create_group(set,"$(field)")
            attrs(GRP)["dim"] = ([nx,ny])
            attrs(GRP)["res"] = ([Δx,Δy])
            if field ∈ excl
                if field == "z"
                    grp  = create_group(GRP,"$field")
                    write(grp,"data",ic[Symbol(field)])
                end
            else
                for k ∈ 1:length(checks)
                    name = "$(field)$(lpad(k,round(Int,1+log10(length(checks))),"0"))"
                    grp = create_group(GRP,name)
                end
            end
        end
    end
    return path::String
end
"""
    dumping(meshD::NamedTuple, t::Real, ns::Tuple{Any})

"""
function dumping(mpi::NamedTuple,startup::String,mesh::Mesh,t::Real,k::Number; deflate::Int=3)
    # unpack
    h5file = load(startup)["h5file"]
    # unpack & cache
    Nx,Ny = Int(mpi.globdim[1]),Int(mpi.globdim[2])
    ix,iy = mpi.index[:i],mpi.index[:j]
    nx,ny = length(ix),length(iy)
    # set args for h5open
    if mpi.status
        h5args,chunk_size = ("r+", mpi.comm, mpi.info),(nx, ny)
    else
        h5args,chunk_size = ("r+",),(1, 1)
    end
    # open & write to h5
    h5open(h5file, h5args...) do file
        excluded = vcat("z",filter(x -> occursin("max", x), keys(file["set"])))
        for key ∈ filter!(x->x∉excluded,keys(file["set"])) 
            GRP  = file["set"][key]   
            grp  = GRP[keys(GRP)[k]]
            data = Array(getfield(mesh,Symbol(key)))
            if size(data,1) != mesh.nx && size(data,2) != mesh.ny
                dims          = (Nx,Ny,3)
                dset          = create_dataset(grp,"data",datatype(Float32),dataspace(dims))
                dset[ix,iy,:] = Float32.(data[2:end-1,2:end-1,:])
            else
                dims        = (Nx,Ny)
                dset        = create_dataset(grp,"data",datatype(Float32),dataspace(dims))
                dset[ix,iy] = Float32.(data)
            end
            attrs(grp)["elapsed"] = t 
        end
    end
    return nothing
end

function dumping(mpi::NamedTuple,startup::String,mesh::Mesh,t::Real,::Val{:final}; deflate::Int=3)
    # unpack
    h5file = load(startup)["h5file"]
    # unpack & cache
    Nx,Ny = Int(mpi.globdim[1]),Int(mpi.globdim[2])
    ix,iy = mpi.index[:i],mpi.index[:j]
    nx,ny = length(ix),length(iy)
    # set args for h5open
    if mpi.status
        args,chunk_size = ("r+", mpi.comm, mpi.info),(nx, ny)
    else
        args,chunk_size = ("r+",),(1, 1)
    end
    h5open(h5file,args...) do file
        for (k,key) ∈ enumerate(filter(x -> occursin("max", x), keys(file["set"])))
            grp  = file["set"]["$(key)"]
            data = Array(getfield(mesh,Symbol(key)))
            if key == "Umax"
                dims          = (Nx,Ny,3)
                dset          = create_dataset(grp,"data",datatype(Float32),dataspace(dims))
                dset[ix,iy,:] = Float32.(data)
            else
                dims          = (Nx,Ny)
                dset          = create_dataset(grp,"data",datatype(Float32),dataspace(dims))
                dset[ix,iy]   = Float32.(data)
            end  
            attrs(grp)["elapsed"] = t 
        end
    end
    return nothing
end

function dumping!(h5src::String,instr)
    out = nothing
    h5open(h5src,"r") do file
        dim = attrs(file["set"]["z"])["dim"]
        res = attrs(file["set"]["z"])["res"]
        key = attrs(file["set"]     )["list"]
        key = key[1:end-1]
        out = Dict(
            :dim   => dim,
            :res   => res,
            :xc    => read(file["where"]["coordinate"]["xc"]),
            :yc    => read(file["where"]["coordinate"]["yc"]),
            :z     => read(file["set"]["z"]["data"]),
        )
        t = []
        for (k,field) ∈ enumerate(vcat(filter(x -> !occursin("max", x), instr[:save].what)))
            FID = file["set"][field]
            serie = keys(FID)
            temp  = []
            for (l,key) ∈ enumerate(serie)
                data = read(FID[key]["data"])
                if l == 1
                    temp = data
                else
                    temp = cat(temp,data,dims=length(size(data))+1)
                end
                push!(t,attrs(FID[key])["elapsed"])
            end
            out[Symbol(field)] = temp
        end
        out[:t] = sort(unique(t))
    end
    return out
end

function solidstatesync(ic::NamedTuple,globloc::Val{:global}; deflate::Int=3, path::String = joinpath(info.sys.out,"ic.h5"))
    # log message
    @info """solidstatesync from root proc: 
    - storing $(typeof(globloc).parameters[1]) initial conditions into $(basename(path))
    """
    # existence check    
    if isfile(path)
        rm(path, force=true)
    end
    # open & write to h5
    h5open(path,"cw") do file
        ic_group  = create_group(file,"ic")   
        for field ∈ collect(keys(ic))
            if isa(ic[field],AbstractArray)
                write(ic_group,string(field),ic[field])
            end
        end
    end
    return path::String
end
function solidstatesync(mpi::NamedTuple,startup::String,mesh::Mesh{T1,T2,D},globloc::Val{:local}) where {T1,T2,D}
    # log message
    @info """solidstatesync from ($(mpi.me)/$(mpi.comm_size)): 
    - accumulating $(typeof(globloc).parameters[1]) initial conditions into $(basename(load(startup)["ic"]))
    """
    # init  
    ix,iy = mpi.index[:i],mpi.index[:j]
    nx,ny = length(ix),length(iy)
    # set args for h5open
    if mpi.status
        args,chunk_size = ("r+", mpi.comm, mpi.info),(nx, ny)
    else
        args,chunk_size = ("r+",),(1, 1)
    end
    h5open(load(startup)["ic"],args...) do file
        ic = file["ic"]
        # host ------------------------------------------------------------------------------------------------------------------------------------------  
        ic["t"][:]        = mesh.T
        # device ----------------------------------------------------------------------------------------------------------------------------------------
        ic["z"][ix,iy]    = mesh.z[2:end-1,2:end-1]

        ic["h"][ix,iy]    = mesh.U[2:end-1,2:end-1,1]
        ic["qx"][ix,iy]   = mesh.U[2:end-1,2:end-1,2]
        ic["qy"][ix,iy]   = mesh.U[2:end-1,2:end-1,3]

        ic["H"][ix,iy]    = mesh.Umax[:,:,1]
        ic["Qx"][ix,iy]   = mesh.Umax[:,:,2]
        ic["Qy"][ix,iy]   = mesh.Umax[:,:,3]
        ic["Qmax"][ix,iy] = mesh.Qmax
        ic["tmax"][ix,iy] = mesh.tmax
        ic["n"][:]        = mesh.n
        ic["K"][:]        = mesh.K
        ic["p"][ix,iy]    = mesh.p
        ic["F"][ix,iy]    = mesh.dF
        ic["St"][ix,iy]   = mesh.St
    end
    return load(startup)["ic"]::String
end

function convert_h5_to_nt(path::String)
    temp = Dict{Symbol,Any}()
    h5open(path,"r") do file
        ic = file["ic"]
        for field ∈ keys(ic)
            #println(field)
            temp[Symbol(field)] = read(ic[field])
        end
    end
    return ic = (; temp...) 
end