using Suppressor,Statistics

export metanalysis

"""
    prepare_simulations!(L, nel, fid, paths) -> Int

Prepare mesh and material points for all simulation configurations.

# Returns
- `nsim`: Total number of simulations prepared
"""
function prepare_simulations!(L, nel, fid, paths)
    deforms, trfrs = ["finite", "infinitesimal"], ["std", "tpic", "apic"]
    sim, nsim = 1, length(deforms) * length(trfrs)
    
    @info "Running metanalysis for $(length(L))d slump problem"
    prog = make_progress(nsim; desc="Preparing simulation(s)")
    for (i, deform) ∈ enumerate(deforms)
        for (j, trsfr) ∈ enumerate(trfrs)
            fwrk = (;
                deform = deform,
                trsfr = trsfr,
                C_pf = 1.0,
                musl = true,
                locking = true,
                damping = 0.1
            )
            plot = (;
                status = false,
                freq = 1.0,
                dpi = 500,
                what = [(; mpts=(name="epII", cblim=(0.0, 1.5)),),],
            )
            @suppress begin
                _ = ic_slump(L, nel; fid="$fid/sim_$sim", fwrk=fwrk, plot=plot)
            end
            next!(prog; desc="Preparing simulation(s) $sim/$nsim...")
            sim += 1
        end
    end
    return nsim
end

"""
    execute_simulations(meta_path, nsim) -> Vector

Execute all prepared simulations and return outputs.

# Returns
- `outputs`: Vector of simulation output objects
"""
function execute_simulations(meta_path, nsim)
    batch = filter(d -> isdir(joinpath(meta_path, d)), readdir(meta_path))
    isempty(batch) && error("No simulation directories found in $meta_path")
    @info "Available simulations for metanalysis in $(trunc_path(meta_path)):"
    
    outputs = []
    prog = make_progress(nsim; desc="Executing simulation(s)")
    for (k, sim) ∈ enumerate(batch)
        sim_path = joinpath(meta_path, sim)
        for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
            @suppress begin
                out = elastoplasm!(joinpath(sim_path, file); workflow=[elastodynamic!, elastoplastic!])
                push!(outputs, out.simulation)
            end
        end
        next!(prog; desc="Executing simulation(s) $k/$nsim...")
    end
    return outputs
end

"""
    postprocess_fields(outputs, paths) -> Nothing

Extract fields, compute statistics, and generate plots for selected fields.
"""
function postprocess_fields(outputs, path)
    @info "Postprocessing metanalysis results..."
    nsim = length(outputs)
    
    # Interactive field selection
    field_keys  = sort(collect(keys(get_mpts_variable_config())))
    field_names = [get_mpts_variable_config()[k].name for k in field_keys]
    menu        = MultiSelectMenu(field_names, pagesize=length(field_names))
    choices     = request("Select field(s) to analyze (Space to select, Enter to confirm):", menu)
    selection   = [field_keys[i] for i in choices]
    isempty(selection) && error("No fields selected")

    # Load reference data
    reference = load(outputs[1])
    npts      = size(reference["ic/mpts"].x, 2)
    
    # Process each selected field
    for field in selection
        get = get_mpts_variable_config()[field]
        @info "Processing field: $(get.name)..."
        
        # Extract field data
        X, Y, D = extract_field_data(outputs, get, npts, nsim)
        
        # Compute statistics
        stats = compute_statistics(X, Y, D)
        
        # Generate and save plots
        plot_field_statistics(stats, get, reference, path, field, nsim)
    end
    return nothing
end

"""
    extract_field_data(outputs, field_info, npts, nsim) -> Tuple

Extract field data from all simulations.
"""
function extract_field_data(outputs, get, npts, nsim)
    X = Matrix{Float64}(undef, npts, nsim)
    Y = Matrix{Float64}(undef, npts, nsim)
    D = Matrix{Float64}(undef, npts, nsim)
    
    prog = make_progress(nsim; desc="Extracting $(get.name)")
    for (k, output) ∈ enumerate(outputs)
        jldopen(output,"r+") do file
            X[:, k] = vec(file["ic/mpts"].x[1, :])
            Y[:, k] = vec(file["ic/mpts"].x[2, :])
            D[:, k] = get.data(file["ic/mpts"])
        end
        next!(prog; desc="Extracting $(get.name) $k/$nsim...")
    end
    finish!(prog)
    
    return X, Y, D
end

"""
    compute_statistics(X, Y, D) -> NamedTuple

Compute mean and standard deviation statistics.
"""
function compute_statistics(X, Y, D)
    return (;
        X_mean = vec(mean(X, dims=2)),
        Y_mean = vec(mean(Y, dims=2)),
        D_mean = vec(mean(D, dims=2)),
        D_std = vec(std(D, dims=2)),
    )
end

"""
    plot_field_statistics(stats, field_info, reference, paths, field_name, nsim) -> Nothing

Generate and save plots for field statistics.
"""
function plot_field_statistics(stats, field_info, reference, path, field_name, nsim)
    @info "Plotting averaged $(field_self.name) with standard deviation..."
    
    instr = reference["cfg/instr"]
    mesh = reference["ic/mesh"]
    ms = instr.plot.dpi * (mesh.prprt.L[1] / mesh.prprt.L[1]) / (mesh.prprt.nel[1] * 2)
    
    # Common plot settings
    common = (;
        xlim = (minimum(getindex.(mesh.x, 1)), maximum(getindex.(mesh.x, 1))),
        ylim = (minimum(getindex.(mesh.x, 2)), maximum(getindex.(mesh.x, 2))),
        size = instr.plot.dpi .* (mesh.prprt.L ./ mesh.prprt.L[1]),
    )
    
    # Define plot configurations with automatic color limits
    D_mean_min, D_mean_max = extrema(stats.D_mean)
    D_std_max = maximum(stats.D_std)
    
    plot_specs = [
        (data=stats.D_mean, color=:viridis, clim=(D_mean_min, D_mean_max),
         label=field_self.label, title="Average $(lowercase(field_self.name)) " * L"\leftangle" * field_self.label * L"\rightangle_{n=%$nsim}"),
        (data=stats.D_std, color=:plasma, clim=(0, D_std_max),
         label=L"\sigma(" * field_self.label * L")", title="Standard deviation " * L"\sigma(" * field_self.label * L")"),
    ]
    
    # Generate plots
    config_plot()
    gr(legend=true, markersize=ms, markershape=:circle, markerstrokewidth=0.75)
    plots = map(plot_specs) do spec
        plot(stats.X_mean, stats.Y_mean;
            seriestype = :scatter,
            marker_z = spec.data,
            xlabel = L"$x-$direction" * " [m]",
            ylabel = L"$z-$direction" * " [m]",
            label = spec.label * " [-]",
            color = spec.color,
            clim = spec.clim,
            xlim = common.xlim,
            ylim = common.ylim,
            title = spec.title,
            aspect_ratio = 1,
            size = common.size
        )
    end
    
    fig = display(plot(plots...; layout=(length(plots), 1), size=(common.size[1], common.size[2] * length(plots))))
    save_plot((; file=joinpath(path, "$((length(mesh.prprt.nel)-1))d_av_$(lowercase(field_name)).png")))
end

"""
    metanalysis(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., field::String="epII") -> Nothing

Perform metanalysis across multiple slump simulations with varying deformation frameworks and transfer schemes.

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier (default: script filename).
- `field::String`: (Optional) Field to average and plot (default: "epII").

# Returns
- `Nothing`: Saves averaged field plot to disk.

# Example
```julia
metanalysis([64.0, 16.0], [40, 10]; fid="slump_meta")
```

# Notes
- Tests all combinations of deformation frameworks (finite/infinitesimal) and transfer schemes (std/tpic/apic)
- Computes and plots ensemble average of the specified field across all simulations
"""
function metanalysis(L, nel; fid::String=first(splitext(basename(@__FILE__))), field::String="epII")
    @assert length(L) == length(nel) "L and nel must have same length"
    root = mkpath(joinpath(self.sys.out, fid))
    # Prepare and execute simulations if not already done
    if isempty(readdir(root))
        # 1. Preparation
        nsim = prepare_simulations!(L, nel, fid, root)
        # 2a. Calculation
        meta_path = joinpath(self.sys.out, fid)
        outputs = execute_simulations(meta_path, nsim)
    else
        # 2b. Load existing outputs
        outputs = String[]
        for sim ∈ filter(d -> isdir(joinpath(root, d)), readdir(root))
            sim_path = joinpath(root, sim)
            for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
                push!(outputs, joinpath(sim_path, file))
            end
        end
    end

    # 3. Postprocess results
    menu, msg = RadioMenu(["Yes", "No"], pagesize=2), "Do metanalysis ?"
    while [true, false][request(msg, menu)]
        postprocess_fields(outputs, root)
        msg  = "Redo metanalysis ?"
    end

    return nothing
end
function metanalysis(root::String, field::String="epII")
    outputs = String[]
    for sim ∈ filter(d -> isdir(joinpath(root, d)), readdir(root))
        sim_path = joinpath(root, sim)
        for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
            push!(outputs, joinpath(sim_path, file))
        end
    end
    menu, msg = RadioMenu(["Yes", "No"], pagesize=2), "Do metanalysis ?"
    while [true, false][request(msg, menu)]
        postprocess_fields(outputs, root)
        msg  = "Redo metanalysis ?"
    end

    return nothing
end