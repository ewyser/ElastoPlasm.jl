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
                push!(outputs, out)
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
function postprocess_fields(outputs, paths)
    @info "Postprocessing metanalysis results..."
    nsim = length(outputs)
    
    # Interactive field selection
    field_keys = sort(collect(keys(VARIABLE_EXTRACTORS)))
    field_names = [VARIABLE_EXTRACTORS[k].name for k in field_keys]
    menu = MultiSelectMenu(field_names, pagesize=length(field_names))
    choices = request("Select field(s) to analyze (Space to select, Enter to confirm):", menu)
    selected_fields = [field_keys[i] for i in choices]
    isempty(selected_fields) && error("No fields selected")
    @info "Selected fields: $(join([VARIABLE_EXTRACTORS[f].name for f in selected_fields], ", "))"
    
    # Load reference data
    reference = load(outputs[1].simulation)
    ref_mesh = reference["ic/mesh"]
    npts = size(reference["ic/mpts"].x, 2)
    
    # Process each selected field
    for selected_field in selected_fields
        field_info = VARIABLE_EXTRACTORS[selected_field]
        @info "Processing field: $(field_info.name)..."
        
        # Extract field data
        X, Y, D = extract_field_data(outputs, field_info, npts, nsim)
        
        # Compute statistics
        stats = compute_statistics(X, Y, D)
        
        # Generate and save plots
        plot_field_statistics(stats, field_info, reference, paths, selected_field, nsim)
    end
    return nothing
end

"""
    extract_field_data(outputs, field_info, npts, nsim) -> Tuple

Extract field data from all simulations.
"""
function extract_field_data(outputs, field_info, npts, nsim)
    X = Matrix{Float64}(undef, npts, nsim)
    Y = Matrix{Float64}(undef, npts, nsim)
    D = Matrix{Float64}(undef, npts, nsim)
    
    prog = make_progress(nsim; desc="Extracting $(field_info.name)")
    for (k, output) ∈ enumerate(outputs)
        jldopen(output.simulation,"r+") do file
            X[:, k] = vec(file["ic/mpts"].x[1, :])
            Y[:, k] = vec(file["ic/mpts"].x[2, :])
            D[:, k] = field_info.extract(file["ic/mpts"])
        end
        next!(prog; desc="Extracting $(field_info.name) $k/$nsim...")
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
function plot_field_statistics(stats, field_info, reference, paths, field_name, nsim)
    @info "Plotting averaged $(field_info.name) with standard deviation..."
    
    instr = reference["cfg/instr"]
    mesh = reference["ic/mesh"]
    ms = instr.plot.dpi * (mesh.prprt.L[1] / mesh.prprt.L[1]) / (mesh.prprt.nel[1] * 2)
    
    # Common plot settings
    common = (;
        xlim = (minimum(mesh.x[1, :]), maximum(mesh.x[1, :])),
        ylim = (minimum(mesh.x[2, :]), maximum(mesh.x[2, :])),
        size = instr.plot.dpi .* (mesh.prprt.L ./ mesh.prprt.L[1]),
    )
    
    # Define plot configurations with automatic color limits
    D_mean_min, D_mean_max = extrema(stats.D_mean)
    D_std_max = maximum(stats.D_std)
    
    plot_specs = [
        (data=stats.D_mean, color=:viridis, clim=(D_mean_min, D_mean_max),
         label=field_info.label, title="Average " * field_info.name * " " * L"\leftangle" * field_info.label * L"\rightangle_{n=%$nsim}"),
        (data=stats.D_std, color=:plasma, clim=(0, D_std_max),
         label=L"\sigma(" * field_info.label * L")", title="Standard deviation " * L"\sigma(" * field_info.label * L")_{n=%$nsim}"),
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
    save_plot((; file=joinpath(paths[:plot], "$(mesh.prprt.dim)d_av_$(field_name).png")))
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
    
    # 1. Preparation
    paths = set_paths(fid, info.sys.out; interactive=false)
    nsim = prepare_simulations!(L, nel, fid, paths)
    
    # 2. Calculation
    meta_path = joinpath(info.sys.out, fid)
    outputs = execute_simulations(meta_path, nsim)
    
    # 3. Postprocessing
    postprocess_fields(outputs, paths)
    return nothing
end