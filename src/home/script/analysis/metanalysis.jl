using Suppressor,Statistics

export metanalysis

"""
    prepare_simulations!(L, nel, fid, paths) -> Int

Prepare mesh and material points for all simulation configurations.

# Returns
- `nsim`: Total number of simulations prepared
"""
function prepare_simulations!(L, nel, fid, paths)
    basis, deforms, trfrs = ["mlsmpm","bsmpm","gimpm","smpm"], ["finite", "infinitesimal"], ["std", "tpic", "apic"]
    sim, nsim = 1, length(basis) * length(deforms) * length(trfrs)
    
    @info "Running metanalysis for $(length(L))d slump problem"
    prog = make_progress(nsim; desc="Preparing simulation(s)")
    
    for (h, shp) ∈ enumerate(basis)
        for (i, deform) ∈ enumerate(deforms)
            for (j, trsfr) ∈ enumerate(trfrs)
                basis = (;
                    which = shp,
                    how = nothing,
                )
                strain = (;
                    deform = deform,
                )
                transfer = (;
                    trsfr = trsfr,
                    C_pf = 1.0,
                    musl = true,
                )
                stab = (;
                    locking = false,
                    damping = 0.1
                )
                plot = (;
                    status = false,
                    freq = 1.0,
                    dpi = 500,
                    what = [(; mpts=(name="epII", cblim=(0.0, 1.5)),),],
                )
                @suppress begin
                    _ = slump_problem(L, nel; fid="$fid/sim_$sim", strain=strain, transfer=transfer, stab=stab, plot=plot)
                end
                next!(prog; desc="Preparing simulation(s) $sim/$nsim...")
                sim += 1
            end
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
                out = elastoplasm!(joinpath(sim_path, file); workflows=[elastodynamic!, elastoplastic!])
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
    npts      = length(reference["ic/problem"].mpts.x)
    
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
        jldopen(output,"r") do file
            mpts    = file["ic/problem"].mpts
            X[:, k] = getindex.(mpts.x, 1)
            Y[:, k] = getindex.(mpts.x, 2)
            D[:, k] = get.data(mpts) .* get.scale
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
    @info "Plotting averaged $(field_info.name) with standard deviation..."
    
    solver = reference["cfg/solver"]
    mesh = reference["ic/problem"].mesh
    mpts = reference["ic/problem"].mpts

    # Common plot settings
    common = (;
        xlim = (minimum(getindex.(mesh.x, 1)), maximum(getindex.(mesh.x, 1))),
        ylim = (minimum(getindex.(mesh.x, 2)), maximum(getindex.(mesh.x, 2))),
        size = solver.plot.dpi .* (mesh.prprt.L ./ mesh.prprt.L[1]),
    )

    # Marker size: px-per-physical-unit (limiting axis, aspect-ratio-aware) times mean mp spacing
    Δx, Δy  = common.xlim[2]-common.xlim[1], common.ylim[2]-common.ylim[1]
    ppu     = min(common.size[1]/Δx, common.size[2]/Δy)
    spacing = 2.0*sum(ℓ -> ℓ[1], mpts.ℓ₀)/length(mpts.ℓ₀)
    ms      = ppu*spacing

    # Define plot configurations with automatic color limits
    D_mean_min, D_mean_max = extrema(stats.D_mean)
    D_std_max = maximum(stats.D_std)
    
    plot_specs = [
        (data=stats.D_mean, color=:viridis, clim=(D_mean_min, D_mean_max),
         label=field_info.label * " in " * field_info.unit, title="Average $(lowercase(field_info.name)) " * L"\leftangle" * field_info.label * L"\rightangle_{n=%$nsim}"),
        (data=stats.D_std, color=:plasma, clim=(0, D_std_max),
         label=L"\sigma(" * field_info.label * L")" * " in " * field_info.unit, title="Standard deviation " * L"\sigma(" * field_info.label * L")"),
    ]
    
    # Generate plots
    config_plot()
    gr(legend=true, markershape=:circle, markerstrokewidth=0.75)
    plots = map(plot_specs) do spec
        plot(stats.X_mean, stats.Y_mean;
            seriestype = :scatter,
            marker_z = spec.data,
            markersize = ms,
            xlabel = L"$x-$direction" * " [m]",
            ylabel = L"$z-$direction" * " [m]",
            label = spec.label,
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
    collect_sim_outputs(root::String) -> Vector{String}

Collect the `.jld2` output paths for every simulation directory already present
under `root`.
"""
function collect_sim_outputs(root::String)
    outputs = String[]
    for sim ∈ filter(d -> isdir(joinpath(root, d)), readdir(root))
        sim_path = joinpath(root, sim)
        for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
            push!(outputs, joinpath(sim_path, file))
        end
    end
    return outputs
end

"""
    run_postprocess_loop(outputs, root) -> Nothing

Repeatedly prompt "Do/Redo metanalysis ?" and call `postprocess_fields` until the
user declines.
"""
function run_postprocess_loop(outputs, root)
    menu, msg = RadioMenu(["Yes", "No"], pagesize=2), "Do metanalysis ?"
    while [true, false][request(msg, menu)]
        postprocess_fields(outputs, root)
        msg  = "Redo metanalysis ?"
    end
    return nothing
end

"""
    metanalysis(L::Vector{Float64}, nel::Vector{Int64}; fid::String=...) -> Nothing

Perform metanalysis across multiple slump simulations with varying deformation frameworks and transfer schemes.

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier (default: script filename).

# Returns
- `Nothing`: Saves averaged field plot to disk.

# Example
```julia
metanalysis([64.0, 16.0], [40, 10]; fid="slump_meta")
```

# Notes
- Tests all combinations of deformation frameworks (finite/infinitesimal) and transfer schemes (std/tpic/apic)
- Computes and plots ensemble average of the field(s) selected interactively in `postprocess_fields`
"""
function metanalysis(L, nel; fid::String=first(splitext(basename(@__FILE__))))
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
        outputs = collect_sim_outputs(root)
    end

    # 3. Postprocess results
    run_postprocess_loop(outputs, root)

    return nothing
end
function metanalysis(root::String)
    outputs = collect_sim_outputs(root)
    run_postprocess_loop(outputs, root)

    return nothing
end