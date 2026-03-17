# Variable extraction configuration function with default colorbar limits
export get_mpts_variable_config, get_variable_plot_options

# Variable extraction functions for metanalysis
get_epII(mpts) = @views vec(mpts.s.ϵpII[1, :])
get_P(mpts)    = @views -vec(mean(mpts.s.σᵢ, dims=1)) / 1e3
get_J(mpts)    = @views vec(mpts.J)
get_v(mpts)    = @views vec(sqrt.(mpts.s.v[1, :].^2 .+ mpts.s.v[2, :].^2)) 
get_Δu(mpts)   = @views vec(sqrt.(mpts.s.u[1, :].^2 .+ mpts.s.u[2, :].^2)) 
get_coh0(mpts) = @views vec(mpts.s.c₀./1e3)
get_phi0(mpts) = @views vec(mpts.s.ϕ₀.*180.0/π)

"""
    get_mpts_variable_config() -> Dict{String, NamedTuple}

Return a dictionary of material point variable plotting configurations, including extraction functions, labels, units, colorbar settings, and limits.
"""
function get_mpts_variable_config()
    return Dict(
        "epII" => (
            data=get_epII, 
            label=L"\epsilon_{\mathrm{II}}^{\mathrm{acc}}",  
            unit = "[-]",
            name="Plastic strain",
            cb=:viridis,
            cblim=(0.0, 1.5)
        ),
        "P" => (
            data=get_P,    
            label=L"p",                                     
            unit = "[kPa]",
            name="Pressure",
            cb=:viridis,
            cblim=(0.0, 150)
        ),
        "J" => (
            data=get_J,    
            label=L"J",                                     
            unit = "[-]",
            name="Deformation determinant",
            cb=:viridis,
            cblim=(0.5, 1.5)
        ),
        "Δu" => (
            data=get_Δu,
            label=L"\Delta u",  
            unit = "[m]",
            name="Displacement magnitude",
            cb=:viridis,
            cblim=(0.0, 5.0)
        ),
        "v" => (
            data=get_v,
            label=L"v",  
            unit = "[m/s]",
            name="Velocity magnitude",
            cb=:viridis,
            cblim=(0.0, 5.0)
        ),    
        "coh0" => (
            data=get_coh0,
            label=L"c_{0}",  
            unit = "[kPa]",
            name="Initial cohesion",
            cb = :vik,
            cblim=(10.0, 30.0)
        ),
        "phi0" => (
            data=get_phi0,
            label=L"\phi_{0}",  
            unit = "[deg.]",
            name="Initial friction angle",
            cb = :viridis,
            cblim=(20.0, 40.0)
        )
    )
end

"""
    get_mesh_variable_config() -> Dict{String, NamedTuple}

Return a dictionary of mesh variable plotting configurations, including extraction functions, labels, units, colorbar settings, and limits.
"""
function get_mesh_variable_config()
    return Dict(
        "m" => (
            data=get_epII, 
            label=L"m_{s}",  
            unit = "[kg]",
            name="Solid mass",
            cb=:viridis,
            cblim=(100.0, 1000.0)
        ),
    )
end

"""
    get_variable_plot_options() -> Vector{NamedTuple}

Generate kwargser-compatible plot options from both material point and mesh variable configurations.

# Returns
- Vector of named tuples with all fields from get_mpts_variable_config() and get_mesh_variable_config(), plus the variable key.

# Example
opts = get_variable_plot_options()
# Returns [(;mpts=(key="epII", data=..., label=..., cb=..., cblim=...)), (;mesh=(key="m", ...)), ...]
"""
function get_variable_plot_options()
    mpts = [(;mpts=merge((;key=k), v)) for (k, v) ∈ get_mpts_variable_config()]
    mesh = [(;mesh=merge((;key=k), v)) for (k, v) ∈ get_mesh_variable_config()]
    return vcat(mpts,mesh)
end