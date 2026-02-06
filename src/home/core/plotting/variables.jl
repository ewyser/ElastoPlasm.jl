# Variable extraction configuration dictionary with default colorbar limits
export MPTS_VAR, get_variable_plot_options

# Variable extraction functions for metanalysis
get_epII(mpts) = @views vec(mpts.s.ϵpII[1, :])
get_P(mpts)    = @views -vec(mean(mpts.s.σᵢ, dims=1)) / 1e3
get_J(mpts)    = @views vec(mpts.J)
get_v(mpts)    = @views vec(sqrt.(mpts.s.v[1, :].^2 .+ mpts.s.v[2, :].^2)) 
get_Δu(mpts)   = @views vec(sqrt.(mpts.s.u[1, :].^2 .+ mpts.s.u[2, :].^2)) 
get_coh0(mpts) = @views vec(mpts.s.c₀./1e3)
get_phi0(mpts) = @views vec(mpts.s.ϕ.*180.0/π)

# Variable extraction configuration dictionary with default colorbar limits
const MPTS_VAR = Dict(
    "epII" => (
        data=get_epII, 
        label=L"\epsilon_{\mathrm{II}}^{\mathrm{acc}}"*" [-]",  
        name="Plastic strain",
        cb=:viridis,
        cblim=(0.0, 1.5)
    ),
    "P" => (
        data=get_P,    
        label=L"p"*" [kPa]",                                     
        name="Pressure",
        cb=:viridis,
        cblim=(0.0, 150)
    ),
    "J" => (
        data=get_J,    
        label=L"J",                                     
        name="Deformation determinant"*" [-]",
        cb=:viridis,
        cblim=(0.5, 1.5)
    ),
    "Δu" => (
        data=get_Δu,
        label=L"\Delta u"*" [m]",
        name="Displacement magnitude",
        cb=:viridis,
        cblim=(0.0, 5.0)
    ),
    "v" => (
        data=get_v,
        label=L"v"*" [m/s]",
        name="Velocity magnitude",
        cb=:viridis,
        cblim=(0.0, 5.0)
    ),    
    "coh0" => (
        data=get_coh0,
        label=L"c_0(x_p)"*" [kPa]",
        name="Initial cohesion",
        cb = :vik,
        cblim=(10.0, 30.0)
    ),
    "phi0" => (
        data=get_phi0,
        label=L"$\phi_0(x_p)$"*" [deg]",
        name="Initial friction angle",
        cb = :viridis,
        cblim=(20.0, 40.0)
    )
)
const MESH_VAR = Dict(
    "m" => (
        data=get_epII, 
        label=L"m_{s}"*" [kg]",  
        name="Nodal solid mass",
        cb=:viridis,
        cblim=(100.0, 1000.0)
    ),
)

"""
    get_variable_plot_options() -> Vector

Generate kwargser-compatible plot options from MPTS_VAR.

# Returns
- Vector of named tuples with all fields from MPTS_VAR plus the variable key

# Example
```julia
opts = get_variable_plot_options()
# Returns [(;mpts=(key="epII", data=..., label=..., cb=..., cblim=...)), ...]
```
"""
function get_variable_plot_options()
    mpts = [(;mpts=merge((;key=k), v)) for (k, v) ∈ MPTS_VAR]
    mesh = [(;mesh=merge((;key=k), v)) for (k, v) ∈ MESH_VAR]
    return vcat(mpts,mesh)
end