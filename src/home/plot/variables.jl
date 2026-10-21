# Variable extraction configuration function with default colorbar limits
export get_mpts_variable_config, get_variable_plot_options

# Variable extraction functions for metanalysis. All return raw (SI) values;
# unit conversion for display is applied afterwards via each entry's `scale`.
get_epII(mpts) = getindex.(mpts.s.ϵpII, 1)
# mean pressure of the full 3×3 stress, positive in compression (in 2D, σzz included)
get_P(mpts)    = getfield.(mpts.s.σᵢ, :p)
get_J(mpts)    = mpts.J
get_n(mpts)    = mpts.n
get_u(mpts)    = [sqrt(sum(abs2, u)) for u in mpts.s.u]     # magnitude, any dimension
get_v(mpts)    = [sqrt(sum(abs2, v)) for v in mpts.s.v]
get_coh0(mpts) = [cmp.c₀ for cmp in mpts.s.cmp]
get_phi0(mpts) = [cmp.ϕ₀ for cmp in mpts.s.cmp]
get_domain(mpts) = mpts.ℓ
get_T(mpts)    = mpts.t.T

"""
    get_mpts_variable_config() -> Dict{String, NamedTuple}

Return a dictionary of material point variable plotting configurations, including extraction functions, labels, units, scaling factors, colorbar settings, and limits.

`data` always returns raw (SI) values; `scale` is the factor applied to those values
(e.g. `1/1e3` for Pa → kPa, `180/π` for rad → deg) to get the quantity in `unit` for display.
`cblim` is expressed in the scaled (display) unit.
"""
function get_mpts_variable_config()
    return Dict(
        "epII" => (
            data=get_epII,
            label=L"\epsilon_{\mathrm{II}}^{\mathrm{acc}}",
            unit = "[-]",
            scale = 1.0,
            name="plastic strain",
            cb=:viridis,
            cblim=(0.0, 1.5)
        ),
        "P" => (
            data=get_P,
            label=L"p",
            unit = "[kPa]",
            scale = 1/1e3,
            name="pressure",
            cb=:viridis,
            cblim=(0.0, 150)
        ),
        "J" => (
            data=get_J,
            label=L"J",
            unit = "[-]",
            scale = 1.0,
            name="deformation determinant",
            cb=:viridis,
            cblim=(0.5, 1.5)
        ),
        "n" => (
            data=get_n,
            label=L"n",
            unit = "[-]",
            scale = 1.0,
            name="porosity",
            cb=:vik,
            cblim=(-1.0, 1.0)
        ),
        "u" => (
            data=get_u,
            label=L"u",
            unit = "[m]",
            scale = 1.0,
            name="displacement magnitude",
            cb=:viridis,
            cblim=(0.0, 0.25)
        ),
        "v" => (
            data=get_v,
            label=L"v",
            unit = "[m/s]",
            scale = 1.0,
            name="velocity magnitude",
            cb=:viridis,
            cblim=(0.0, 5.0)
        ),
        "coh0" => (
            data=get_coh0,
            label=L"c_{0}",
            unit = "[kPa]",
            scale = 1/1e3,
            name="initial cohesion",
            cb = :vik,
            cblim=(10.0, 30.0)
        ),
        "phi0" => (
            data=get_phi0,
            label=L"\phi_{0}",
            unit = "[deg.]",
            scale = 180.0/π,
            name="initial friction angle",
            cb = :viridis,
            cblim=(20.0, 40.0)
        ),
        "domain" => (
            data  = get_domain,   # unused by the geometry render path (see what_plot_field)
            label = L"\ell",
            unit  = "[-]",
            scale = 1.0,
            name  = "domain",
            cb    = nothing,
            cblim = nothing,
        ),
        "T" => (
            data=get_T,
            label=L"T",
            unit = "[K]",
            scale = 1.0,
            name="temperature",
            cb=:thermal,
            cblim=(20.0, 293.15)
        ),
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
            scale = 1.0,
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