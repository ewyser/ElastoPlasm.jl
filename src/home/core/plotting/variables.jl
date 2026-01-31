export VARIABLE_EXTRACTORS, get_variable_plot_options

# Variable extraction functions for metanalysis
extract_epII(mpts) = vec(mpts.s.ϵpII[1, :])
extract_P(mpts)    = -vec(mean(mpts.s.σᵢ[1:2, :], dims=1)) / 1e3
extract_J(mpts)    = vec(mpts.J)

# Variable extraction configuration dictionary with default colorbar limits
const VARIABLE_EXTRACTORS = Dict(
    "epII" => (extract=extract_epII, label=L"\epsilon_{\mathrm{II}}^{\mathrm{acc}}", name="plastic strain",            cblim=(0.0, 1.5)),
    "P"    => (extract=extract_P,    label=L"p",                                     name="pressure",                  cblim=nothing),
    "J"    => (extract=extract_J,    label=L"J",                                     name="deformation determinant",  cblim=(0.5, 1.5)),
)

"""
    get_variable_plot_options() -> Vector

Generate kwargser-compatible plot options from VARIABLE_EXTRACTORS.

# Returns
- Vector of named tuples in the format `(;mpts=(name=..., cblim=...))`

# Example
```julia
opts = get_variable_plot_options()
# Returns [(;mpts=(name="epII", cblim=(0.0,1.5))), ...]
```
"""
function get_variable_plot_options()
    return [(;mpts=(name=k, cblim=v.cblim)) for (k, v) in VARIABLE_EXTRACTORS]
end
