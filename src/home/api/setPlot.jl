"""
    configPlot(titf=12, gf=12, tickf=10, lf=10, lw=2, fs=:box, l=nothing, g=false)

Configures the default styling parameters for plots, allowing customization of various plot elements such as fonts, line widths, and grid visibility.

# Arguments
- `titf::Int`: Font size for the plot title. Controls the size of the text used for the plot title. (default: `12`).
- `gf::Int`: Font size for guide elements, such as axis labels and legends. Determines the size of the text used for these elements. (default: `12`).
- `tickf::Int`: Font size for tick labels. Specifies the size of the text used for the axis tick labels. (default: `10`).
- `lf::Int`: Font size for the legend text. Adjusts the size of the text displayed in the legend. (default: `10`).
- `lw::Int`: Line width for plot lines. Defines the thickness of lines used in the plot. (default: `2`).
- `fs::Symbol`: Frame style for the plot. This controls the appearance of the plot's frame, such as `:box` for a boxed frame or `:none` for no frame. (default: `:box`).
- `l::Union{Nothing, String}`: Label for the plot. Allows setting a custom label for the plot. If set to `nothing`, no label is displayed. (default: `nothing`).
- `g::Bool`: Grid display option. Enables or disables grid lines on the plot. If set to `true`, grid lines will be shown; if `false`, they will be hidden. (default: `false`).

# Description
The `configPlot` function sets default parameters that dictate the visual style of plots generated in your Julia session. By adjusting these parameters, users can tailor the look and feel of their plots to meet specific presentation or publication standards. It is recommended to call this function before generating any plots to ensure that the desired settings are applied consistently across all visualizations.

This function is particularly useful in customizing the plot's aesthetics in a simple and centralized manner, eliminating the need to repeatedly specify styling options for individual plots.

# Example
```julia
julia> configPlot(titf=14, gf=12, tickf=8, lf=12, lw=3, fs=:box, l="My Plot", g=true)
```
"""
function configPlot(titf=12,gf=12,tickf=10,lf=10,lw=2,fs=:box,l=nothing,g=false)
    default(
        fontfamily  = "Computer Modern",
        titlefont   = titf, 
        guidefont   = gf,  
        tickfont    = tickf, 
        legendfont  = lf,
        linewidth   = lw,
        framestyle  = fs,
        label       = l,
        grid        = g,
    )
    return nothing
end
export configPlot


