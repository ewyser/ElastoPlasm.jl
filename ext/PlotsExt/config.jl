export config_plot

"""
    config_plot(; titf=12, gf=12, tickf=10, lf=10, lw=2, fs=:box, l=nothing, g=false)

Configure the default styling parameters for plots, allowing customization of fonts, line widths, frame style, labels, and grid visibility.

# Keyword Arguments
- `titf`: Font size for the plot title (default: `12`).
- `gf`: Font size for guide elements, such as axis labels and legends (default: `12`).
- `tickf`: Font size for tick labels (default: `10`).
- `lf`: Font size for legend text (default: `10`).
- `lw`: Line width for plot lines (default: `2`).
- `fs`: Frame style for the plot, e.g., `:box` or `:none` (default: `:box`).
- `l`: Label for the plot, or `nothing` for no label (default: `nothing`).
- `g`: Grid display option; `true` to show grid, `false` to hide (default: `false`).

# Returns
- `Nothing`. Configures plot defaults.

# Example
```julia
config_plot(titf=14, gf=12, tickf=8, lf=12, lw=3, fs=:box, l="My Plot", g=true)
```
"""
function config_plot(; titf=12, gf=12, tickf=10, lf=10, lw=2, fs=:box, l=nothing, g=false)
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
