# 3D figures take their height from the y width, not the z height

**Status: open; fix proposed, not applied.**

`bake` (`src/home/plot/bake.jl`) and `plot_initial_fields` (`src/home/plot/initial.jl`) size figures as
`dims = dpi .* (L ./ L[1])` and plot with `size=(dims[1], dims[2])`. In 3D, `L = [Lx, Ly, Lz]` and the
plot shows x–z, so the height should come from `L[end]`, not `L[2]`. The slump hides it (`Ly = Lz`); a
drawn 20 × 10 m domain extruded to a 4 m width gives a 500×100 figure instead of 500×250. Fix:
`dims = dpi .* ([L[1], L[end]] ./ L[1])` in both places.

Separately reported, not reproduced: in the live GR window (not the saved PNGs) the aspect ratio looked
inverted for a drawn problem.
