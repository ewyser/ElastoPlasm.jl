"""
    mpts_populate(props; ni=2) -> Matrix

The candidate material points of a mesh (`MeshProperties`): `ni`×`ni` per element, at sub-cell
centres, as a D×N matrix, the vertical (last) axis fastest, then x, then y.
"""
function mpts_populate(props::MeshProperties{T1,T2,D}; ni::Integer=2) where {T1,T2,D}
    axs = [collect(props.xB[d,1]+(0.5*props.h[d]/ni):props.h[d]/ni:props.xB[d,2]-(0.5*props.h[d]/ni)) for d ∈ 1:D]
    pts = Iterators.product(axs[D], axs[1:D-1]...)
    return reshape(stack((t[2:end]..., t[1]) for t ∈ pts), D, :)
end
