"""
    get_geom(nel::Vector{T1}, L::Vector{T2}, instr) where {T1,T2}

Compute mesh geometry parameters based on the number of elements, domain size, and basis type.

# Arguments
- `nel::Vector{T1}`: Number of elements in each spatial direction.
- `L::Vector{T2}`: Length of the domain in each spatial direction.
- `instr`: Instruction dictionary containing basis information.

# Returns
- `ndim::T1`: Number of spatial dimensions.
- `nn::T1`: Number of nodes.
- `h::Vector{T2}`: Element size in each direction.

# Example
```julia
ndim, nn, h = get_geom([10, 10], [1.0, 1.0], Dict(:basis => Dict(:which => "bsmpm")))
```
"""
function get_geom(nel::Vector,L::Vector,instr)
    # unpack arithmetic precision
    T1,T2 = first(instr[:dtype].T0),last(instr[:dtype].T0) 
    if instr[:basis][:which] == "bsmpm"
        ndim,nn,h = length(L),4^length(L),min.(L./nel,L./4)
    else
        ndim,nn,h = length(L),4^length(L),(L./nel)
    end
    return (; ndim = T1(ndim), nn = T1(nn), h =T2.(h), nel = T1.(nel), L = T2.(L))
end