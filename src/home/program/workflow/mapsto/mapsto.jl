function init_mapsto(dim::Number,instr::Dict) 
    if instr[:fwrk][:deform] == "finite"
        kernel0 = transform(CPU())
    else
        kernel0 = nothing
    end
    kernel2 = euler(CPU())
    if instr[:fwrk][:trsfr] == "musl"
        if dim == 1
            kernel1 = flip_1d_p2n(CPU())
        elseif dim == 2
            kernel1 = flip_2d_p2n(CPU())
        elseif dim == 3
            kernel1 = flip_3d_p2n(CPU())
        end
        kernel3 = flip_nd_n2p(CPU())
        kernel3a = augm_momentum(CPU())
        kernel3b = augm_velocity(CPU())
        kernel3c = augm_displacement(CPU())
        return Dict(:map  => (;σᵢ! = kernel0, p2n! = kernel1 , solve! = kernel2 , n2p! = kernel3, ), 
                    :augm => (;p2n! = kernel3a, solve! = kernel3b, Δu!  = kernel3c,),)
    elseif instr[:fwrk][:trsfr] == "tpic"
        if dim == 2
            kernel1 = tpic_2d_p2n(CPU())
            kernel3 = pic_nd_n2p(CPU())
        elseif dim == 3
            kernel1 = tpic_3d_p2n(CPU())
            kernel3 = pic_nd_n2p(CPU())
        end
        return Dict(:map  => (;σᵢ! = kernel0, p2n! = kernel1, solve! = kernel2, n2p! = kernel3, ),)
    else
        return throw(ArgumentError("$(instr[:fwrk][:trsfr]) is an unsupported transfer scheme"))
    end    
end
function mapsto(mp::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2},dt::T2,instr::Dict) where {T1,T2}
    # maps material point to node
    p2n(mp,mesh,g,instr)
    # solve Eulerian momentum equation
    solve(mesh,dt,instr)
    # maps back solution to material point
    n2p(mp,mesh,dt,instr)
    if instr[:fwrk][:trsfr] == "musl"
        augm(mp,mesh,dt,instr)
    end
    return nothing
end



















#=
@views function mapstoN!(mp,mesh,g)
    # initialize nodal quantities
    mesh.mᵢ  .= 0.0
    mesh.p  .= 0.0
    mesh.oobf.= 0.0
    # mapping back to mesh
    for dim ∈ 1:mesh.dim
        lk = ReentrantLock()
        @threads for p ∈ 1:mp.nmp
            # accumulation
            lock(lk) do 
                if dim == 1 
                    mesh.mᵢ[mp.p2n[:,p]].+= mp.ϕ∂ϕ[:,p,1].*mp.m[p] 
                end
                mesh.p[  mp.p2n[:,p],dim].+= mp.ϕ∂ϕ[:,p,1].*(mp.m[p]*mp.v[p,dim])
                mesh.oobf[mp.p2n[:,p],dim].+= mp.ϕ∂ϕ[:,p,1].*(mp.m[p]*g[dim]      )
                mesh.oobf[mp.p2n[:,p],dim].-= mp.V[p].*(mp.B[dim:mesh.dim:end,:,p]*mp.σ[:,p]) 
            end
        end
    end
    return nothing
end
=#