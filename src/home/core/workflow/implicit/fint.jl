#=
function ϵᵢⱼ(bᵢⱼ)
    λ,n = eigen(bᵢⱼ,sortby=nothing)
    return T2(0.5).*(n*diagm(log.(λ))*n')
end
ϵᵢⱼ = ϵᵢⱼ(bᵢⱼ)
=#


@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:OneDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms, Ω = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx   = mpts.s.σᵢ[1,p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,D,B,E,R},mesh::Mesh{T1,T2,D},g::Vector{T2},dt::T2,Del) where {T1,T2,D<:TwoDimension,B<:AbstractBasis,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        # compute velocity & displacement gradients
        ∇vᵢⱼ = zeros(T2,2,2)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end                                     
            for i ∈ 1:mesh.prprt.dim  
                for j ∈ 1:mesh.prprt.dim
                    ∇vᵢⱼ[i,j]+= ∂N[j]*mesh.s.v[i,no]
                end
            end
        end
        
        # compute incremental deformation and update
        ΔFᵢⱼ = mpts.δᵢⱼ+(dt.*∇vᵢⱼ)
        ωᵢⱼ  = T2(0.5).*(∇vᵢⱼ-∇vᵢⱼ')

        # calculate elastic strains & spins
        ϵᵢⱼ = T2(0.5).*(ΔFᵢⱼ+ΔFᵢⱼ').-mpts.δᵢⱼ[:,:] 
        ωᵢⱼ = T2(0.5).*(∇vᵢⱼ-∇vᵢⱼ')
        
        # update cauchy stress tensor
        σᵢ   = copy(mpts.s.σn[:,p])
        σJᵢⱼ = mutate(σᵢ,T2(1.0),:tensor)
        σJᵢⱼ = σJᵢⱼ*ωᵢⱼ'+σJᵢⱼ'*ωᵢⱼ
        σᵢ  += Del*mutate(ϵᵢⱼ,T2(2.0),:voigt).+mutate(σJᵢⱼ,T2(1.0),:voigt)
        
        # buffered values
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx, σyy, σxy = σᵢ[1], σᵢ[2], σᵢ[3]
        
        # calculate oobf
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no] -= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end

        # save deformation gradient and cauchy stress
        mpts.s.Fᵢⱼ[:,:,p] .= ΔFᵢⱼ*mpts.s.Fn[:,:,p]
        mpts.s.σᵢ[:,p]    .= σᵢ
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,D,B,E,R},mesh::Mesh{T1,T2,D},g::Vector{T2},dt::T2,Kc::T2,Gc::T2) where {T1,T2,D<:TwoDimension,B<:AbstractBasis,E<:FiniteElasticity,R<:AbstractRheology}
    mp = @index(Global)
    if mp ≤ mpts.nmp
        # Get previous deformation gradient and logarithmic strain tensor
        Fn = SMatrix{2,2,T2}(mpts.s.Fn[:, :, mp])
        ϵn = SMatrix{2,2,T2}(mpts.s.ϵn[:, :, mp])

        # Compute velocity & displacement gradients
        ∇vᵢⱼ = zeros(MMatrix{2,2,T2})
        for nn ∈ 1:mesh.prprt.nn
            # Buffering & compute basis functions on-the-fly
            no, _, ∂N = basis(mpts, mesh, mp, nn)
            if iszero(no) continue end                                     
            for i ∈ 1:mesh.prprt.dim  
                for j ∈ 1:mesh.prprt.dim
                    ∇vᵢⱼ[i,j]+= ∂N[j]*mesh.s.v[i,no]
                end
            end
        end
        
        # Compute incremental deformation, update current deformation gradient and get the deformation determinant
        ΔFᵢⱼ = mpts.δᵢⱼ+(dt.*∇vᵢⱼ)
        Fᵢⱼ  = ΔFᵢⱼ*Fn
        J⁻¹  = T2(1.0)/det(Fᵢⱼ)

        # Compute trial left Cauchy-Green tensor
        λ,n = eigen(ϵn)
        bᵢⱼ = ΔFᵢⱼ*(n*diagm(exp.(T2(2.0).*λ))*n')*ΔFᵢⱼ'
        
        # Compute trial logarithmic strain tensor
        λ,n = eigen(bᵢⱼ)
        ϵᵢⱼ = T2(0.5).*(n*diagm(log.(λ))*n')
        
        # Compute volumetric strain
        ϵV  = (ϵᵢⱼ[1,1] + ϵᵢⱼ[2,2]) / T2(3.0)

        # Calculate pressure (positive in compression)
        P   = - T2(3.0) * Kc * ϵV

        # Calculate Krichhoff stress tensor
        τxx = T2(2.0) * Gc * (ϵᵢⱼ[1,1] - ϵV) - P
        τyy = T2(2.0) * Gc * (ϵᵢⱼ[2,2] - ϵV) - P
        τxy = T2(2.0) * Gc *  ϵᵢⱼ[1,2]

        # Buffered values
        ms , Ω        = mpts.s.ρ₀[mp]*mpts.Ω₀[mp], mpts.Ω[mp]*det(ΔFᵢⱼ)
        σxx, σyy, σxy = J⁻¹*τxx, J⁻¹*τyy, J⁻¹*τxy
        
        # Accumulate oobf
        for nn ∈ 1:mesh.prprt.nn
            # Buffering & compute basis functions on-the-fly
            no, N, ∂N = basis(mpts, mesh, mp, nn)
            # Transform shape function derivatives to current configuration
            ∂N = ∂N' *inv(ΔFᵢⱼ)
            # Assemble oobf contributions from material point level to mesh level
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no] -= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end

        # Store deformation gradient, logarithmic strain and Kirchhoff stress
        mpts.s.Fᵢⱼ[:,:,mp] .= Fᵢⱼ
        mpts.s.ϵᵢⱼ[:,:,mp] .= ϵᵢⱼ
        mpts.s.τᵢ[:,mp]    .= [τxx, τyy, τxy]
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D<:ThreeDimension}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx, σyy, σzz = mpts.s.σᵢ[1,p], mpts.s.σᵢ[2,p], mpts.s.σᵢ[3,p]
        σyx, σzy, σzx = mpts.s.σᵢ[6,p], mpts.s.σᵢ[4,p], mpts.s.σᵢ[5,p]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * ( ∂N[1] * σxx + ∂N[2] * σyx + ∂N[3] * σzx)
            @atom mesh.s.oobf[2,no] -= Ω * ( ∂N[1] * σyx + ∂N[2] * σyy + ∂N[3] * σzy)
            @atom mesh.s.oobf[3,no] -= Ω * ( ∂N[1] * σzx + ∂N[2] * σzy + ∂N[3] * σzz) - N * (ms * g[3])
        end
    end
end
