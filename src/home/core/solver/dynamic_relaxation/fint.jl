#=
function ϵᵢⱼ(bᵢⱼ)
    λ,n = eigen(bᵢⱼ,sortby=nothing)
    return T2(0.5).*(n*diagm(log.(λ))*n')
end
ϵᵢⱼ = ϵᵢⱼ(bᵢⱼ)
=#


@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,1},mesh::Mesh{T1,T2,1},basis::Basis{T1,T2,1},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms, Ω = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx   = get_vector(mpts.s.σᵢ[p])[1]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx) - N * (ms * g[1])
        end
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,2,CM,TM,TV,TS,ST},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},g::Vector{T2},dt::T2,Del) where {T1,T2,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        # compute velocity & displacement gradients
        ∇v = zeros(MMatrix{2,2,T2})
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, _, ∂N = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            for i ∈ 1:(length(mesh.prprt.nel)-1)
                for j ∈ 1:(length(mesh.prprt.nel)-1)
                    ∇v[i,j] += ∂N[j]*mesh.s.v[no][i]
                end
            end
        end
        ∇vᵢⱼ = SMatrix(∇v)
        
        # compute incremental deformation and update
        ΔFᵢⱼ = SMatrix{2,2,T2}(I) + dt .* ∇vᵢⱼ
        ωᵢⱼ  = T2(0.5) .* (∇vᵢⱼ - ∇vᵢⱼ')

        # calculate elastic strains & spins
        ϵᵢⱼ = T2(0.5) .* (ΔFᵢⱼ + ΔFᵢⱼ') .- SMatrix{2,2,T2}(I)
        ωᵢⱼ = T2(0.5) .* (∇vᵢⱼ - ∇vᵢⱼ')
        
        # update cauchy stress tensor
        σᵢ   = get_vector(mpts.s.σn[p])
        σJᵢⱼ = mutate(σᵢ, T2(1.0), :tensor)
        σJᵢⱼ = σJᵢⱼ*ωᵢⱼ'+σJᵢⱼ'*ωᵢⱼ
        σᵢ   = σᵢ + typeof(σᵢ)(Del * mutate(ϵᵢⱼ, T2(2.0), :voigt) .+ mutate(σJᵢⱼ, T2(1.0), :voigt))
        
        # buffered values
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σxx, σyy, σxy = σᵢ[1], σᵢ[2], σᵢ[3]
        
        # calculate oobf
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no] -= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end

        # save deformation gradient and cauchy stress
        mpts.s.Fᵢⱼ[p] = ΔFᵢⱼ * mpts.s.Fn[p]
        mpts.s.σᵢ[p]    = CauchyStress(σᵢ)
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,2,CM,TM,TV,TS,ST},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},g::Vector{T2},dt::T2,Kc::T2,Gc::T2) where {T1,T2,CM,TM,TV,TS,ST<:LogarithmicStrain}
    mp = @index(Global)
    if mp ≤ mpts.nmp
        # Compute velocity & displacement gradients
        ∇vᵢⱼ = zeros(MMatrix{2,2,T2})
        for nn ∈ 1:mesh.prprt.nn
            # Buffering & compute basis functions on-the-fly
            no, _, ∂N = eval_basis(mpts, mesh, basis, mp, nn)
            if iszero(no) continue end
            ∇vᵢⱼ += ∂N * mesh.s.v[no]'
        end
        # Compute incremental deformation, update current deformation gradient and get the deformation determinant
        ΔFᵢⱼ = SMatrix{2,2,T2}(I) + dt .* SMatrix(∇vᵢⱼ)
        Fᵢⱼ  = ΔFᵢⱼ*mpts.s.Fn[mp]
        J⁻¹  = T2(1.0)/det(Fᵢⱼ)

        # Compute trial left Cauchy-Green tensor
        λ,n = eigen(mpts.s.ϵn[mp])
        bᵢⱼ = ΔFᵢⱼ*(n*diagm(exp.(T2(2.0)*λ))*n')*ΔFᵢⱼ'
        
        # Compute trial logarithmic strain tensor
        λ,n = eigen(bᵢⱼ)
        ϵᵢⱼ = T2(0.5)*(n*diagm(log.(λ))*n')
        
        # Compute volumetric strain
        ϵV  = tr(ϵᵢⱼ) / T2(3.0)

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
            no, N, ∂N = eval_basis(mpts, mesh, basis, mp, nn)
            # Transform shape function derivatives to current configuration
            ∂N = ∂N' *inv(ΔFᵢⱼ)
            # Assemble oobf contributions from material point level to mesh level
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * (∂N[1] * σxx + ∂N[2] * σxy)
            @atom mesh.s.oobf[2,no] -= Ω * (∂N[1] * σxy + ∂N[2] * σyy) - N * (ms * g[2])
        end

        # Store deformation gradient, logarithmic strain and Kirchhoff stress
        mpts.s.Fᵢⱼ[mp] = Fᵢⱼ
        mpts.s.ϵᵢⱼ[mp] = LogarithmicStrain(SMatrix{2,2,T2,4}(ϵᵢⱼ))
        mpts.s.τᵢ[mp]  = KirchhoffStress(SVector{3,T2}(τxx, τyy, τxy))
    end
end

@kernel inbounds = true function oobf_assembly(mpts::Point{T1,T2,3},mesh::Mesh{T1,T2,3},basis::Basis{T1,T2,3},g::Vector{T2}) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp
        ms , Ω        = mpts.s.ρ[p]*mpts.Ω[p], mpts.Ω[p]
        σ             = get_vector(mpts.s.σᵢ[p])
        σxx, σyy, σzz = σ[1], σ[2], σ[3]
        # (was `mpts.s.σᵢ[6,p]`/`[4,p]`/`[5,p]` — matrix-style indexing into what has been
        #  a Vector of per-particle statics for a while now, i.e. dead on arrival; the
        #  typed-tensor port forces it to be written correctly)
        σyx, σzy, σzx = σ[6], σ[4], σ[5]
        for nn ∈ 1:mesh.prprt.nn
            no, N, ∂N = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            @atom mesh.s.oobf[1,no] -= Ω * ( ∂N[1] * σxx + ∂N[2] * σyx + ∂N[3] * σzx)
            @atom mesh.s.oobf[2,no] -= Ω * ( ∂N[1] * σyx + ∂N[2] * σyy + ∂N[3] * σzy)
            @atom mesh.s.oobf[3,no] -= Ω * ( ∂N[1] * σzx + ∂N[2] * σzy + ∂N[3] * σzz) - N * (ms * g[3])
        end
    end
end
