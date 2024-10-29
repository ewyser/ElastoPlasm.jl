@views @kernel inbounds = true function euler(meD,Δt,η)
    no = @index(Global)
    for dim ∈ 1:meD.nD
        if no≤meD.nno[end] 
            if iszero(meD.mn[no])
                m = 0.0                                                        #(1,)                       
            else
                m = (1.0/meD.mn[no])*meD.bc[no,dim]                            #(1,)
            end
            meD.Dn[no,dim] = η*norm(meD.oobf[no,:])*sign(meD.pn[no,dim]*m)     #(2,)
            meD.fn[no,dim] = meD.oobf[no,dim]-meD.Dn[no,dim]                   #(2,)
            meD.an[no,dim] = meD.fn[no,dim]*m                                  #(2,)
            meD.vn[no,dim] = (meD.pn[no,dim]+Δt*meD.fn[no,dim])*m              #(2,)   
        end
    end
end