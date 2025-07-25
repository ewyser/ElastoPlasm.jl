function getinfo(L,nel)
    ndim = length(L)
    nn = 4^ndim
    if ndim == 1
        L   = L
        h   = L/nel
    elseif ndim == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif ndim == 3
        L   = [L[1],L[2],ceil(L[3])]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    else 
        throw(error("UnsupportedMeshGeometry: dim > 3"))
    end
    return L,h,ndim,nn
end
function get_coords(ndim,nn,L,h;ghosts::Vector=[0.0])
    if ndim == 1
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        x   = collect(first(x0):h[1]:last(x0))
        nno = [length(x),length(x)] 
        nel = [nno[1]-1 ,nno[1]-1 ]
        x   = vcat(vec(x))
    elseif ndim == 2
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        z0  = [0.0-ghosts[2],L[2]+ghosts[2]]
        x,z = collect(first(x0):h[1]:last(x0)),collect(first(z0):h[2]:last(z0))
        nno = [length(x),length(z),length(x)*length(z)  ] 
        nel = [nno[1]-1 ,nno[2]-1 ,(nno[1]-1)*(nno[2]-1)]

        x   = repeat((reshape(x,1     ,nno[1])),nno[2],1     )
        z   = repeat((reshape(z,nno[2],1     )),     1,nno[1])
        x   = vcat(vec(x)',vec(z)')
    elseif ndim == 3
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        y0  = [0.0-ghosts[2],L[2]+ghosts[2]]
        z0  = [0.0-ghosts[3],L[3]+2.0*h[3]+ghosts[3]]
        x   = collect(first(x0):h[1]:last(x0)) 
        y   = collect(first(y0):h[2]:last(y0)) 
        z   = collect(first(z0):h[3]:last(z0))   
        nno = [length(x),length(y),length(z),length(x)*length(y)*length(z)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]

        x   = repeat((reshape(x,1     ,nno[1],1     )),nno[3],1     ,nno[2])
        y   = repeat((reshape(y,1     ,1     ,nno[2])),nno[3],nno[1],1     )
        z   = repeat((reshape(z,nno[3],1     ,1     )),1     ,nno[1],nno[2]) 
        x   = vcat(vec(x)',vec(y)',vec(z)')
    end
    return x,nel,nno
end