function getinfo(L,nel)
    nD = length(L)
    nn = 4^nD
    if nD == 1
        L   = L
        h   = L/nel
    elseif nD == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif nD == 3
        L   = [L[1],L[2],ceil(L[3])]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    else 
        err_msg = "dim(L = ($(L)))>3: unsupported mesh geometry"
        throw(error(err_msg))
    end
    return L,h,nD,nn
end
function getcoords(nD,nn,L,h;ghosts::Vector=[0.0])
    if nD == 1
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        xn  = collect(first(x0):h[1]:last(x0))
        nno = [length(xn),length(xn)] 
        nel = [nno[1]-1,nno[1]-1    ]
        #xt  = vcat([1],[2],3*ones(Int32,nno[1]-4),[4],[1])
        return xn,nel,nno
    elseif nD == 2
        x0,z0 = [0.0-ghosts[1],L[1]+ghosts[1]],[0.0-ghosts[2],L[2]+2.0*h[2]+ghosts[2]]
        xn,zn = collect(first(x0):h[1]:last(x0)),collect(first(z0):h[2]:last(z0))
        nno   = [length(xn),length(zn),length(xn)*length(zn)] 
        nel   = [nno[1]-1  ,nno[2]-1  ,(nno[1]-1)*(nno[2]-1)]
        #xt    = vcat([1],[2],3*ones(Int32,nno[1]-4),[4],[1])
        zt    = vcat([1],[2],3*ones(Int32,nno[2]-4),[4],[1])
        xn    = hcat(vec(repeat(xn',nno[2],1)),vec(repeat(zn,1,nno[1])))
        #xt    = hcat(vec(repeat(xt',nno[2],1)),vec(repeat(zt,1,nno[1])))
        return xn,nel,nno
    elseif nD == 3
        x0 = [0.0-ghosts[1],L[1]+ghosts[1]]
        y0 = [0.0-ghosts[2],L[2]+ghosts[2]]
        z0 = [0.0-ghosts[3],L[3]+2.0*h[3]+ghosts[3]]
        xn  = collect(first(x0): h[1]:last(x0) ) 
        yn  = collect(first(y0): h[2]:last(y0) ) 
        zn  = collect(last(z0) :-h[3]:first(z0))
        #xt  = repeat([3],length(xn))
        #yt  = repeat([3],length(yn))
        #zt  = repeat([3],length(zn))
        #xt[1]     = yt[1]     = zt[1]     = 1
        #xt[2]     = yt[2]     = zt[2]     = 2
        #xt[end-1] = yt[end-1] = zt[end-1] = 4
        #xt[end  ] = yt[end  ] = zt[end  ] = 1      
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        x   = hcat(vec(xn),vec(yn),vec(zn))
        #xt  = (xt'.*ones(Int64,nno[3],1     ))     .*ones(Int64,1,1,nno[2])
        #zt  = (     ones(Int64,nno[1],1     )'.*zt).*ones(Int64,1,1,nno[2])
        #yt  = (     ones(Int64,nno[3],nno[1]))     .*reshape(yt,1,1,nno[2])
        #xt  = hcat(vec(xt),vec(yt),vec(zt))
        return xn,nel,nno
    end
end