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
        err_msg = "dim(L = ($(L))) > 3: unsupported mesh geometry"
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
    elseif nD == 2
        x0    = [0.0-ghosts[1],L[1]+ghosts[1]]
        z0    = [0.0-ghosts[2],L[2]+2.0*h[2]+ghosts[2]]
        x,z   = collect(first(x0):h[1]:last(x0)),collect(first(z0):h[2]:last(z0))
        nno   = [length(x),length(z),length(x)*length(z)] 
        nel   = [nno[1]-1  ,nno[2]-1  ,(nno[1]-1)*(nno[2]-1)]

        x     = reshape(x,1        ,length(x))
        z     = reshape(z,length(z),1        )
        xn    =  repeat(x,length(z),1        )
        zn    =  repeat(z,        1,length(x))
    elseif nD == 3
        x0 = [0.0-ghosts[1],L[1]+ghosts[1]]
        y0 = [0.0-ghosts[2],L[2]+ghosts[2]]
        z0 = [0.0-ghosts[3],L[3]+2.0*h[3]+ghosts[3]]
        x   = collect(first(x0):h[1]:last(x0)) 
        y   = collect(first(y0):h[2]:last(y0)) 
        z   = collect(first(z0):h[3]:last(z0))   
        nno = [length(x),length(y),length(z),length(x)*length(y)*length(z)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]

        x  = reshape(x,1         ,length(x),1        )
        y  = reshape(y,1        ,1         ,length(y))
        z  = reshape(z,length(z),1         ,1        )

        xn = repeat(x,length(z),1        ,length(y))
        yn = repeat(y,length(z),length(x),1        )
        zn = repeat(z,        1,length(x),length(y)) 
    end

    if nD == 1
        xn = hcat(vec(xn))
    elseif nD == 2
        xn = hcat(vec(xn),vec(zn))
    elseif nD == 3
        xn = hcat(vec(xn),vec(yn),vec(zn))
    end

    return xn,nel,nno
end