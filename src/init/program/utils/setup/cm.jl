function D(E,ν,nD)
    Kc,Gc = E/(3.0*(1.0-2.0*ν)),E/(2.0*(1.0+ν))                                # bulk & shear modulus               [Pa]
    if nD == 1
        D  = [ 
            Kc+4/3*Gc 0.0;
            0.0       Gc ;
            0.0       0.0]
    elseif nD == 2
        D  = [ 
            Kc+4/3*Gc Kc-2/3*Gc 0.0 ;
            Kc-2/3*Gc Kc+4/3*Gc 0.0 ;
            0.0       0.0       Gc  ]
    elseif nD == 3
        D  = [ 
            Kc+4/3*Gc Kc-2/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
            Kc-2/3*Gc Kc+4/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
            Kc-2/3*Gc Kc-2/3*Gc Kc+4/3*Gc 0.0 0.0 0.0;
            0.0       0.0       0.0       Gc  0.0 0.0;
            0.0       0.0       0.0       0.0 Gc  0.0;
            0.0       0.0       0.0       0.0 0.0 Gc ;]
    end
    return Kc,Gc,D
end
function cm(dim,instr; E::Real=1.0e6,ν::Real=0.3,ρ0::Real= 2700.0)
    # independant physical constant          
    K,G,Del = D(E,ν,dim)                                                  # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    c       = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    Hp      = -30.0e3                                                           # softening modulus
    # constitutive model param.
    cmParam = (;
        cmType   = instr[:plast][:constitutive], 
        nonlocal = instr[:nonloc],
        E   = E, 
        ν   = ν, 
        Kc  = K, 
        Gc  = G, 
        Del = Del, 
        Hp  = Hp, 
        c0  = c0,
        cr  = cr,
        ϕ0  = ϕ0,
        ϕr  = ϕr,
        ρ0  = ρ0,
        c   = c,
    )
    return cmParam
end