"""
    get_elastic_stiffness(E, ν, ndim) -> Kc, Gc, D

Compute the bulk modulus, shear modulus, and elastic stiffness matrix for a given Young's modulus, Poisson's ratio, and spatial dimension.

# Arguments
- `E`: Young's modulus (Pa).
- `ν`: Poisson's ratio.
- `ndim`: Number of spatial dimensions (1, 2, or 3).

# Returns
- `Kc`: Bulk modulus.
- `Gc`: Shear modulus.
- `Del`: Elastic stiffness matrix (Voigt notation).

# Example
```julia
Kc, Gc, D = get_elastic_stiffness(1.0e6, 0.3, 2)
```
"""
function get_elastic_stiffness(E::T2,ν::T2,ndim::T1) where {T1,T2}
    Kc,Gc = E/(3.0*(1.0-2.0*ν)),E/(2.0*(1.0+ν))  
    i ,I  = 1.0-ν              ,(1.0-2.0*ν)/2.0                              # bulk & shear modulus               [Pa] 
    if ndim == 1
        D  = [ 
            i   0.0;
            0.0 I  ;
        ]
    elseif ndim == 2
        D  = [ 
            i   ν   0.0 ;
            ν   i   0.0 ;
            0.0 0.0 I   ;
        ]
        
    elseif ndim == 3
        D  = [ 
            i   ν   ν   0.0 0.0 0.0;
            ν   i   ν   0.0 0.0 0.0;
            ν   ν   i   0.0 0.0 0.0;
            0.0 0.0 0.0 I   0.0 0.0;
            0.0 0.0 0.0 0.0 I   0.0;
            0.0 0.0 0.0 0.0 0.0 I  ;
        ]
    end
    Del = E/((1.0+ν)*(1.0-2.0*ν)).*D
    return T2(Kc),T2(Gc),T2.(Del)
end

"""
    Del(Kc, Gc, n, ::Val{dim}) -> Matrix{T2}

Compute the elastic stiffness matrix for a given bulk modulus, shear modulus, porosity, and dimension.

# Arguments
- `Kc`: Bulk modulus.
- `Gc`: Shear modulus.
- `n`: Porosity.
- `dim`: Dimension as Val{:1}, Val{:2}, or Val{:3}.

# Returns
- `Matrix{T2}`: Elastic stiffness matrix.

# Example
```julia
Del(Kc, Gc, n, Val(:2))
```
"""
function Del(Kc::T2,Gc::T2,n::T2,::Val{:1}) where {T2}
    n      = max(n,1e-2)                                # small fix for small porosity value             
    Kc′,Gc = (Kc/n),Gc                                  # bulk & shear modulus [Pa]
    Del    = [                                          #                
        Kc′+4/3*Gc 0.0;                                 #
        0.0        Gc ;                                 #
        0.0        0.0]                                 # elastic stiffness matrix [Pa]
    return T2.(Del)
end
function Del(Kc::T2,Gc::T2,n::T2,::Val{:2}) where {T2}
    n      = max(n,1e-2)                                # small fix for small porosity value             
    Kc′,Gc = (Kc/n),Gc                                  # bulk & shear modulus [Pa]
    Del    = [                                          #
        Kc′+4/3*Gc Kc′-2/3*Gc 0.0 ;                     #
        Kc′-2/3*Gc Kc′+4/3*Gc 0.0 ;                     #
        0.0        0.0        Gc  ]                     # elastic stiffness matrix [Pa]
    return T2.(Del)
end
function Del(Kc::T2,Gc::T2,n::T2,::Val{:3}) where {T2}
    n      = max(n,1e-2)                                # small fix for small porosity value             
    Kc′,Gc = (Kc/n),Gc                                  # bulk & shear modulus [Pa]
    Del    = [                                          #
        Kc′+4/3*Gc Kc′-2/3*Gc Kc′-2/3*Gc 0.0 0.0 0.0;   #
        Kc′-2/3*Gc Kc′+4/3*Gc Kc′-2/3*Gc 0.0 0.0 0.0;   #
        Kc′-2/3*Gc Kc′-2/3*Gc Kc′+4/3*Gc 0.0 0.0 0.0;   #
        0.0        0.0        0.0        Gc  0.0 0.0;   #
        0.0        0.0        0.0        0.0 Gc  0.0;   #
        0.0        0.0        0.0        0.0 0.0 Gc ;]  # elastic stiffness matrix [Pa]
    return T2.(Del)
end

"""
    setup_cmpr(mesh::Mesh{T1,T2}; E::T2=1.0e6, ν::T2=0.3, ρ0::T2=2700.0) -> NamedTuple

Set up the constitutive model parameters for a simulation, including elastic and plastic properties, based on mesh and instruction dictionary.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh object containing dimension information.
- `E::T2=1.0e6`: (Optional) Young's modulus (Pa).
- `ν::T2=0.3`: (Optional) Poisson's ratio.
- `ρ0::T2=2700.0`: (Optional) Initial density (kg/m³).

# Returns
- `NamedTuple`: Constitutive model parameters, including elastic, plastic, and nonlocal properties.

# Example
```julia
cmp = setup_cmpr(mesh; E=2.0e6, ν=0.25, ρ0=2500.0)
println(cmp.Kc)  # Bulk modulus
```
"""
function setup_cmpr(mesh::Mesh{T1,T2}; E::T2=T2(1.0e6),ν::T2=T2(0.3),ρ0::T2= T2(2700.0)) where {T1,T2}
    # independant physical constant          
    K,G,Del = get_elastic_stiffness(E,ν,mesh.prprt.dim)                         # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    c       = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    Hp      = -30.0e3                                                           # softening modulus
    # constitutive model param.
    cmp = (;
        # solid phase
        E   = T2(E), 
        ν   = T2(ν), 
        Kc  = T2(K), 
        Gc  = T2(G), 
        Del = T2.(Del), 
        Hp  = T2(Hp), 
        c0  = T2(c0),
        cr  = T2(cr),
        ϕ0  = T2(ϕ0),
        ϕr  = T2(ϕr),
        ρ0  = T2(ρ0),
        c   = T2(c),
        # thermal phase
        specific_heat_capacity = T2(10.0),  # J/(kg·K)
        thermal_conductivity   = T2(3000.0),     # W/(m·K)
        initial_temperature    = T2(293.15),  # K, 20°C in Kelvin
    )
    return cmp::NamedTuple
end


    #c   ::Vector{T2} # specific heat capacity vector
    #k   ::Vector{T2} # thermal conductivity vector
    #q   ::Matrix{T2} # heat flux array
    #T   ::Vector{T2} # temperature vector