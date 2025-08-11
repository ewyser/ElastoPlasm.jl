"""
    D(E, ν, ndim) -> Kc, Gc, D

Compute the bulk modulus, shear modulus, and elastic stiffness matrix for a given Young's modulus, Poisson's ratio, and spatial dimension.

# Arguments
- `E`: Young's modulus (Pa).
- `ν`: Poisson's ratio.
- `ndim`: Number of spatial dimensions (1, 2, or 3).

# Returns
- `Kc`: Bulk modulus.
- `Gc`: Shear modulus.
- `D`: Elastic stiffness matrix (Voigt notation).

# Example
```julia
Kc, Gc, D = D(1.0e6, 0.3, 2)
```
"""
function D(E,ν,ndim)
    Kc,Gc = E/(3.0*(1.0-2.0*ν)),E/(2.0*(1.0+ν))                                # bulk & shear modulus               [Pa]
    if ndim == 1
        D  = [ 
            Kc+4/3*Gc 0.0;
            0.0       Gc ;
            0.0       0.0]
    elseif ndim == 2
        D  = [ 
            Kc+4/3*Gc Kc-2/3*Gc 0.0 ;
            Kc-2/3*Gc Kc+4/3*Gc 0.0 ;
            0.0       0.0       Gc  ]
    elseif ndim == 3
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

"""
    setup_cmpr(mesh::Mesh{T1,T2}, instr::NamedTuple; E::T2=1.0e6, ν::T2=0.3, ρ0::T2=2700.0) -> NamedTuple

Set up the constitutive model parameters for a simulation, including elastic and plastic properties, based on mesh and instruction dictionary.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh object containing dimension information.
- `instr::NamedTuple`: Simulation instruction dictionary (must include plasticity and nonlocal settings).
- `E::T2=1.0e6`: (Optional) Young's modulus (Pa).
- `ν::T2=0.3`: (Optional) Poisson's ratio.
- `ρ0::T2=2700.0`: (Optional) Initial density (kg/m³).

# Returns
- `NamedTuple`: Constitutive model parameters, including elastic, plastic, and nonlocal properties.

# Example
```julia
cmp = setup_cmpr(mesh, instr; E=2.0e6, ν=0.25, ρ0=2500.0)
println(cmp.Kc)  # Bulk modulus
```
"""
function setup_cmpr(mesh::Mesh{T1,T2}; E::T2=T2(1.0e6),ν::T2=T2(0.3),ρ0::T2= T2(2700.0)) where {T1,T2}
    # independant physical constant          
    K,G,Del = D(E,ν,mesh.dim)                                                   # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    c       = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    Hp      = -30.0e3                                                           # softening modulus
    # constitutive model param.
    cmp = (;
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
    )
    return cmp::NamedTuple
end