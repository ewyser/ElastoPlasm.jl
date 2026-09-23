# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Voigt-space vol/dev projection constants (dimension-generic, stateless — not stored per-particle state)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export StressStrainStiffness

struct StressStrainStiffness{NSTR,T}
    vol::SMatrix{NSTR,NSTR,T}
    dev::SMatrix{NSTR,NSTR,T}
    function StressStrainStiffness{3,T}() where {T}
        vol = SMatrix{3,3,T}(
            1, 1, 0,
            1, 1, 0,
            0, 0, 0
        )
        dev = SMatrix{3,3,T}(
             2/3, -1/3, 0.0,
            -1/3,  2/3, 0.0,
             0.0,  0.0, 0.5
        )
        new{3,T}(vol, dev)
    end
    function StressStrainStiffness{6,T}() where {T}
        vol = SMatrix{6,6,T}(
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        dev = SMatrix{6,6,T}(
             2/3, -1/3, -1/3, 0.0, 0.0, 0.0,
            -1/3,  2/3, -1/3, 0.0, 0.0, 0.0,
            -1/3, -1/3,  2/3, 0.0, 0.0, 0.0,
             0.0,  0.0,  0.0, 0.5, 0.0, 0.0,
             0.0,  0.0,  0.0, 0.0, 0.5, 0.0,
             0.0,  0.0,  0.0, 0.0, 0.0, 0.5
        )
        new{6,T}(vol, dev)
    end
end
# sig = Kc * (Del.vol * eps) + 2 * Gc * (Del.dev * eps), where eps = [exx,eyy,2exy] in 2D or [exx,eyy,ezz,2exy,2eyz,2exz] in 3D
