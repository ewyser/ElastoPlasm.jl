function get_version()
    return string(Pkg.project().version)
end
function get_vals(mesh,mp,it,ηmax,ηtot,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("nel,np",(round(Int64,mesh.nel[1]*mesh.nel[2]),mp.nmp)),
            ("iteration(s)",it),
            ("ηmax,ηtot",(ηmax,ηtot)),
            (symb*" t/T",cmpl)]
    return vals
end
function msg(message)
    message = "└ "*message
    try
        return printstyled(message,color=:red,bold=true,blink=true)
    catch
        return printstyled(message,color=:blink)
    end
end
