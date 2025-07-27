function get_version()
    return string(Pkg.project().version)
end
function get_vals(mesh,mp,it,ηmax,ηtot)
    # save vals
    vals = [("nel,np",(round(Int64,prod(mesh.nel[1:end-1])),mp.nmp)),
            ("iteration(s)",it),
            ("ηmax,ηtot",(ηmax,ηtot))]
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
