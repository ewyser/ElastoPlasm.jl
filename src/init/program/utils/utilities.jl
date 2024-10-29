function getVersion()
    return string(Pkg.project().version)
end
function getvals(meD,mpD,it,ηmax,ηtot,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("nel,np",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),
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
