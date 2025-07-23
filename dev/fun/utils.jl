function mdGenerate(file="function.md",path="docs/src")
    # concatenate in a single stirng the list all function declaration in api
    libs = Dict()
    for (k,lib) ∈ enumerate(keys(elastoPlasm.info.sys.lib))
        lib_list = undef
        for (l,val) ∈ enumerate(lib)
            if l == 1
                lib_list = "\"$(val)\","
            else
                lib_list = lib_list*"\"$(val)\","
            end
        end
        push!(libs,lib=>lib_list)
    end
    # check if path already exists
    if isfile(joinpath(path,file))
        rm(joinpath(path,file))
    end
    # open a .md file and writes corresponding doc
    open(joinpath(path,file),"w") do f
        content = 
        """
        ```@meta
        CollapsedDocStrings = true
        ```

        # Functions

        ## Public API
        ```@autodocs
            Modules = [elastoPlasm]
            Pages   = [$(libs["api"])]
        ```
        ## Internals
        ```@autodocs
            Modules = [elastoPlasm]
            Pages   = [$(libs["program"])]
        ```
        """
        write(f,content)
    end
    return file
end