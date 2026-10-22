"""
    welcome_log(; greeting::String="Welcome to ϵlastσPlasm 👻 v$(get_version())")

Prints a styled welcome message to the console, highlighting "Welcome" and vertical bars in green and bold.

# Arguments
- `greeting::String`: The greeting message to display at the top (default: "Welcome to ϵlastσPlasm 👻 v$(get_version())").

# Returns
- `Nothing`. Prints the welcome message to the console.

# Example
```julia
welcome_log()
welcome_log(greeting="Hello from ElastoPlasm!")
```
"""
function welcome_log(; greeting::String="Welcome to ϵlastσPlasm 👻 v$(get_version())", showcase::String = "slumping") 
    printstyled("┌ $greeting\n", color=:green, bold=true)
    printstyled("│", color=:green, bold=true); println(" New comer ? Try $(showcase) out")
    if showcase == "slumping"
        printstyled("│", color=:green, bold=true); println("   L,nel = [64.1584,64.1584/4.0],[40,10];")
        printstyled("│", color=:green, bold=true); println("   jld2  = setup_simulation(slump_case(L,nel;cli()...));")
        printstyled("└", color=:green, bold=true); println("   execute!(jld2);")
    elseif showcase == "heating"
        printstyled("│", color=:green, bold=true); println("   L,nel = [16.0,16.0],[80,80];")
        printstyled("│", color=:green, bold=true); println("   jld2  = setup_simulation(thermal_case(L,nel));")
        printstyled("└", color=:green, bold=true); println("   execute!(jld2; workflows = [thermodynamic!]);")
    elseif showcase == "collapsing"
        printstyled("│", color=:green, bold=true); println("   nel  = [92,22];")
        printstyled("│", color=:green, bold=true); println("   jld2 = setup_simulation(collapse_case(nel));")
        printstyled("└", color=:green, bold=true); println("   execute!(jld2);")
    else
        printstyled("└", color=:green, bold=true); println("   ...$(showcase) ?!? \e[5m¯\\_(ツ)_/¯\e[0m")
    end

    return nothing
end

"""
    ic_log(mesh, mpts, time, solver) -> String

The summary `setup_simulation` logs: precision, elements, material points and the time stages
(`time.t`, gravity ramp `tg`, elastic stage `te`, plastic stage `tep`).
"""
function ic_log(mesh::Mesh,mpts::Point,time::Time,solver::S) where {T1<:Integer,T2<:Real,D, S<:AbstractSolver{T1,T2,D}}
    # build the list of constant log lines
    logs = [
        "Summary: ",
        "- $(solver.dtype.precision)",
        "- elements: $(mesh.prprt.nel[end])",
        "- material points: $(mpts.nmp)", 
        "- simulation time ∈ $(time.t) s:",
    ]
    # add optional lines
    if isa(time.tg,AbstractFloat)
        push!(logs, "   - gravity ramp-up: $(time.tg ) s")
    end
    if isa(time.te,AbstractFloat)
        push!(logs, "   - elastodynamic  : $(time.te ) s")
    end
    if isa(time.tep,AbstractFloat)
        push!(logs, "   - elastoplastic  : $(time.tep) s")
    end
    return join(logs,"\n")::String
end

"""
    workflow_log(solver; msg="elastodynamic") -> String

What `execute!` logs before each workflow (`msg`, its name): version, threads, solution, elastic law,
transfer scheme, basis, and F-bar or non-local regularization when on (the latter for plastic workflows).
"""
function workflow_log(solver::S; msg::String="elastodynamic") where {T1<:Integer,T2<:Real, D, S<:AbstractSolver{T1,T2,D}}
    # build the list of log lines
    logs = [
        "Launching ϵlastσPlasm 👻 v$(get_version()):",
        "└ $(nthreads()) active thread(s)",
        "- solver: $(solver.solution)",
        "- $(solver.material.elastic) material",
        "- $(solver.basis.trsfr) mapping scheme",
        "- $(solver.basis.which) calculation cycle",
    ]
    # add optional lines only if the corresponding flags are true
    if solver.stab.locking
        push!(logs, "- F-bar locking mitigation")
    end
    if solver.nonloc.status && occursin("plastic",msg)
        push!(logs, "- non-local plastic regularization")
    end
    return join(logs,"\n")
end
