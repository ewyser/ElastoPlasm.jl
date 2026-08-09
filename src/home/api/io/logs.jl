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
        printstyled("│", color=:green, bold=true); println("   L,nel  = [64.1584,64.1584/4.0],[40,10];")
        printstyled("│", color=:green, bold=true); println("   jld2   = ic_slump(L,nel;cli()...);")
        printstyled("└", color=:green, bold=true); println("   out    = elastoplasm(jld2; workflow = [elastodynamic!,elastoplastic!]);")
    elseif showcase == "collapsing"
        printstyled("│", color=:green, bold=true); println("   plot    = (; status=true, freq=1.0, what=[\"sigxx\"], dims=(500.0,250.0) );")
        printstyled("│", color=:green, bold=true); println("   nel     = [5,10];")
        printstyled("│", color=:green, bold=true); println("   ic, cfg = ic_collapse(nel,0.0,1.0e4,80.0,50.0; plot);")
        printstyled("└", color=:green, bold=true); println("   out     = collapse(ic, cfg);")
    elseif showcase == "heating"
        printstyled("│", color=:green, bold=true); println("   ic, cfg = ic_thermal();")
        printstyled("└", color=:green, bold=true); println("   out     = thermal(ic, cfg);")
    
    else
        printstyled("└", color=:green, bold=true); println("   ...$(showcase) ?!? \e[5m¯\\_(ツ)_/¯\e[0m")
    end

    return nothing
end

"""
    ic_log(mesh, mpts, time) -> String

Generates a summary log string for the initial simulation setup, including mesh, material points, and time parameters.

# Arguments
- `mesh`: Mesh object containing element information.
- `mpts`: Material point object containing number of points.
- `time`: Named tuple with time parameters (`t`, `te`, `tg`, `tep`).

# Returns
- `String`: Multi-line summary of the simulation setup.

# Example
```julia
summary = ic_log(mesh, mpts, time)
println(summary)
```
"""
function ic_log(mesh::Mesh,mpts::Point,time::Time,solver::S) where {T1<:Integer,T2<:Real,D<:AbstractDimension, S<:AbstractSolver{T1,T2,D}}
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
    elastoplasm_log(instr::NamedTuple; msg::String="elastodynamic") -> String

Generates a summary log string describing the current simulation configuration for ElastoPlasm.

# Arguments
- `instr::NamedTuple`: Instruction/configuration named tuple containing simulation options.
- `msg::String`: (Optional) Workflow name, default is "elastodynamic".

# Returns
- `String`: Multi-line string summarizing the simulation setup.

# Example
```julia
logstr = elastoplasm_log(instr)
println(logstr)
```
"""
function elastoplasm_log(solver::S; msg::String="elastodynamic") where {T1<:Integer,T2<:Real, D<:AbstractDimension, S<:AbstractSolver{T1,T2,D}}
    # build the list of log lines
    logs = [
        "Launching ϵlastσPlasm 👻 v$(get_version()):",
        "└ $(nthreads()) active thread(s)",
        "- solver: $(typeof(solver))",
        "- $(solver.fwrk.deform) strain formulation",
        "- $(solver.fwrk.trsfr) mapping scheme",
        "- $(solver.basis.which) calculation cycle",
    ]
    # add optional lines only if the corresponding flags are true
    if solver.fwrk.locking
        push!(logs, "- F-bar locking mitigation")
    end
    if solver.nonloc.status && occursin("plastic",msg)
        push!(logs, "- non-local plastic regularization")
    end
    return join(logs,"\n")
end

"""
    exit_log(message::String)

Print a styled exit message to the console, using green, bold, and blinking text if supported by the terminal.

# Arguments
- `message::String`: The message to display at program exit.

# Returns
- `Nothing`. Prints the message to the console.

# Example
```julia
exit_log("Simulation finished successfully.")
```
"""
function exit_log(message::String)
    try
        return printstyled("└ $message",color=:green,bold=true,blink=true)
    catch
        return printstyled("└ $message",color=:blink)
    end
end