export elastoplastic!,elastodynamic!,thermodynamic!              

"""
    elastodynamic!(mpts::Point{T1,T2,E,R}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle, updating material points and mesh.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastodynamic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D,S}) where {T1,T2,D,E,R,S<:DynamicRelaxationSolver}
    @warn "Elastodynamic workflow is not implemented for DynamicRelaxationSolver. Please use ExplicitSolver instead."
    return nothing
end  
"""
    elastoplastic!(mpts::Point{T1,T2,E,R}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastoplastic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle with elastoplastic update.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastoplastic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D,S}) where {T1,T2,D,E,R,S<:DynamicRelaxationSolver}
    @warn "Elastoplastic workflow is not implemented for DynamicRelaxationSolver. Please use ExplicitSolver instead."
    return nothing
end
