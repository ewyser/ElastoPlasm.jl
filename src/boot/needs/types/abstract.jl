#=
Supertype Problem definition
=#

abstract type AbstractGeometry end

abstract type AbstractEulerian end
abstract type AbstractCartesianMesh{T1, T2}  <: AbstractEulerian end
abstract type AbstractMeshPhase{T1, T2}      <: AbstractCartesianMesh{T1,T2} end

abstract type AbstractLagrangian end
abstract type AbstractMaterialPoint{T1, T2} <: AbstractLagrangian end
abstract type AbstractMaterialPointPhase{T1, T2} <: AbstractMaterialPoint{T1,T2} end
abstract type AbstractElasticity{T1, T2} end

abstract type AbstractConstitutiveModel{T2, D} end

abstract type AbstractTime end

abstract type AbstractProblem{T1,T2,D,SP<:AbstractMaterialPointPhase{T1,T2}} end


#=
Supertype Basis and Solver definition
=#

abstract type AbstractBasis end

abstract type AbstractSolver{T1<:Integer,T2<:Real,D} end
