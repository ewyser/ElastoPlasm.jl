abstract type AbstractSolver{T1<:Integer,T2<:Real,D} end

abstract type AbstractTime end





abstract type AbstractGeometry end

abstract type AbstractBasis end

abstract type AbstractEulerian end
abstract type CartesianMesh{T1, T2}  <: AbstractEulerian end
abstract type UniformMesh{T1, T2}    <: CartesianMesh{T1, T2} end
abstract type NonUniformMesh{T1, T2} <: CartesianMesh{T1, T2} end
abstract type MeshPhase{T1, T2}      <: CartesianMesh{T1,T2} end

abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end
abstract type MaterialPointPhase{T1, T2} <: MaterialPoint{T1,T2} end
abstract type AbstractElasticity{T1, T2} end
abstract type AbstractRheology{T1, T2} end