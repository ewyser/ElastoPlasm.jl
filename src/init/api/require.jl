"""
    require(in::Symbol)

Generates and returns a Dictionnary of reference instruction

# Examples

```jldoc
julia> require(:instr)
Dict{Symbol, Any} with 10 entries:
  :dtype   => Float64
  :nonloc  => (cond = false, ls = 1.0)
  :GRF     => false
  :shpfun  => :bsmpm
  :trsfr   => :mUSL
  :vollock => true
  :plast   => (false, "DP")
  :plot    => (cond = true, freq = 1.0, what = ["epII"])
  :perf    => true
  :fwrk    => :finite

julia> 
```

# Further specification
Admissible keywords, by-default value and purpose are presented below when ```instr = require(:instr)```
- ```:dtype   ```, # set the arithmetic precision
- ```:basis   ```, # define shapefunction type
- ```:fwrk    ```, # set the deformation framework
- ```:trsfr   ```, # set the mapping scheme
- ```:vollock ```, # set volumetric locking mitigation strategy
- ```:GRF     ```, # activate Gaussian Random Field generator
- ```:plast   ```, # set plasticity onset and plastic flow law
- ```:nonloc  ```, # set non-local regularization
- ```:plot    ```, # option for plotting capabilities
- ```:perf    ```, # set performance mode

"""
function require(in::Symbol=:instr)
    if in == :instr
        instr = Dict(
            :bits    => 64,
            :basis   => (;which="bsmpm",
                        how=nothing,
                        ghost=false,
                        ),
            :fwrk    => "finite",
            :trsfr   => "mUSL",
            :vollock => true,
            :GRF     => false,
            :plast   => (false,"DP"),
            :nonloc  => (;cond=true,
                        ls=0.5,),
            :plot    => (;cond=true,
                        freq=1.0,
                        what=["epII"],
                        dims=(500.0,250.0),
                        ),
            :perf    => true,
        )
        return instr
    else
        error("$(in) is an unsupported symbol")
        return nothing
    end
end
export require