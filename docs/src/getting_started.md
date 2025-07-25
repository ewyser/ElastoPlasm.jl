## Getting started

### How to ```plasmazing``` ?  

0. (opt.) Get Julia [here](https://julialang.org/downloads/) and follow instructions for installation
1. ```cd``` to the local repo ```./elastoPlasm ``` 
2. Launch Julia (on macOS, drag & drop ```start_macOS.sh``` in the terminal)
```julia
manuwyser@mBp elastoPlasm.jl % julia --project=. --threads=1
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.5 (2021-12-19)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```
3. Enter pkg mode ``` ] ```, then ```activate .``` the project ```elastoPlasm``` and ```instantiate``` its environment and related packages. You can ```st``` to check package status.
```julia
(elastoPlasm) pkg> st
Project elastoPlasm v0.3.5
Status `~/Dropbox/Jobs/git/elastoPlasm.jl/Project.toml`
  [6e4b80f9] BenchmarkTools v1.5.0
⌃ [052768ef] CUDA v5.4.3
  [f67ccb44] HDF5 v0.17.2
⌃ [63c18a36] KernelAbstractions v0.9.22
  [b964fa9f] LaTeXStrings v1.3.1
⌃ [91a5bcdd] Plots v1.40.5
  [92933f4c] ProgressMeter v1.10.2
  [295af30f] Revise v3.6.0
  [37e2e46d] LinearAlgebra
  [44cfe95a] Pkg v1.9.2
  [3fa0cd96] REPL
  [9a3f8284] Random
  [2f01184e] SparseArrays
Info Packages marked with ⌃ have new versions available and may be upgradable.

(elastoPlasm) pkg> 

```
4. Once ```elastoPlasm``` has been correctly instantiated, you can ```using elastoPlasm```

```juliaREPL
julia> using elastoPlasm
[ Info: Precompiling elastoPlasm [ff84b311-7f94-4593-ba3b-b1d6d2226836]
┌ Info: elastoPlasm: sucessful superInclude()
│ 	✓ init/scripts
│ 	✓ init/misc
│ 	✓ init/fun
└ 	✓ init/api

julia> names(elastoPlasm)
6-element Vector{Symbol}:
 :allocCheck
 :config_plot
 :e2eTest
 :elastoPlasm
 :slump
 :ϵp23De!

julia> 
julia> 
```

5. Input parameters: ```L``` is a vector of the Eulerian mesh dimensions $\ell_{x,z|x,y,z}$[^1], the argument ```nel``` defines the elements along the $x$ dim., ```varPlot``` defines state variable to be displayed in plot (```"P"``` for pressure, ```"du"``` for displacement or ```"epII"``` for plastic strain), and ```cmType``` defines the constitutive model being used. 

6. Optional ```kwargs``` are: ```shpfun={:smpm,:gimpm,:bsmpm}``` defining shape functions, ```fwrk={:infinitesimal,:finite}``` defining the deformation framework, ```trsf={:picflipUSL,:mUSL,:tpicUSL}``` defining the mapping scheme, ```vollock={true|false}``` is a boolean controlling volumetric locking corrections using the $\Delta\bar{F}$ method (see [1,2]), and ```GRF={true|false}``` is a boolean switching on/off gaussian random cohesion fields (see [3,4]). An illustrative example with ```kwargs``` is given below
```julia
ϵp23De([64.1584,12.80],40,"P","MC";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true)
```

7. Outputs (figs, gif, etc.) are saved in the folder ```./docs/out/ ```

8. (**final note**): ```ϵp23De``` automatically infers whether the two- or three-dimensional solver should be compiled, solely based on the dimension of ```L```, *i.e.,* ```length(L)={2|3}```.

[^1]: Two- and three-dimensional geometries are now supported since ```v0.3.0```.