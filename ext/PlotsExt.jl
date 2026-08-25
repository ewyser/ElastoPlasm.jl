__precompile__(false)

module PlotsExt

using ElastoPlasm,Plots,LaTeXStrings,Base.Threads,JLD2

#include(joinpath(dirname(@__FILE__),"PlotsExt","path-to-functions.jl"))


self.ui.plot = true

end