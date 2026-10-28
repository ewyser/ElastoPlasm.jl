# EasyDraw (lib/EasyDraw), the package the shapes, zones and drawing window come from: its own tests, run here so
# that one run of this suite covers both
include(joinpath(dirname(dirname(@__DIR__)), "lib", "EasyDraw", "test", "runtests.jl"))
