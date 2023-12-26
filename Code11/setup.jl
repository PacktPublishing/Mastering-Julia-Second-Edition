#!/usr/local/bin/julia
#
isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))

const PKGS = ["ArgParse","Getopt","BenchmarkTools","Example","Plots","Statistics","Traceur","Debugger","Revise","StatProfilerHTML"]

using Pkg
println("Installing required packages ...")

for p in PKGS
  isinstalled(p) || Pkg.add(p)
end
println("Done!")
