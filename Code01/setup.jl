#!/usr/local/bin/julia
#
# Note: Assume OSX and an alias setup to julia in /usr/local/bin
#
isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))

const PKGS = ["BenchmarkTools","UnicodePlots","PyCall","PyPlot","IJulia","Pluto","PlutoUI"]

using Pkg
println("Installing required packages ...")
for p in PKGS
  isinstalled(p) || Pkg.add(p)
end
println("Done!")
