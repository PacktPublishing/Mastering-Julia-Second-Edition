#!/usr/local/bin/julia
#
isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))
const PKGS = ["CSV","DelimitedFiles","JLD","HDF5","LightXML","Query","Dates",
      "Statistics","StatsBase","DataFrames","TimeSeries","Plots","PyPlot","XLSX",
      "RDatasets","KernelDensity","HypothesisTests","GLM"]

using Pkg
Pkg.activate(".")

println("Installing required packages ...")
for p in PKGS
  isinstalled(p) || Pkg.add(p)
end
println("Done!")
