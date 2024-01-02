#!/usr/bin/env julia
#
# Note: Script to work with MacOS or Linux
#
isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))

const PKGS = ["Images", "ImageView", "Colors", "ImageFiltering", "MosaicViews",
              "ImageEdgeDetection", "ImageTransformations"]

using Pkg
println("Installing required packages ...")
for p in PKGS
  isinstalled(p) || Pkg.add(p)
end
println("Done!")
