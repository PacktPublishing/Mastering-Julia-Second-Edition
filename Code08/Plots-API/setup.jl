#!/usr/bin/env julia
#
# Note: Assume MacOS/Linux
#

# Delete any existing TOML files

println("\nDeleting any TOML files ...")
rm("./Project.toml",force=true)
rm("./Manifest.toml",force=true)

# Create new TOML files

isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))
const PKGS = ["Plots","GR","IJulia"]

using Pkg
Pkg.activate(".")

println("Installing required packages ...")
for p in PKGS
  isinstalled(p) || Pkg.add(p)
end

println("\nListing installed packages ...")
Pkg.status()
