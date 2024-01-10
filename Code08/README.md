# Code Contents for Chapter 8

#### The following files are present in this code directory.

**Chp08.jl** -- Listings for basic examples which can be run in the REPL or uploaded to VS Code.

**Chp08.ipynb** -- Jupyter notebook containing selected worked listings.

---

### Folders with further more advanced examples

- [Plots API](Plots-API)
- [Gadfly](Gadfly) 
- [Plotly](Plotly)
- [Statistics Plots](StatPlots)
- [Makie Framework](Makie)
- [Images Package](ImageProcs)

---

### TOML files and Setup script
 
**Project.toml** - Project TOML file.

**Manifest.toml** - Manifest TOML file.

**setup.jl** - Installation script for basic Julia examples.

Use the *setup.jl* script to create new project files, deleting any existing ones.

This can either be run from the command line as *julia setup.jl* or using *include("setup.jl")* in the REPL.

The first time the directory is activated it is advised to <u>instantiate</u> it to bring the installed packages uptodate with the current version of Julia.
