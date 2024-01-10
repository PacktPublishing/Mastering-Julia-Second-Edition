# Code Contents for Chapter 10

#### The following files are present in this code directory.

**Chp10.jl** -- Listings for all code which can be run in the REPL or uploaded to VS Code.

**Chp10.ipynb** -- Jupyter notebook containing worked listings.

**echos.jl** -- Standalone echo server

**echoargs.jl** -- Script to cho arguments from the command line.

---

### TOML files and Setup script

**Project.toml** - Project TOML file.

**Manifest.toml** - Manifest TOML file.

**setup.jl** - Installation script for required Julia modules.

Use the *setup.jl* script to create new project files, deleting any existing ones.

This can either be run from the command line as *julia setup.jl* or using *include("setup.jl")* in the REPL.

The first time the directory is activated it is advised to <u>instantiate</u> it to bring the installed packages uptodate with the current version of Julia.
