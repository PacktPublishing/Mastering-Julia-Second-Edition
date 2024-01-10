# Code Contents for Chapter 2

#### The following files are present in this code directory.

**Chp02.jl** -- Listings for all code which can be run in the REPL or uploaded to VS Code.

**Chp02.ipynb** -- Jupyter notebook containing worked listings.

**JSet.pluto.jl** -- Pluto workbook to display Juliaset without saving to disk.

**jmain.jl** -- Standalone scrupt to create Juliaset

_( jset.jl, pgmfile.jl files required by jmain.jl )_

---

### TOML files and Setup script
 
**Project.toml** - Project TOML file.

**Manifest.toml** - Manifest TOML file.

**setup.jl** - Installation script for required Julia modules.

Use the *setup.jl* script to create new project files, deleting any existing ones.

This can either be run from the command line as *julia setup.jl* or using *include("setup.jl")* in the REPL.

The first time the directory is activated it is advised to <u>instantiate</u> it to bring the installed packages uptodate with the current version of Julia.


