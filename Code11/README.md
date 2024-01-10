# Code Contents for Chapter 11
#### The following files are present in this code directory.

**Chp11.jl** -- Listings for all code which can be run in the REPL or uploaded to VS Code.

**ftop.jl** -- Stand-alone script for finding largest files in a folder and its subfolders.

**getopt.jl** -- Script for testing option passing on the command line.

**panto.jl** -- File of the _Pantomime_ routines. 

**queens.jl** -- Function to solve the _Queens_ problem.

---

**Funky** -- Folder for the Funky project.

**statprof** -- Folder created by Stat HTML Profiler

---
### TOML files and Setup script

**Project.toml** - Project TOML file.

**Manifest.toml** - Manifest TOML file.

**setup.jl** - Installation script for required Julia modules.

Use the *setup.jl* script to create new project files, deleting any existing ones.

This can either be run from the command line as *julia setup.jl* or using *include("setup.jl")* in the REPL.

The first time the directory is activated it is advised to <u>instantiate</u> it to bring the installed packages uptodate with the current version of Julia.
