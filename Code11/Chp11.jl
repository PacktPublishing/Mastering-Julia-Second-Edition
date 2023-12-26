[ SAMPLE Config file ]

#= Packages which I usually need and are in Base so add them here. =#
using Pkg, Printf, Random, REPLHistory

#= I use this to check for installed packages in my setup.jl, 
it is a little more verbose since the haskey() function was deprecated =#  

 isinstalled(pkg::String) =
  any(x -> x.name == pkg && x.is_direct_dep,
                     values(Pkg.dependencies()))
#= Julia has some pseudo Unix commands such as cd, pwd, 
   mkdir, rm, but these are missing – these can, of course,
   be run in shell mode =#
ls(dir::String = ".") = run(`ls -l`)
ls(dir::String = ".") = run(`ls -l $dir`)
cat(fname::String )   = run(`cat $fname`)
more(fname::String )  = run(`more $fname`)
# A little superfluous perhaps 
pwup(x, p::Number)  = x*p
pwup(x::Array)  = x.*p
sq(x)     = pwup(x, 2)
sqroot(x) = pwup(x, 0.5)
cb(x)     = pwup(x, 3)
cbroot(x) = pwup(x, 1//3)
# Defined in the book, so I have aded them here.
mad2(a,b,c) = a*b + c
systime()   = ccall((:time, "libc"), Int32, ());
# A macro to activate and instantiate the current directory
macro PkgSetup()
  Pkg.activate(".")
  Pkg.instantiate()
end
# … and one to execute an ‘immediate’ if statement. 
macro iif(cond,doit)
  if (eval(cond)) doit end
end
# I miss this older statement, so defined it here.
function linspace(x0::Real, x1::Real, N::Integer)
  @assert (x0 < x1) && (N > 0)
  h = (x1 - x0) / N
  return collect(x0:h:x1)
end
#= Another way to activate folders, automatically when 
   Julia starts by defining a new environment variable
   in the user’s shell start up file =#
if haskey(ENV,"JULIA_ACTIVATE") &&
        first(uppercase(ENV["JULIA_ACTIVATE"])) == 'Y'
  Pkg.activate(".")
  Pkg.instantiate()
end



# Juliet alias
juliet="julia -i -e 'import Pkg; Pkg.activate(\".\")'"


[ ftop.jl script]

#! /usr/bin/env julia --quiet --depwarn=no
#
const USAGE = "usage: $PROGRAM_FILE -h -d dir -n nfiles"
let
  nargs = size(ARGS)[1]
  dir = pwd()
  nf = 10
  hflag = false
  if nargs > 0
    i = 0
    while i < nargs
      i += 1
      s = ARGS[i]
      if s == "-h" 
        hflag = true
      elseif s == "-d"
        (i < nargs) && begin
          i += 1
          dir = ARGS[i]
        end
      elseif s == "-n"
        (i < nargs) && begin
          i += 1
          nf = ARGS[i]
        end
     end
    end
  end
  if hflag
    println(USAGE)
  else
    efind = 
      `find $dir -type f -iname "*" -exec du -sh "{}" + `
    println("Directory: $dir")
    try
      run(pipeline(efind,`sort -rh`,`head -$nf`))
      println("Done.")
    catch
      println("No files found.")
    end
  end
end


[Getopt.jl script]
#! /usr/bin/env julia
#
using Getopt
short_list = "ho:q"
long_list  = ["help", "quiet", "output="]
println("\nSize of original ARGS array: ", size(ARGS))
println("Short list :: ", short_list)
println("Long list  :: ", long_list, "\n")
for (opt, arg) in getopt(ARGS, short_list, long_list)
   @show (opt, arg)
end
println("Length of modified ARGS array: ", size(ARGS))


[Create Gadfly system image and use it to retrieve iris dataset and displayit]

using PackageCompiler
import Pkg; Pkg.activate(".")
Pkg.add("Gadfly");
Pkg.add("RDatasets")
create_sysimage([:Gadfly, :RDatasets];
                sysimage_path="sys_Gadfly.dylib")
                 
julia> using Gadfly, RDatasets
#= 
To create the plot we are required to push the Gadfly backend
into the O/S memory stack
=# 

pushdisplay(Gadfly.GadflyDisplay());
iris = dataset("datasets", "iris"); first(iris,4)

d = plot(iris, x=:SepalLength, y=:SepalWidth, color=:Species, Geom.point);
display(d)


[Panto.jl source]

"""
Calculate the sum of the series i/(i+1)^2 using the genie function for an integer 
i in the range [1:n].
"""
function aladdin(n::Integer)
  @assert n > 0
  s = 0.0
  for i in 1:n
    s += genie(i,2)
  end
  return s
end
"""
Compute the value of the expression x/(x+1)^k where x is a numeric and k is a 
(non-complex) number.
""" 
genie(x,k) = x/(x+1)^k
const N_THIEVES = 40; # Ali Baba has 40 thieves in the fairy story.
"""
Compute and store the items using Aladdin's genie storing each partial sum of the 
series i/(i+1)^2 in an array upto a value of 40.

So the preferable count to choose is a multiple of 40 in order to capture the 
final value sum of the series in the array.
"""
function alibaba(n::Integer)
   @assert n >= N_THIEVES
   s = 0.0
   k = n/N_THIEVES
#= 
We could define a fixed array but for only 41 items in since the Ali Baba band 
has 40 thieves, so push values instead 
=#
   a = Float64[]
   push!(a,s)
   for i in 1:n
     s += genie(i,2)
     (mod(i, k) == 0) && push!(a,s)
   end
   return a
end

# panto() function just directs the users to the contents and their help. 
function panto()
  helptxt = """Help is available on each of the individual pantomime characters: 
alibabi, aladdin, genie."""
  println(helptxt);
end

# Test out the alibaba function
y = alibaba(40);
round(y[end],digits=6)

# Time it with a larger number of trials
using BenchmarkTools
@benchmark alibaba(10^6)

# Compute values of sum for 40, 400, 4000, 40000 trials
for i in 1:4
  k = 4*10^i
  X = alibaba(k)
 @show (k, round(X[end], digits=5))
end

# Calcuate timing for one large number of trials
@benchmark(alibaba(4*10^9), samples=1, evals=1)


[Traceur]

# Modify alibabi to highlight step with uses most computing time

nn_thieves = 40;
julia> function alibaba_1(n::Integer)
  @assert n >= nn_thieves
  s = 0.0
  k = n/nn_thieves
  a = Float64[]
  push!(a,s)
  for i in 1:n
    s += genie(i,2)
    (mod(i, k) == 0) && push!(a,s)
  end
  return a
end

using Traceur
@trace alibaba_1(40);


# Use StatProfile to determin where bottle next is and redo alibabi
# This is later in the chapter, assume the panto function are included as separate file

using BenchmarkTools, StatProfilerHTML

include("panto.jl")
@profilehtml alibaba(10^6);

function alibaba_2(n::Integer)
  @assert n > 0   
  s = 0.0
  a = Array{Float64}(undef, 41)
  a[1] = s
  for i in 0:39
    for j = 1:n
      k = i*n + j
      s += genie(k,2)
    end
    a[i+2] = s
  end
  return a
end

# Check both version give same results and then time them

b1 = alibaba(120)
b2 = alibaba_2(3);   # Argument is 3 since it is multiplied by 40

@assert b1[end] == b2[end]

@btime b1 = alibaba(40*10^5);
@btime b2 = alibaba_2(10^5);

[ Source for the Funky module]

module Funky
using HTTP, CSV, DataFrames, TimeSeries, IndexedTables
using Printf, Dates, Statistics
const PUNCTS = 
       [' ','\n','\t',','.',',',':',';','!','?','\'','"'];
include("ftop.jl")
include("queens.jl")
include("quandl.jl")
include("panto.jl")

"""
basel(N::Integer)
Compute the sum of the Basel sequence
Returns a real number for positive values of N.
The solution was finally derived by Euler in 1734 and is equal to π^2 /6.
# Examples
julia> basel(10^6)
1.64493306684877
"""
function basel(N::Integer)
  @assert N > 0
  s = 0.0
  for i = 1:N
    s += 1.0/float(i)^2
  end
  return s
end

basel_c = ccall((:basel,"libmyfuns"), Cdouble, (Clong,), n)

horner_c = ccall((:horner,"libmyfuns.dylib"), Cdouble,
                (Cdouble, Ptr{Cdouble}, Clong), x, aa, n)

macro bmk(fex, n::Integer)
  quote
    let s = 0.0
    if $(esc(n)) > 0
      val = $(esc(fex))
      for i = 1:$(esc(n))
        local t0 = Base.time_ns()
        local val = $(esc(fex))
        s += Base.time_ns() - t0
      end
      return s/($(esc(n)) * 10e9)
    else
      Base.error("Number of trials must be positive")
    end
  end
end

export basel, fac, fib, hailstone, Kempner
export isdate, wordcount, filter, ftop, quandl
export @traprun, @bmk, @until

end












