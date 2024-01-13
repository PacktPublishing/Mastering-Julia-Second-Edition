
using Pkg
Pkg.update()

ipks = Pkg.installed()
pkgs = ["PyCall","RCall","RDatasets","LibCURL","JavaCall","DistributedArrays","SymPy"]
for p in pkgs
    !haskey(ipks,p) && Pkg.add(p)
end

using PyCall, PyPlot, Plots, RCall, RDatasets, JavaCall 
using DistributedArrays, SymPy, LibCURL
using Printf, Random, Distributed

# Simple call to a C-routine
# Can be used to set the random number seed using the system clock

systime()   = ccall((:clock,"libc"),Int32,())
randomize() = Random.seed!(systime())

systime()

# A 'mad' example (a*b + c) 
# The routine is called fma in library: libc
#
mad(a,b,c) = ccall((:fma,"libc"),Float64,(Float64,Float64,Float64),a,b,c)
mad(3.1,5.2,7.4)

# Use the system library to generate some random numbers
[ccall((:rand, "libc"), Int32, ()) for i = 1:6]

# Call a FORTRAN routine from LAPACK to compute the dot product
# between two arrays.
# FORTRAN passes scalar arguments by reference
# Note: the vectors arepassed by reference already. 

function compute_dot(DX::Vector{Float64},DY::Vector{Float64})
    @assert length(DX) == length(DY)
    n = length(DX)
    incx = incy = 1
    dotprod = ccall((:ddot, "libLAPACK"),
                   Float64,
                   (Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, 
                   Ptr{Float64}, Ptr{Int64}),
                   Ref(n), DX, Ref(incx), DY, Ref(incy))
    return dotprod
end

# Test it out
aa = [rand() for i = 1:1000]
bb = [rand() for i = 1:1000]
compute_dot(aa,bb)    # => Should be about 250

# Get the user's home directory ENV["HOME"]
ptr = ccall((:getenv, "libc"), Ptr{UInt8},(Ptr{UInt8},),"HOME")

# My modules are in Julia/MyModules - put this on the LOAD PATH
#
myHome = unsafe_string(ptr)
push!(LOAD_PATH,string(myHome,"/Julia/MyModules"))

#= 
The system library call has a third parameter, which when set to zero will create a new variable but not overwrite an existing one.
=#
# Define the function to replace existing variable
evset(var::String, val::String) = 
     ccall((:setenv,"libc"),Clong,
           (Cstring,Cstring,Clong),var,val,1);
 
# The unset routine is quite simple
evunset(evvar::String) = 
     ccall((:unsetenv,"libc"),Clong,(Cstring,),evvar);

#= 
The system library call has a third parameter, which when set to zero will create a new variable but not overwrite an existing one.
=#
# Define the function to replace existing variable
evset(var::String, val::String) = 
     ccall((:setenv,"libc"),Clong,
           (Cstring,Cstring,Clong),var,val,1);
 
# The unset routine is quite simple
evunset(evvar::String) = 
     ccall((:unsetenv,"libc"),Clong,(Cstring,),evvar);


# Set an environment variable PACKT_HOME ...
# ... and check it
packt = evset("PACKT_HOME",
            string(myHome,"/Users/malcolm/PacktPub"));

julia> ENV["PACKT_HOME"]
"/Users/malcolm/PacktPub"

# Now unset it, verify it is so.
julia> evunset("PACKT_HOME");
julia> ENV["PACKT_HOME"]
ERROR: KeyError: key "PACKT_HOME" not found


#= 

// Basel function in C

#include<stdio.h>
#include<stdlib.h>

double basel(int N) {

  double s = 0.0L;
  int i;
  double x;

  if (N < 1) return s;
  for (i = 1; i <= N; i++) {
    x = 1.0L/((double) i);
    s += x*x;
  }
  return s;
}

// Horner's method

#include<math.h>

double horner(double x, double aa[], long n) {
  long i;
  double s = aa[n-1];
  if (n > 1) { 
    for (i = n-2; i >= 0; i--) {
      s = s*x + aa[i];
    }
  }
  return s;
}

// Build a dynamic library (on OSX) as:
//
// clang -c basel.c horner.c
// libtool -dynamic basel.o horner.o -o libmyfuns.dylib  
//         -lSystem -macosx_version_min 10.13
//
// sudo cp libmyfuns.dylib /usr/local/lib

=#

run(`nm -g libmyfuns.dylib`)

basel = ccall((:basel,"libmyfuns"),Float64,(Int64,),10000000)

# Time the function for 10^7 loops
@elapsed ccall((:basel,"libmyfuns"),Float64,(Int64,),10000000)

# C version of Horners method
# Note use of Cdouble, Clong rather than Float64, Int64 etc.
#
x = 2.1;
aa = [1.0, 2.0, 3.0, 4.0, 5.0];

ccall((:horner,"libmyfuns.dylib"),Cdouble,
              (Cdouble,Ptr{Cdouble},Clong),x,aa,length(aa))

# Compute PI in 'C', usual Monte Carlo method
# We ca put the C code in a multiline string

C_code = """
#include <stddef.h>
#include <stdlib.h>

double c_pi(long n) {
    long k = 0L;
    float rmax = (float) RAND_MAX;
    for (long i = 0L; i < n; i++) {
        float x = ((float) rand())/rmax;
        float y = ((float) rand())/rmax;
        if ((x*x + y*y) < 1.0) {
          k++;
        }
    }
    return 4.0*((double) k)/((double) n);
}
"""

# Get a temporary name to create a library
const Clib = tempname()   # ... make a temporary file

# compile to a shared library by piping C_code to gcc
# (works only if you have gcc installed):
using Libdl
tmplib = string(Clib,".",dlext)
open(`gcc -fPIC -O3 -msse3 -xc -shared -o $tmplib -`, "w") do f
    print(f, C_code) 
end

# define a Julia function that calls the C function:
c_pi(N::Int64) = ccall(("c_pi", Clib), Float64, (Clong,), N)

using Random
randomize();
c_pi(1000000)

include <julia.h>
#include <stdio.h>
#include <math.h>

// Only define this once if in an executable ...
// (i.e. not in a shared library) ...
// ... and if we want the fastest code.
JULIA_DEFINE_FAST_TLS()

int main(int argc, char *argv[])
{
/* required: setup the Julia context */
  jl_init();

/* run Julia commands */
  jl_function_t *fnc1 = jl_get_function(jl_base_module, "exp");
  jl_function_t *fnc2 = jl_get_function(jl_base_module, "sin");
  jl_value_t* arg1 = jl_box_float64(-0.3);
  jl_value_t* arg2 = jl_box_float64(3.0);
  jl_value_t* ret1 = jl_call1(fnc1, arg1);
  jl_value_t* ret2 = jl_call1(fnc2, arg2);

/* unbox and setup final result */
  double retD1 = jl_unbox_float64(ret1);
  double retD2 = jl_unbox_float64(ret2);
  double retD3 = retD1*retD2;
  printf("sin(3.0)*exp(-0.3) from Julia API: %e\n", retD3); 
  fflush(stdout);

/* Allow Julia time to cleanup pending write requests etc. */
  jl_atexit_hook(0);
  return 0;
}


# Can use a nice command script in 'share/julia/julia-config.jl
./julia-config.jl 
usage: julia-config [--cflags | --ldflags |
                     --ldlibs | --allflags]

JULIA_HOME = \
/Applications/Julia-1.0.app/Contents/Resources/julia; export JULIA_HOME

ls $JULIA_HOME
LICENSE.md bin etc include lib share


cc jltest.c -o jltest -std=gnu99 \
      -I$JULIA_HOME/include/julia \
      -DJULIA_ENABLE_THREADING=1 -fPIC \
      -L$JULIA_HOME/lib \
      -Wl,-rpath,$JULIA_HOME/lib \
      -Wl,-rpath,$JULIA_HOME/lib/julia \
      -ljulia

./jltest
sin(3.0)*exp(-0.3) from Julia API: 1.045443e-01

using LibCURL
url = "http://LondonJulia.org/mastering-julia.html"

# init a curl handle
curl = curl_easy_init();

# set the URL and request to follow redirects
curl_easy_setopt(curl, CURLOPT_URL, url);
curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);

# setup the callback function to recv data
function curl_write_cb(curlbuf::Ptr{Nothing}, s::Csize_t, n::Csize_t, p_ctxt::Ptr{Nothing})
  sz = s * n
  data = Array{UInt8}(undef,sz)
  ccall(:memcpy, Ptr{Nothing}, (Ptr{Nothing}, Ptr{Nothing}, UInt64), data, curlbuf, sz)
  println(String(data))
  sz::Csize_t
end

c_curl_write_cb = 
  @cfunction(curl_write_cb, Csize_t, 
             (Ptr{Nothing}, Csize_t, Csize_t, Ptr{Nothing}));
curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, c_curl_write_cb);

# execute the query
res = curl_easy_perform(curl)

println("curl url exec response : ", res)

# retrieve HTTP code
http_code = Array{Clong}(undef,1)
curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, http_code)
println("httpcode : ", http_code)

# release handle
curl_easy_cleanup(curl)

using PyCall

@pyimport numpy.random as nr
aa = nr.rand(4,5)      # aa is a Julia array generated with Numpy

# So an array slice has no overhead(i.e. in Julia)
aa[2:3,2:3]

# Using the SciPy package
@pyimport scipy.optimize as so
so.ridder(x -> x*cos(x),1,π)

@pyimport scipy.integrate as si
si.quad(x -> x*sin(x),1,π)

# ... and plotting with matplotlib
@pyimport matplotlib.pyplot as plt
x = range(0,stop=10,length=1000); 
y = sin.(3*x + 4*cos.(2*x));

plt.plot(x, y, color="red", linewidth=2.0, linestyle="--")

plt.show()

using SymPy

# We need to work with a special type: a symbol (in the SymPy sense)
u = symbols("u")
x = symbols("x", real=true)
y1, y2 = symbols("y1, y2", positive=true)
alpha = symbols("alpha", integer=true, positive=true)

typeof(x)

# We can solve some algebraic equations
solve(u^2 + 1)

# Perform functional expansions
p = expand(prod([sin(x^(-i)) for i in 1.0:1.0:5.0]))

# SymPy can use the Plot API to display graphs
x = symbols("x")
using Plots
pyplot()

plot(p, 0.0, 1.0)

# For a more complex example, which we have seen before, ...
# ... the sum of the Basel series
#
i, n = symbols("i, n")
sn = Sum(1/i^2, (i, 1, n))
doit(sn)

limit(doit(sn), n, oo)

# An initial valued problem
#
y = SymFunction("y")
a, x = symbols("a,x")
eqn = y'(x) - 3*x*y(x) - 1

# Solve this equation for x0 = 0 ...
# and currying the function y0 -> a 
#
x0, y0 = 0, a
out = dsolve(eqn, x, (y, x0, y0))

# Resolve the curried function by substuting a value for 'a'
out |> subs(a,2.1)

using RCall

# Perform a t-test on a set of 1000 normally distributed variates
#
x = randn(1000);
R"t.test($x)"

# Optimise function using BFGS (Broyden–Fletcher–Goldfarb–Shannon§§§§§≠±±) algorithm
#
f(x) = 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x + 40
R"optim(0, $(x -> f(x)), method='BFGS')"

@rlibrary MASS

 geyser = rcopy(R"MASS::geyser")  

round(sum(geyser[:waiting])/sum(geyser[:duration]),digits=5)

rcall(:summary, geyser[:waiting])

rcall(:summary, geyser[:duration])

# Use R to grab some financial data
# and plot it.

R"""
library(data.table)
library(scales)
library(ggplot2)
        
link <- "fin_data.csv"
dt <- data.table(read.csv(link))
dt[, date := as.Date(date)]
        
# create indexed values
dt[, idx_price := price/price[1], by = ticker]
"""

R"""
ggplot(dt, aes(x = date, y = idx_price, color = ticker)) + 
       geom_line() + theme_bw() + 
       xlab("Date") + ylab("Pricen(Indexed 2000 = 1)") + 
       scale_color_discrete(name = "Company")
"""

using JavaCall
JavaCall.init(["-Xmx128M"])     # Initialise and ask for additional memory

#=
 May set the current directory as the classpath in the init() call
 viz. JavaCall.init(["-Xmx128M -Djava.class.path=$(@__DIR__)"])

 Ignore the (possible) segmentation fault
 See: http://juliainterop.github.io/JavaCall.jl/faq.html

 This works inside in a Notebook, despite the segmentation fault (on OSX)
 In the REPL, julia may nned to be started with a --handle-signals=no
 option to disable Julia's signal handler. 
 [This may cause issues with handling ^C in Julia programs.]
=#

# JavaCall needs to modify Julia variables before using them in functions  
a = JString("Hello, Blue Eyes")

# The JString as a single field (:ptr)
fieldnames(typeof(a))

a.ptr

# The pointer can now be used in a Java function

b = ccall(JavaCall.jnifunc.GetStringUTFChars, Ptr{UInt8},
            (Ptr{JavaCall.JNIEnv}, Ptr{Nothing}, Ptr{Nothing}),
                                 JavaCall.penv, a.ptr, C_NULL)

unsafe_string(b)

# Similar to PyCall, we can use JavaCall to perform some mathematics

jlm = @jimport "java.lang.Math"
jcall(jlm, "exp", jdouble, (jdouble,), pi*pi/6.0)

# ... and to access Internet libraries
jnu = @jimport java.net.URL
jurl = jnu((JString,), "http://LondonJulia.org/mastering-julia.html")
www = jcall(jurl, "getHost", JString,())

# Finally let's import the ArrayList class, ...
# ... define an instance and add some items to it
JArrayList = @jimport(java.util.ArrayList)
a = JArrayList(())
jcall(a, "add", jboolean, (JObject,), "'Twas ")
jcall(a, "add", jboolean, (JObject,), "brillig, ")
jcall(a, "add", jboolean, (JObject,), "and ")
jcall(a, "add", jboolean, (JObject,), "the ")
jcall(a, "add", jboolean, (JObject,), "slithy ")
jcall(a, "add", jboolean, (JObject,), "toves,")

# Now iterate thru' the array and push it on a Julia array ...
# ... converting the bit type to an (unsafe) string.
#
t = Array{Any, 1}()
for i in JavaCall.iterator(a)
  push!(t, unsafe_string(i))
end

# Evaluate the  result.
join(t)


# Get a webpage using curl or wget
# Notice the backticks
# First check that curl (or wget) is available, ...
# ... otherwise you will need to install it. 
#
cmd = `which curl`
typeof(cmd)

#=
 Commands now need to be run, i.e. they are not executed immediately
 Since commands run as separate tasks, it is usually preferable to 
 suppress the output of the run() function

 The task output will goto STDOUT (by default) but can be captured
 and used by the Julia program
=#
run(cmd)

# We can use curl to get the webpage ...
# .. and can do it in one step
proc = run(`curl "http://LondonJulia.org/mastering-julia.html"`);

# Look at Base:process.jl for definition of Julia process structure

#=
mutable struct Process <: AbstractPipe
    cmd::Cmd
    handle::Ptr{Cvoid}
    in::IO
    out::IO
    err::IO
    exitcode::Int64
    termsignal::Int32
    exitnotify::Condition
    closenotify::Condition
    . . .
    . . .
end
=#

proc.exitcode   # Zero is normal exit status

cd(string(ENV["HOME"],"/PacktPub/Alice"))
pwd()  #  Needs to be where the chapter 5 files are

# Define a function to use the wordcount ('wc') utility
wc(f) = isfile(f) && run(`wc $f`)

# Count the number of occurences of 'beaver" in the Hunting of the Snark,
# piping the output to a file

txtfile = "hunting-the-snark.txt";
logfile = "hunting-the-snark.log";

run(pipeline(`grep -i beaver $txtfile`,stdout=logfile));

wc(logfile);

# Do this again for the bellman but append the output to the log file.

run(pipeline(`grep -i bellman $txtfile`,stdout=logfile,append=true));
wc(logfile);

#=
  Do we have any lines with both beaver and bellman
  Note in the 'grep' we ignored case (-i) but will stiil find
  plurals and possibly punctuation, viz.: Beavers.
  I'll attend to this later.
=#
run(pipeline(`grep -bellman $txtfile`,`grep -i beaver`));

function ftop(dir=".", nf=20)
  @assert nf > 0
  efind =  `find $dir -type f -iname "*" -exec du -sh "{}" + `
  try
    run(pipeline(efind,`sort -rh`,`head -$nf`))
    println("Done.")
  catch
    println("No files found.")
  end
end

dd = ENV["HOME"]*"/PacktPub"
ftop(dd,10)


# Code to run from the command line
# Source: ftop-main.jl

#=
# See if any arguments have been passed
# If so process them otherwise set the defaults
#
nargs = size(ARGS)[1]
if nargs > 2
  println("usage: ftop [dir [nfiles]]")
  exit(-1)
else
  dir = (nargs > 0) ? ARGS[1] : "."
  nf  = (nargs > 1) ? ARGS[2] : 20
end

# Main 'find' command
efind =  `find $dir -type f -iname "*" -exec du -sh "{}" + `

# Get and print the directory to be searched
cwd = string(pwd(),"/",dir)
println("Directory: $cwd")

# Then run the pipeline
try
  run(pipeline(efind,`sort -rh`,`head -$nf`))
  println("Done.")
catch
  println("No files found.")
end
=#



#=
 Perl has falled out of fashion with the rise of Python but still
 but still remains one of the best methods for data munging.

  Unix distros and OSX (normally) have Perl avaiable but in Windows
  it needs to installed and on the executable path.

  Julia performance in handling string is not one of its greatest strenghs
  so munging large files can successfully done using Perl.

  Note: 
  Julia has introduced an analytical engine (JuliaDB) to tackle the
  processing of large datasets, employing some clever memory management
  techniques and we will discuss this in the next chapter

=#

# On OSX and Linux, there is a word list in /usr/share/dict
# This command outputs an palindromes of 6 or more letters
# (I am not try to teach Perl but most of the syntax is easily followed)
cmd = `perl -nle 'print if $_ eq reverse && length > 5' /usr/share/dict/words`
run(cmd);

# Here is a work-around using temporary files
# Capture the output of palindrome command and 

tmpfile = mktemp(tempdir());   # Alternative to tempfile() call
cmd = `perl -nle 'print if $_ eq reverse && length > 5' /usr/share/dict/words`
run(pipeline(cmd; stdout=tmpfile)); 
dmp = read(tmpfile);
run(`rm -f $tmpfile`);  # Good idea to remove the temporary file

# Dump is a byte array
(eltype(dmp),length(dmp))

# Convert the byte array to a string
# Note the carriage returns 
ss = String(dmp)

# Remove the trailing \n and then split into words
sa = split(chomp(ss),"\n")
typeof(sa)

# PERL is good (and quick) for processing strings
# A large pipeline (probably not on Windows)
# Get top 10 words in the Hunting of the Snark, ignoring blank lines
# I'll leave it the the reader ti work out what the command do!
#
cd(ENV["HOME"]*"/PacktPub/Alice")

fl = "hunting-the-snark.txt";
c1 = `perl -pne 'tr/[A-Z]/[a-z]/' $(fl)`;
c2 = `perl -ne 'print join("\n", split(/\s+/,$_));print("\n")'`;

run(pipeline(c1,c2,`sort`,`grep -ve '^$'`,`uniq -c`,`sort -rn`,`head -10`));

# Look at occurences of Bellman
#
run(pipeline(c1,c2,`sort`,`grep -ve '^$'`,`uniq -c`,`sort -rn`,`grep -i Bellman`));

# Suspect there is a problem with punctuation marks
# Here are the bottom 10 words.
#
run(pipeline(c1,c2,`sort`,`grep -ve '^$'`,`uniq -c`,`sort -rn`,`tail -10`));

# So we need to add an extra task in the pipe
# If we do this through a function we can pass any poem
#
function munge(fl)
  c1 = `perl -pne 'tr/[A-Z]/[a-z]/' $(fl)`;
  c2 = `perl -pne 's/()[[:punct:]]//g'`;
  c3 = `perl -ne 'print join("\n", split(/\s+/,$_));print("\n")'`
  c4 = `grep -ve '^$'`
  read(pipeline(c1,c2,c3,`sort`,c4,`uniq -c`,`sort -rn` ))
end

text = munge("hunting-the-snark.txt");
(eltype(text),length(text))

lines = split(String(text),"\n")
n = length(lines)

s2 = [split(lines[i]) for i = n-9:n]

pwd()

# We can use the process in/out channels to capture the I/O
# (rev.pl is in the parent folder)
jabber = "jabberwocky.txt";
proc = open(`../rev.pl $jabber`,"r+");
close(proc.in);

proc.out

poem = readlines(proc.out);
close(proc.out)

poem

# We are not limited to just Perl(5)
# Get Perl6 from https://rakudo.org/files/
# Put it on the path OR setup a symbolic link to the binary: 
# viz: perl6 -> /Applications/Rakudo/bin/perl6

run(`which perl6`);

# Use Perl6 to find the line of the greatest length
#
run(`perl6 -e 'my $max=""; 
      for (lines) {$max = $_ if .chars > $max.chars};
      END { $max.say }' hunting-the-snark.txt`);

# These techniques can be used with any command processor
# This turns all characters after 40th spot to RED
#
run(`ruby -e 'w = $*.shift; $<.each { |l| puts "#{l}\e[31m#{l.chop!.slice!(w.to_i..-1)}\e[0m" }' 40 hunting-the-snark.txt`);


# And also with Python
# Encode the jabberwocky as Base 64

cd(ENV["HOME"]*"/PacktPub/Alice")
f1="jabberwocky.txt"
f2="jabberwocky.b64"
b64encode = `python -c 'import base64,sys; base64.encode(open(sys.argv[1],"rb"),open(sys.argv[2],"wb"))' $f1 $f2`
run(b64encode);

# Check it
run(`cat $f2`);


# Now reverse the process
#
run(`base64 --decode $f2`); 


cd(string(ENV["HOME"],"/PacktPub/Chp05")); # You may have to change this

# I got an error after exiting the following process
# So have wrapped it in a try/catch block
try run(`wc $(readdir())`) catch end

# Define a macro to trap any run(CMD) errors
macro traprun(c)
  quote
    if typeof($(esc(c))) == Cmd
      try
        run($(esc(c)))
     catch
     end
    end
  end
end

# Filter on a regular expression
# This will get rid of the directory warnings above
function filter(pat::Regex, dir=".")
  a = Any[]
  for f in readdir(dir)
    occursin(pat,f) && push!(a, f)
  end
  return a
end

# ... and find all the files ending in '.txt'


@traprun `wc $(filter(r"\.txt$"))`;

# Find all the Jupyter notebooks 
cd("..")
@traprun `find "/Users/malcolm/PacktPub" -name Chp\*.ipynb`;

# Location of an Apache access log

logf = "access_log";
run(`/usr/bin/wc -l $logf`);

# First 5 lines - using head command
# Can use the '$' prefix inside the backticks
run(`head -5 $logf`)

# Using PERL this is not a problem
# Do a "tail -5" as a one-liner to show last 5 lines

cmd = `perl -ne 'push @a, $_; @a = @a[@a-5..$#a]; END { print @a }' $logf`
run(cmd);

;cat hcount.pl

# Show the top 10 IP addresses in the weblog

println("Top Ten Hitters ($logf)");
@time run(`./hcount.pl $logf`);
