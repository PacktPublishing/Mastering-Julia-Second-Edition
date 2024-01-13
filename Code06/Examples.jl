# First either activate project or run setup.jl

import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames, RDatasets, CSV, JLD, HDF5, Query, Match
using Plots, PyPlot, GR, StatsBase, Distributions
using GLM, HypothesisTests, KernelDensity, Clustering

using Statistics, Printf, Dates, DelimitedFiles

 write(stdout, "Help ME!!!\n")

write(stderr,"Error messages should come here.")

#= 
read(stdin) will not work in a Jupyter notebook

It is possible to use a terminal

Use File/Open and choose New->Terminal from the pulldowmn menu
Run the Julia REPL

TIP: 
On my system (OSX) I have setup the following symbolic link in my 
~/bin folder (which is on the executable path)

julia -> /Applications/Julia-1.0.app/Contents/Resources/julia/bin/julia

Try some of the following:- (CR => carriage return, ^D => control-D)

read(stdin,Char)        - type a letter + CR
read(stdin,String)      - type a string + ^D (CR is part of the string)

read(stdin,Int32)       - type a number + CR
read(stdin,UInt32)      - type same number + CR
readline(stdin)         -   "    "    "       "


# However readline will work in Julia
a = readline()

# Because readline() works in Jupyter this is OK.

function getInt()
    s = chomp(readline())
    try
        parse(eltype(1),s)
    catch ex 
        println(ex, "\nCan't convert $s to an integer")       
    end
end

getInt()

# As is this

function getFloat()
    s = chomp(readline())
    try
        parse(eltype(1.1),s)
    catch ex
        println(ex, "\nCan't convert $s to a float")
    end
end

getFloat()


bb = Array{UInt8}(undef,0)     # This will not work in a notebook

#=

Opening a file has to main syntax:

open(filename::AbstractString; keywords...)
open(filename::AbstractString, [mode::AbstractString])

If neither keywords of a mode string are given then the default is to open the file for a read

=#

#-------------------------------------------------------------------------

# Discuss diffence between s1*s2 vs String(s1,s2)

flip = ENV["HOME"]*"/MJ2/Alice/jabberwocky.txt";
isfile(flip)


# Classic paradigm
# This prints the first 4 lines of the Jabberwocky 
# We need the let/end because of the new scoping rules
let
    fin = open(flip)
    k = 0;
    while !eof(fin)
        ln = readline(fin)
        k = k+1
        (k > 4) ? break : println(ln)
    end
    close(fin)
end


# Because file processing is a common process Julia has an alternative syntax
# Note jabber4.txt is just the first verse of the Jabberwocky, save us having to count 
# and break to reduce the output

flp4 = ENV["HOME"]*"/MJ2/Alice/jabber4.txt";

open(flp4) do pn4
  while !eof(pn4)
     println(readline(pn4))
  end
end


# The file handle is not defined after the loop, so does not need to be closed.
julia> eof(pn4)

#= 
Conveniently we can read all a file in to a variable without explicitly opening it
Note that at the data is returned as a byte array (the file may be a binary)
=#
# For a text file we need to convert it

jb4 = String(read(flp4))

# NB: This needs to be split of "\n"s to reproduce the familiar first verse

# There is an additional format:  open(function, filename)

capitalize(f::IOStream) = chomp(uppercase(String(read(f))));
split(open(capitalize, flp4),"\n")

# We need the chomp() to remove the "\n" in the last line, 
# otherwise a 5-element array would be returned with the last element "".


# In the previous chapter we used Perl to reverse a lines of a poem.
# Here is a native Julia version

open( flp4) do pn4
    while !eof(pn4)
        println(reverse(chomp(readline(pn4))))
    end
end

# Again we need to chomp the line, reverse it and use println
# Otherwise the "\n" would come at the front

# We can also emulate the 'wc' or or wc() approaches from previously
# In this routine only the words or counted
# Actually this is the most difficult, counting lines and characters are easy

# We will split on white space and also pinctation marks

const PUNCTS =  [' ','\n','\t','-','.',',',':',';','!','?','\'','"'];

# The routine returns a hash (Dict) with the words found as keys and the count as the values,

function wordcount(text)
  wds = split(lowercase(text), PUNCTS; keepempty = false)
  d = Dict()
  for w = wds
     d[w] = get(d,w,0)+1
  end
  return d
end


# Test it on our 4 line Jabberwocky
wordcount(String(read(flp4)))

# Most words only occur once but 'and' & 'the' have a higher frequency

# Wordcount returns a dictionary of the words in a file
#
# So if we collect the values (i.e the counts) and sum them ...
# ... this gives a total for the file.

using Printf
AliceDir = ENV["HOME"]*"/MJ2/Alice/";  # Note the convenient trailing '/'

# Filter to look just a the '.txt' files in the Alice directory
# We can collect all the value sof the Dict in an array and sum it for the total i the file
#=
for fname in filter!(r"\.txt$", readdir(AliceDir))
  open(AliceDir*fname) do f
    n = sum(collect(values(wordcount(String(read(f))))))
    @printf "%s: %d\n" fname n
  end
end
=#

for fname in readdir(AliceDir)
  if match(r"\.txt$", fname) != nothing
    open(AliceDir*fname) do f
      n = sum(collect(values(wordcount(String(read(f))))))
      @printf "%s: %d\n" fname n
    end
  end
end


# Look at some of the characters in the Hunting of the Snark
# All of the ones in the boat begin with a 'B'

snarkDict = wordcount(String(read(AliceDir*"hunting-the-snark.txt")))

wds = ["baker","banker","barrister","beaver","bellman","boots","butcher"];
for w in wds
  @printf "%12s  => %4d\n" w snarkDict[w]
end

# Set the correct directory and open the Juliaset
# Read the first line for the "magic" number

cd(ENV["HOME"]*"/MJ2/Chp06"); 
img = open("Files/juliaset.pgm");
magic = chomp(readline(img))

#=
Check is is a PGM file
If so pick up the next line to get the image size
Open an 'output' file
This will be another PGM file of the same size ...
... so we can write 'magic' and 'params' to the output
=#
if magic == "P5"
  out = open("Files/jsetinvert.pgm", "w");
  println(out, magic);
  params = chomp(readline(img));  # => "800 400 255"
  println(out, params);
# Params splits into strings, we need integers
  (wd, ht, pmax) = parse.(Int64,split(params));
# Create a byte array and read ALL the image data in one call.
  np = wd*ht;
  buf = Array{UInt8,1}(undef,np)
  readbytes!(img, buf, np);
# Invert the gra
  bufX = [UInt8(255 - buf[i]) for i = 1:np]
  write(out,bufX))
  close(out);
else
  error("Not a NetPBM grayscale file")
end

close(img)

# Define a pseudocolor filter
#
function pseudocolor(pix)
  if pix < 64
    pr = UInt8(0); pg = UInt8(0); pb = UInt8(4*pix)
  elseif pix < 128
    pr = UInt8(0); pg = UInt8(min(4*(pix - 64),255)); pb = UInt8(255)
  elseif pix < 192
    pr = UInt8(0); pg = UInt8(255); pb = UInt8(min(4*(192 - pix),255))
  else
    pr = UInt8(min(4*(pix - 192),255))
    pg = UInt8(min(4*(256 - pix),255))
    pb = UInt8(0)
  end
  return (pr, pg, pb) 
end

pseudocolor(0xa1)

cd(ENV["HOME"]*"/MJ2/Chp06"); 
img = open("Files/juliaset.pgm");
magic = chomp(readline(img))

#=
Same logical as before
Except:-
1. Changing magic number from P5 => P6
2. Applying the pseudocolor filter
3. Writing 3 bytes for each (r,g,b)
=#

if magic == "P5"
  out = open("Files/jsetcolor.ppm", "w");
  println(out, "P6");
  params = chomp(readline(img));  # => "800 400 255"
  println(out, params);
  (wd, ht, pmax) = parse.(Int64,split(params));
  np = wd*ht;
  buf = Array{UInt8,1}(undef,np)
  readbytes!(img, buf, np);
  for j = 1:np
    (r,g,b) = pseudocolor(buf[j]);
    write(out,UInt8(r)); write(out,UInt8(g)); write(out,UInt8(b))
  end
  close(out);
else
  error("Not a NetPBM grayscale file")
end

close(out)
close(img)

#-------------------------------------------------------------------------

# Start by working with some Apple stock prices
using CSV, Statistics, Printf

cd(ENV["HOME"]*"/MJ2/Chp06/")
aaplcsv = "Files/aapl.csv"; isfile(aaplcsv)

# The CSV.File routine returns the schema
# Notice the unions as the data may have missing values
aapl = CSV.File(aaplcsv)

# Look at the fields comprising the datastructure.
fieldnames(typeof(aapl))

# We are at the first position
aapl.currentrow

# Print the first 5 rows
k = 0;
for r in aapl
    @printf "%s : %.2f\n" r.Date (r.Close - r.Open)
    global k = k + 1
    if k > 5 break end
end

# Check that we are now at the 6th row.
aapl.currentrow

# It is possible to 'pipe' the CSV.File structure to a Dataframe
# Useful to sort the dataframe in place
#
using DataFrames
df = aapl |> DataFrame

sort!(df)

# Dataframes can be queried with the "Queryverse"
# Of which more later
using Query, Dates

x = @from i in df begin
    @where i.Date >= Date(2013,12,20)
    @select {i.Date, i.Open, i.High, i.Low, i.Close}
    @collect DataFrame
end

#-----------------------------------------------------------------------

using DelimitedFiles

cd(ENV["HOME"]*"/MJ2/Chp06/")
ukhptsv = "Files/UKH-Prcs.tsv"; isfile(ukhptsv)

(ukhpData,ukhpHead) = readdlm(ukhptsv, '\t'; header=true)
ukhpHead

(ukd1, ukd2) = size(ukhpData)

typeof(ukhpData)

using PyPlot

t = collect(1:ukd1);

for i in 1:10
  plot(t,ukhpData[:,i])
end

aaplcsv = "Files/aapl.csv"; isfile(aaplcsv)

(aaplData,aaplHead) = readdlm(aaplcsv, ','; header=true)
aaplHead

aaplData[1:10, 1:6]

using Dates, PyPlot

d0 = Date("2000-01-01")
aaplDate  = reverse(Date.(aaplData[:,1]))
aaplOpen = reverse(Float64.(aaplData[:,2]))
aaplClose = reverse(Float64.(aaplData[:,5]))

const NAAPL = length(aaplDate)
aapl_days = zeros(Int64, NAAPL)

d0 = aaplDate[1]
dt = [(aaplDate[i] - d0).value for i = 1:NAAPL]

plot(dt, aaplClose)
title("Apple Stock - Closing Price")


aaplDifs = [(aaplClose[i] - aaplOpen[i]) for i = 1:NAAPL]
plot(dt,aaplDifs)
title("Apple Stock - Daily Change")

using JLD, HDF5

cd(ENV["HOME"]*"/MJ2/Chp06/")
h5file = "Files/mydata.h5"
aa = [u + v*rand() for u = 0.5:0.5:10.0, v = 0.5:0.5:6.0]

h5open(h5file, "w") do f
    write(f, "aa", aa) 
end

# Alternatively, we can say "@write h5file aa" or else h5write(h5file,"aa",aa)

# This can be read back without the h5open statement (similar to h5write)
# Also we can create a slice of the file on-the-fly

bb = h5read(h5file, "aa", (2:3:14, 4:3:10))

#-----------------------------------------------------------------------

aapljld = "Files/aapl.jld";
rm(aapljld, force=true)

jldopen(aapljld, "w") do fid
  write(fid,"aaplDate",aaplDate)
  write(fid,"aaplClose",aaplClose)
  write(fid,"aaplDifs",aaplDifs)
end

isfile(aapljld)

fid = jldopen(aapljld, "r")

aaDate = read(fid, "aaplDate")
aaDifs = read(fid, "aaplDifs")
close(fid)

fid = jldopen("Files/myaapl.jld", "w")
g = g_create(fid, "aapl")
g["aaplDate"]  = aaplDate
g["aaplClose"] = aaplClose
g["aaplDifs"]  = aaplDifs

dump(fid)

g = fid["aapl"]
adf = g["aaplDifs"]
close(fid)

# Alternatively: adf = fid["aapl/aaplDifs"]

using LightXML
xdoc = parse_file("Files/books.xml");

xtop = LightXML.root(xdoc);
println(LightXML.name(xtop));

xdoc

using Printf

for c in child_nodes(xtop)
    if is_elementnode(c)
        e = XMLElement(c)
        t = find_element(e, "title")
        title = content(t)
        genre = attribute(t, "genre")
        @printf "%30s -:- %s\n" title genre
    end
end

using Dates
for c in child_nodes(xtop)
   if is_elementnode(c)
      e = XMLElement(c)
      t = find_element(e, "title")
      genre = attribute(t, "genre")
      if genre == "Computing"
         a = find_element(e,"author")
         p = find_element(e,"price")
         curr = attribute(p, "currency")
         d = find_element(e,"publish_date")
         dc = DateTime(content(d))
         ds = string(day(dc)," ",monthname(dc)," ",year(dc))
         desc = find_element(e,"description")
         println("Title:      ", content(t))
         println("Author:     " ,content(a))    
         println("Date:       " ,ds)
         println("Price:      " ,p ," (", curr, ")")
         println(content(desc),"\n");
      end
   end
end

#-------------------------------------------------------------------------

# https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/00Index.html
# http://juliastats.github.io/StatsBase.jl/latest

using DataFrames, RDatasets

RDatasets.datasets("MASS")

using Statistics, StatsBase

quakes = dataset("datasets", "quakes");

quakes[1:5,:]

mags  = Float64.(quakes[:Mag]);
depth = Float64.(quakes[:Depth]);

cor(mags,depth)

using Plots
gr

scatter(depth,mags)

describe(mags)

mags = Float64.(quakes[:Mag]);
(m1, m2, m3, m4) = map(x -> round(x,digits=4), 
              [mean(mags), std(mags), skewness(mags), kurtosis(mags)]);

using Printf
@printf "Mean: %.4f\nStdV: %.4f\nSkew: %.4f\nKurt: %.4f\n" m1 m2 m3 m4

histmag = fit(Histogram, mags, 4.0:0.1:6.6, closed=:left)

plot(histmag)

using RDatasets

#--------------------------------------------------------------------------------------

mlmf = dataset("mlmRev","Gcsemv");
size(mlmf)

mlmf[1:5,:]

describe(mlmf)

writtenF = collect(skipmissing(mlmf[mlmf.Gender .== "F", :Written]));
writtenM = collect(skipmissing(mlmf[mlmf.Gender .== "M", :Written]));

(μWM, μWF) = round.((mean(writtenM), mean(writtenF)), digits=3)
(σWM, σWF) = round.((std(writtenM), std(writtenF)),digits=3)
(nWM, nWF) = (length(writtenM), length(writtenF))

σW = sqrt((σWM*σWM)/(nWM - 1) + (σWF*σWF)/(nWF - 1))
tt = round(abs(σWM - σWF)/σW , digits=4)    

# p ~ 0.33; 95% ~ 0.06, 90% ~ 0.13
0.9726


courseF = collect(skipmissing(mlmf[mlmf.Gender .== "F", :Course]));
courseM = collect(skipmissing(mlmf[mlmf.Gender .== "M", :Course]));


(μCM, μCF) = round.((mean(courseM), mean(courseF)), digits=3)
(σCM, σCF) = round.((std(courseM), std(courseF)),digits=3)
(nCM, nCF) = (length(writtenM), length(writtenF))

σC = sqrt((σCM*σCM)/(nCM - 1) + (σCF*σCF)/(nCF - 1))
tt = round(abs(σCM - σCF)/σC , digits=4)    

# p ~ 0.33; 95% ~ 0.06, 90% ~ 0.13 


#-------------------------------------------------------------------------

using RDatasets, KernelDensity

mlmf = dataset("mlmRev", "Gcsemv");

df = mlmf[completecases(mlmf[[:Written, :Course]]), :]

macro F64(sym)
    quote
      Float64.(skipmissing(Array($sym)))
    end
end

dc = @F64 df[:Course];
dW = @F64 df[:Written];

kdc = kde(dc)

summarystats(dc)

kdw = kde(dw)

summarystats(dw)

using PyPlot

PyPlot.plot(kdc.x, kdc.density)
PyPlot.plot(kdw.x, kdw.density, linestyle="--")

for subdf in groupby(df, :School)
  (size(subdf)[1] > 40) &&  
    let
     sch = subdf[:School][1]
     msw = mean(subdf[:Written]) 
     msc = mean(subdf[:Course])
     @printf "%10s : %8.4f %8.4f\n" sch msw msc
    end
end

#-------------------------------------------------------------------------

using HypothesisTests

df68107 = mlmf[mlmf[:School] .== "68107", :];
df68107cc = df68107[completecases(df68107[[:Written, :Course]]), :];

df68411 = mlmf[mlmf[:School] .== "68411", :];
df68411cc = df68411[completecases(df68411[[:Written, :Course]]), :];


df68107wri = @F64 df68107cc[:Written];
df68411wri = @F64 df68411cc[:Written];

UnequalVarianceTTest(df68107wri, df68411wri)


df68107cou = @F64 df68107cc[:Course];
df68411cou = @F64 df68411cc[:Course];

UnequalVarianceTTest(df68107cou, df68411cou)

using GLM

dw68411ss =  sort(sample(df68411wri,50));
dw68107ss =  sort(sample(df68107wri,50));
dc68411ss =  sort(sample(df68411cou,50));
dc68107ss =  sort(sample(df68107cou,50));

cor(dw68107ss, dw68411ss)

cor(dc68107ss, dc68411ss)

dwf = convert(DataFrame, hcat(dw68107ss, dw68411ss))
names!(dwf, [:s68107, :s68411])

dwf[1:5,:]

dcf = convert(DataFrame, hcat(dc68107ss, dc68411ss))
names!(dcf, [:s68107, :s68411])

dcf[1:5, :]

lm1 = fit(LinearModel, @formula(s68107 ~ s68411), dwf)
lm2 = fit(LinearModel, @formula(s68107 ~ s68411), dcf)

using TimeSeries, MarketData

# Two years of AAPL data (for 2000 / 2001)

ohlcv[1:5]

#=
struct TimeArray{T,N,D<:TimeType,A<:AbstractArray{T,N}} <: AbstractTimeSeries{T,N,D}
    timestamp::Vector{D}
    values::A   # an AbstractArray{T,N}
    colnames::Vector{Symbol}
    meta::ANY
end

=#

meta(ohlc)    # ohlc does not include the volume (v)
colnames(ohlc)
timestamp(ohlc)[end-2:end]
values(ohlc)[end-2:end]
cl[1:3]   # Also op, hi, lo

using PyPlot

PyPlot.plot(timestamp(cl),values(cl))

# Marketdata used to provide routines for getting data from Yahoo (no-defunct)
# Also from FRED (Federal Reserve Bank of St Louis), for US economic data
# see: https://fred.stlouisfed.org
#
# As with the ohlc dataset fred() returns a default series, the US CPI
# 
cpi =  fred()
vcat(cpi[1:3],cpi[end-2:end])

dsg10 = fred("DGS10")
vcat(dsg10[1:3],dsg10[end-2:end])

t = collect(1:length(dsg10))
y = values(dsg10)
PyPlot.plot(t,y)

# Temporal.jl is an alternative timeseries package, maintained by JuliaQuant
# and integrating with some packages from the group, such as Quandl and Indicators
# which we will meet later in the book
#
# * Motivated by the `xts` package in R and the `pandas` package in Python *
#
using Temporal, MarketData

md = MarketData   # Ensure no name clash on 'ohlc'
aaplts = TS(values(md.ohlc),timestamp(md.ohlc),colnames(md.ohlc))

using Printf
aapl_spread = aaplts[:High] - aaplts[:Low]
@printf "Average spread : %.3f\n" round(sum(aapl_spread)/length(aapl_spread), digits=3)

describe(aapl_spread.values[:,1])

plot(aapl_spread.index,aapl_spread.values[:,1])


