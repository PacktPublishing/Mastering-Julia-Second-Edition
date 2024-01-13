## [Sockets] ##

# Define two sockets and bind them to ports 8011 & 8012

using Sockets
s1 = UDPSocket();
s2 = UDPSocket();
(s1, s2)

bind(s1,ip"127.0.0.1",8011)
bind(s2,ip"127.0.0.1",8012)

#  Send a message between the sockets
using Dates
send(s2,ip"127.0.0.1",8011,string(Dates.now()));

(wip, msg) = recvfrom(s1);
println(String(msg))

# Connect to website on the default port 80
sock = connect("mj2.website",80)
close(sock)

# Get website IP address
getaddrinfo("mj2.website")

# Create echo server
close(s2); close(s1);


## ["Looking-Glass World" Echo Server] ##

#! /usr/bin/env julia -q
#
using ArgParse, Sockets
const ECHO_PORT = 4000
const ECHO_HOST = "localhost"
function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--server", "-s"
    help = "hostname of the echo server"
    default = ECHO_HOST
    "--port", "-p"
    help = "port number running the service"
    arg_type = Int
    default = ECHO_PORT
    "--quiet", "-q"
    help = "run quietly, i.e. with no server output"
    action = :store_true
  end
  return parse_args(s)
end
pa = parse_command
ehost = pa["server"];
eport = pa["port"];
vflag = !pa["quiet"];
pp = (ehost == "localhost" ? "" : "$ehost>");

vflag && println("Listening on port $eport")
server = listen(eport)
while true
  conn = accept(server)
  @async begin
    try
     while true
       s0 = readline(conn)
       s1 = chomp(s0)
      (length(chomp(s1)) > 0) && (s2 = reverse(s1)))
       if s2 == "."
         println("Done.")
         close(conn)
         exit(0)
       else
         write(conn,string(pp,s2,"\r\n"))
       end
      end
     catch err
       println("Connection lost: $err")
     exit(1)
   end
  end
end


## [Quotes server] ##

using HTTP, SQLite, DataFrames
DQ = ENV["HOME"]*"/MJ2/DataSources/Quotes";
db = SQLite.DB("$DQ/quotes.db");

function getquote(db::SQLite.DB)
  sql = "select count(*) as nq from quotes";
  res = DBInterface.execute(db, sql) |> DataFrame
  nq = res[!,:nq][1];
  sql = """select a.autname, q.quotext 
       from quotes q
       join authors a on a.id = q.aid
       where q.id = (1 + abs(random() % $nq))"""
  res = DBInterface.execute(db, sql) |> DataFrame
  auth = res[!,:autname][1]
  if (auth == "missing") 
    auth = "Anon"
  end
  quotext = res[!,:quotext][1]
  (auth, quotext)
end

(auth,quot) = getquote(db)
@printf "%s => %s" auth quot

HTTP.listen("localhost",8082) do http::HTTP.Stream
  (auth, quot) = getquote(db)
  HTTP.setstatus(http, 200)
  HTTP.setheader(http, "Content-Type" => "text/html")
  HTTP.startwrite(http)
  write(http, "<html><body>\n")
  write(http, "$quot<br/>\n")
  write(http, "<i>$auth</i><br/>\n")
  write(http,"</body></html>\n")
end


## [Routing] ##

using HTTP, JSON3, StructTypes, Sockets
mutable struct Pet
    id::Int
    userId::Base.UUID
    type::String
    breed::String
    name::String
    Pet() = new()
end

StructTypes.StructType(::Type{Pet})=StructTypes.Mutable()

const PETS = Dict{Int, Pet}()

const NEXT_ID = Ref(0)
function getNextId()
    id = NEXT_ID[]
    NEXT_ID[] += 1
    return id
end

function addPet(req::HTTP.Request)
    pet = JSON3.read(req.body, Pet)
    pet.id = getNextId()
    PETS[pet.id] = pet
    return HTTP.Response(200, JSON3.write(pet))
end

function getPet(req::HTTP.Request)
    PetId = HTTP.getparams(req)["id"]
    pet= PETS[parse(Int, PetId)]
    return HTTP.Response(200, JSON3.write(pet))
end

PET_SERVER = HTTP.Router();

HTTP.register!(PET_SERVER, "POST", "/pet", addPet);
HTTP.register!(PET_SERVER, "GET", "/pet/{id}", getPet);
const PORT = 7878;
server = HTTP.serve!(PET_SERVER, Sockets.localhost, PORT);

p0 = Pet(); 
p0.type = "cat"; 
p0.breed = "Domestic"; 
p0.name = "Harry";
resp = HTTP.post("http://localhost:$PORT/pet", [], JSON3.write(p0));
resp.status

p1 = Pet();
p1.type = "dog";
p1.breed = "Jack Russell";
p1.name = "Benny";
resp = HTTP.post("http://localhost:$PORT/pet", [], JSON3.write(p1));
resp.status

p2 = Pet();
p2.type = "dog";
p2.breed = "Border Collie";
p2.name = "Jessie";

resp = HTTP.post("http://localhost:$PORT/pet", [], JSON3.write(p2));

resp.status
resp = HTTP.get("http://localhost:$PORT/pet/1")

pp = JSON3.read(resp.body,Pet)
println("$(pp.name) is a $(pp.breed) of $(pp.type)")
close(server)


## [Mux] ##

bye = """
Surely you are not leaving yet?
Don't call me Shirley!
""";

using Mux
@app muxtest = (
  Mux.defaults,
  page(respond("/","<h1>We are off to catch the White Whale</h1>")),
  page("/about",respond("Call me,<b>Ismail</b>")),
  page("/user/:user", 
    req -> "Hello, $(req[:params][:user])!"),
  page("/bye", respond(bye)),
  Mux.notfound())

serve(muxtest)

using AssetRegistry
bacon = ENV["HOME"]*"/MJ2/DataSources/Files/bacon.html";
assb = AssetRegistry.register(bacon)
"/assetserver/243f4cc0f7c3769dcb9b8da88895c744ceeb8bb5-bacon.html"


## [Web Crawler] ##

using HTTP, Gumbo

# Define a function to get a page, check its status, and that it has content.

function fetchpage(url)
  response = HTTP.get(url)
  if response.status == 200 && 
    parse(Int,Dict(response.headers)["Content-Length"]) > 0
    String(response.body)
  else
    ""
  end
end

function extractlinks(elem, links)
  if isa(elem, HTMLElement) &&
         tag(elem) == :a &&
         in("href", collect(keys(attrs(elem))))

    url = getattr(elem, "href")
    startswith(url,"/") && 
         length(url) > 1 && push!(links,url)
  end
  for child in children(elem)
    extractlinks(child, links)
  end
end

content = fetchpage("https://julialang.org")

jinx = String[]
if !isempty(content)
  dom = Gumbo.parsehtml(content)
  extractlinks(dom.root, jinx)
end

display(unique(jinx))


## [Genie] ##

using Genie
route("/looking") do
  "Here's looking at you, kid !!!"
end

up(8888)

using Genie.Renderer.Html
route("/kissing") do
  h2("Just put your lips together and blow.") |> html 
end

using Genie.Renderer.Json
route("/playing") do
   json("Play it again, Sam") |> json
end

resp = HTTP.get("http://localhost:8888/playing");
println(String(resp.body))

Genie.Generator.newapp_webservice("BlueEyes")

[Tasks and remote procedures]

ta = Task(() -> let s = 0; sum(s += i for i in 1:100); println(s) end)

istaskdone(ta)
schedule(ta)


# Define a Fibonacci function
fib(n::Integer) = n < 3 ? 1 : fib(n-1) + fib(n-2);

# ... and create an asynchronous task
chnl = Channel{Int}(1)
m = 8
@async for i=1:m
  put!(chnl, fib(i))
end

for i in 1:m
  println(i, " => ", take!(chnl))
end

close(chnl);

chnl = Channel{Int}(1)
m = 20
@async for i=1:m
  put!(chnl, fib(i))
end

for p in chnl
  if p > 100
    break
  else
    println(p, " ")
    yield()
  end
end
close(chnl)


## [Remote procedures] ##

using Distributed
addprocs(4)      # Add 4 additional processors

# Run this on PID = 2
r1 = remotecall(x -> x^2, 2, π)

# ... and get the result
fetch(r1)

remotecall_fetch(sin, 5, π/4)

using Distributed
addprocs(4)     # Only necessary if not done already

@everywhere function fac(n::Integer)
  @assert n > 0
  return ((n == 1) ? 1 : n*fac(n-1))
end

# Spawn the task over all the workers and calculate a different factorial on each worker.
r = [@spawnat w fac(3 + 2*w) for w in workers()];

#= This set up a 4-element vector, spawning the factorial on each processor for a value 
of 3 + 2 times the worker number,  that is values of: 7, 9, 11, 13 =#

# Return all the results
for id in r 
  @show id, fetch(id)
end

# Just as a check ...
@show(fac(7),fac(9),fac(11),fac(13))


## [ Needles and PI(ns) ] ##

# A sequential version
function needles_seq(n::Integer)
  @assert n > 0
  k = 0
  for i = 1:n
    ρ = rand()
    φ = (rand() * π) - π / 2 # angle at which needle falls
    xr = ρ + cos(φ)/2 # x-location of needle
    xl = ρ - cos(φ)/2
# Count times when needle crosses either x == 0 or x == 1?
    k += (xr >= 1 || xl <= 0) ? 1 : 0
  end
  m = n – k
  return (n / k * 2)
end;

# A parallel version
function needles_par(n)
  k = @distributed (+) for i = 1:n
    ρ = rand()
    φ = (rand() * π) - π / 2
    xr = ρ + cos(φ)/2
    xl = ρ - cos(φ)/2
    (xr >= 1 || xl <= 0) ? 1 : 0
  end
  m = n – k
  return (n / k * 2)
end

needles_seq(1000)
needles_par(250)

@time needles_seq(100*10^6)
@time needles_par(25*10^6)


## [Distributed arrays and Map-Reduce] ##

using Distributed
addprocs(4);

@everywhere using DistributedArrays, Statistics

da = drandn(10,10,10,10)
da[1][1][1][1]      

fieldnames(typeof(da))
da.dims
da.pids
da.indices

for ip in procs(da)
  rp = @spawnat ip mean(localpart(da))
  @show fetch(rp)
end

# Calculate the mean (~ 0.0)
using Printf
let
  μ = 0.0
  for ip in procs(da)
    rp = @spawnat ip mean(localpart(da))
    μ += fetch(rp)
  end
  @printf "\nMean = %.6f\n" μ/length(workers())
end

# ... and the standard deviation (~ 1.0)
let
  σ = 0.0
  for ip in procs(da)
    rp = @spawnat ip std(localpart(da))
    σ += fetch(rp)
  end
  @printf "Standard Deviation = %.6f\n" σ/length(workers())
end

@everywhere using Distributed, DistributedArrays, SpecialFunctions
aa = [rand() for i = 1:100,j = 1:100];

da = distribute(aa);
da.dims
da.indices

# This mapping works ...
@elapsed bb = map(x -> abs(bessely0(x)), aa)

# ... and also so does this just a little more slowly!
@elapsed db = map(x -> abs(bessely0(x)), da)

# Pick an element at random
i = rand(1:100); j = rand(1:100); (i,j)
bb[i,j] == db[i,j]

# . . . and that the overall sum of all the 10000 Bessel functions is the same
abs(sum(bb) - sum(db))


## [Distributed data sources] ##

import Pkg; Pkg.activate(".")
using Distributions, StatsBase, OnlineStats
using DataFrames, DTables, Query, CSV, Printf

DS=ENV["HOME"]*"/MJ2/DataSources/CSV";
players = ["$DS/P1.csv","$DS/P2.csv","$DS/P3.csv"]

d = DTable(CSV.File, players);
tabletype!(d)

df = fetch(d);
dp = DataFrame(df)

dt = CSV.File("$DS/teams.csv") |> DataFrame

dp1 = dp[:,[:surname,:team,:position]];
dt1 = dt[:,[:team,:games,:wins]]; 

dp2 = dp[:,[:surname,:team,:position, :saves]] |>
  @filter(_.position == "goalkeeper" && _.saves >= 16) |> 
  @orderby_descending(_.saves) |>  DataFrame;  

println(dp2[:,[:surname,:team,:saves]])

dp3 = dp |> @groupby(_.team) |>
  @map({Team=key(_), Shots=sum(_.shots)}) |> DataFrame

dj = dt |>  @join(dp3, _.team, _.Team, 
    {_.team, _.goalsFor, __.Shots} ) |> 
    @orderby(_.team) |>  DataFrame;

# Printing this out as a percentage success rate
for r in eachrow(dj)
  print(r[:team], " => ", 
        round(100*r[:goalsFor]/r[:Shots], digits=2),"%\n")
end

pc(x,y) = round(100*x/y, digits=2
count = 0
for r in eachrow(dj)
  count += 1
  (count == 1) && 
     @printf "%10s   %s\n" "Team" "Shot success"
  team = r[:team]
  if team in ["England","Germany","North Korea","USA"]
       goals = r[:goalsFor]
       shots = r[:Shots]
       @printf("%12s %6.2f %%\n", team, pc(goals,shots))
    end
  end

