## [SQLite Quotes Database ] ##

# Script to create the SQL loader

etl.jl quodata.tsv [> loader.sql

using DelimitedFiles

nargs = length(ARGS);
if nargs == 1
  tsvfile = ARGS[1];
else
  println("usage: etl.jl tsvfile");
  exit();
end

# One liner to double up single quotes
escticks(s) = replace(s, "'" => "''");

# Read all file into a matrix (using DelimitedFiles)
# The first dimension is number of lines

qq = readdlm(tsvfile, '\t')
n = size(qq)[1];

# Store all categories in a dictionary
j = 0;
cats = Dict{String,Int64}();

# Main loop to load up the quotes table
for i = 1:n
  cat = qq[i,1];
  if haskey(cats,cat)
    jd = cats[cat];
  else
    global j = j + 1;
    jd = j;
    cats[cat] = jd;
  end
  sql = "insert into quotes values($i,$jd,";
  if (length(qq[i,2]) > 0)
    sql *= string("'", escticks(qq[i,2]), "',");
  else
    sql *= string("null,");
  end
  sql *= string("'", escticks(qq[i,3]), "');");
  println(sql);
end

# Now dump the categories
for cat = keys(cats)
  jd = cats[cat];
  println("insert into categories($jd,'$cat');");
end

# Interfacing with  SQLite

$ > sqlite3

# This is relative from a shell on OSX (or Linux) and ...
# the prompt will depend on your own O/S setup
sqlite> .read build.sql
sqlite> .read load.sql
sqlite> .save quotes.db
sqlite> .ex

$> cat build.sql | sqlite3 quotes.db
$> julia etl.jl quodata.tsv | sqlite3 quotes.db

using SQLite

# SQLiteDB() etc., are not exported, so fully qualify the call
db = SQLite.DB("quotes.db")

# Check which tables are in the database
SQLite,tables(db)

SQLite.columns(db,"quotes")
sql = "select count(*) from quotes";
df = DataFrame(SQLite.Query(db, sql);

df[1]
select * from quotes limit 10

using Feather
Feather.write("QuoSnap01.feather", df)

dfx = Feather.read("QuoSnap01.feather");
size(dfx)

dfx[1,1]

sql = "select q.quoname from quotes q ";
sql *= " where q.author = 'Oscar Wilde'";

# Alternatively, if we do not wish to retain the dataframe
# ... we can just pipe it using the |> operator

SQLite.Query(db,sql) |> DataFrame
sql="select q.quoname,q.author,c.catname from quotes q ";
sql *= "join categories c on q.cid = c.id limit 5";
df = DataFrame(SQLite.Query(db,sql))

## [ SQL Databases ] ##

# MySQL
# Chinook
# ODBC

#!/bin/bash
cat /etc/odbcinst.ini

#!bash
cat odbc.ini

# Run these from a Unix shell
odbcinst -q -d; # => [MYSQL]odbcinst -q -s; # => [Chinook]
isql -v Chinook malcolm mypasswd

using ODBC

# Connect via the valid DSN
dsn = ODBC.DSN("Chinook",usr="malcolm",pwd="mypasswd");
df = ODBC.query(dsn, "select count(*) from Customers");
println("Number of customers: , df[1])Number of customers: 59

sql = "select a.LastName, a.FirstName,";
sql *= " count(b.InvoiceId) as Invs, sum(b.Total) as Amt";
sql *= " from Customer a";
sql *= " join Invoice b on a.CustomerId = b.CustomerId";
sql *= " group by a.LastName having Amt >= 45.00";
sql *= " order by Amt desc;";

df = DataFrame(SQLite.Query(db,sql));

# Looping through the dataframe, to produce better output
using Printf
for i in 1:size(df)[1]
  LastName = df[:LastName][i]
  FirstName = df[:FirstName][i]
  Amt = df[:Amt][i]
  @printf "%10s %10s %10.2f\n" LastName FirstName Amt
end

sql = "select a.LastName, a.FirstName, d.Name as TrackName";
sql *= " from Customer a";
sql *= " join Invoice b on a.CustomerId = b.CustomerId";
sql *= " join InvoiceLine c on b.InvoiceId = c.InvoiceId";
sql *= " join Track d on c.TrackId = d.TrackId";
sql *= " where a.LastName = 'Cunningham' and a.FirstName =

sql *= " limit 5;";
df = DataFrame(SQLite.Query(db,sql));
for i in 1:size(df)[1]
  TrackName = df[:TrackName][i]
  @printf "%s\n" TrackName
end


# (Native) MySQL

df = DataFrame(MySQL.query(conn,"show tables"));
tb = df[:,1]
print("Chinook tables:: ")
for i = 1:size(tb)[1]
  print(tb[i]," ")
end

sql = "select FirstName,LastName,Address,City,State"
sql *= " from Customer where LastName like 'C%'";
DataFrame(MySQL.query(conn, sql))

sql = "select a.FirstName,a.LastName, sum(b.Total) as 'Total spent'";
sql *= " from Customer a";
sql *= " join Invoice b on a.CustomerId = b.CustomerId";
sql *= " group by a.FirstName, a.LastName"
sql *= " having a.LastName like 'C%'"
MySQL.query(conn,sql) |> DataFrame
MySQL.disconnect(conn);

# Using PyCall

# Run is in Python
# The mysql connector has to be installed
#
import mysql.connector as mc
cnx = mc.connect(user="malcolm", password="mypasswd")
csr = cnx.cursor()
qry = "SELECT * FROM Chinook.Genre"
csr.execute(qry)
for vals in csr:
   print(vals)
csr.close()
cnx.close()

using PyCall
@pyimport mysql.connector as mc;
cnx = mc.connect(user="malcolm", password="mypasswd");

query = "SELECT * FROM Chinook.Genre"
csr[:execute](query)
for vals in csr
  id = vals[1]
  genre = vals[2]
  @printf "ID: %2d, %s\n" id genre
end

csr[:close]()
cnx[:close]()
derbytools.ja
# Java and JDBC
# Derby

#! /bin/bash
DH = $HOME/MJ2/db-derby # Save some typing
export JAVA_HOME=$(/usr/libexec/java_home)
export PATH=$PATH:$DH/bin
export CLASSPATH=$DH/lib/derbytools.jar:$DH/lib/derbynet.jar:.

# Check that the setup is OK by using sysinfo sysinfo
# Change to the Quotes directory to pickup the load files

cd $HOME/MJ2/db-derby
startNetworkServer &

cd ../DataSources/Quotes

# Now use the ij utility to look the dataset and ...
# ... run a query to see if the build/load was OK.
ij> connect 'jdbc:derby:Quotes;create=true';
ij> run 'build.sql';
ij> run 'qloader.sql';
ij> select count(*) from quotes;

using DataStructures
s = Stack{String}();
push!(s, @__DIR__);
cd(ENV["HOME"]*"/MJ2/DataSources/Quotes")

ENV["JULIA_COPY_STACKS"]=1;

using JavaCall

CP = "/Users/malcolm/MJ2/db-derby/lib/derbytools.jar";

JavaCall.addClassPath(CP);
JavaCall.addOpts("-Xmx1024M"); 
JavaCall.addOpts("-Xrs");
JavaCall.init()

jsd = @jimport java.sql.DriverManager;
dbURL = "jdbc:derby:Quotes";
conn = nothing;

jsp = @jimport java.sql.PreparedStatement;
jsr = @jimport java.sql.ResultSet;
try
  conn = jcall(jsd,"getConnection",JavaObject,(JString,),dbURL);
  sql  = "select count(*) as K from quotes";
  stmt = jcall(conn,"prepareStatement",jsp,(JString,),sql);
  res  = jcall(stmt,"executeQuery",jsr,());
  k    = jcall(res,"getString",JString,(JString,),"K");
catch e
  println(e);
finally
  if conn != nothing
    jcall(conn,"close",Void,());
  end
end;

println("\nNumber of quotes in database: $k);

# Create a cursor to the Quotes database
# ... and execute a SQL statement
csr = cursor("jdbc:derby:Quotes")
execute!(csr, "select * from my_table;")

# Then iterate over rows
for row in rows(csr)
# access the data in 'row'
end
close(csr);

# An alternative is to return a dataframe directly

# ... then then work with it.
df = JDBC.load(DataFrame, cursor(csr), "select * from quotes")

# This works not only for DataFrames but for any Data.Sink

# PostgreSQL

# !/bin/bash
export PATH = /Applications/Postgres.app/Contents/Versions/latest/bin:$PATH
pg_ctl -D /data/psq -l logfile start
ps -ef | grep postgres # use this to check it has started

#! /bin/bash
# The directory below is where my Chinook files are kept
cd /Users/malcolm/PacktPub/Chp09/Chinook/Dataset
createdb chinook
psql -d chinook -a -f Chinook_PostgreSql.sql

using LibPQ, DataStreams
conn = LibPQ.Connection("dbname=chinook");

res = execute(conn, "SELECT * FROM \"MediaType\"");
res = Data.stream!(res, NamedTuple)

# Alternatively using a single routine : fetch!()
res = fetch!(NamedTuple, execute(conn, "SELECT * as FROM \"MediaType\""));
res[:Name]

sqlx = """
select e."FirstName", e."LastName", count(i."InvoiceId") as "Sales"
from "Employee" as e
join "Customer" as c on e."EmployeeId" = c."SupportRepId"
join "Invoice" as i on i."CustomerId" = c."CustomerId"
group by e."EmployeeId"
""";

# Loop through the NameTuple
qry = fetch!(NamedTuple, execute(conn, sqlx))
using Printf
for i in 1:length(qry)
@printf("%s %s has %d sales.\n",
qry.FirstName[i], qry.LastName[i]. qry.Sales[i])
end

# NoSQL databases
# Key-Value Datastores

# Redis

$ redis-cli -p 99999 \
-h redis-99999.z99.eu-west-1-2.ec2.cloud.redislabs.com

redis> PING

# Set a simple key-value, check it is set and the retrieve it
redis> set me "Malcolm Sherrington"
redis> keys *
redis> get me

const RLHOST =
"redis-99999.z99.eu-west-1-2.ec2.cloud.redislabs.com";

const RLPORT = 99999
const RLAUTH = "aBcde18FGhiJklM77noPQrstuv"
conn = RedisConnection(host=RLHOST, port=RLPORT, password=RLAUTH)

# Check with connection is alive and setup a key-values
ping(conn)

set(conn,"me","Malcolm Sherrington")

# Note that the return status for the SET is the string "OK"
# List the keys
keys(conn,*)

# ... and return the value
get(conn,"me")

using HTTP,CSV
const QURL = "https://www.quandl.com/api/v3/datasets/";
const QAPIKEY = "ABCd1357EfgH2468yZ";

# Return the dataset as a dataframe, in ascending date order
function quandl(qdset::String, apikey::String="")
  url = string(QURL,qdset)
  (length(apikey) > 0) && (url = string(url,"?apikey=",apikey))
  resp = HTTP.get(url);
  if resp.status == 200
    df = CSV.read(IOBuffer(resp.body));
    return sort!(df)
  else
    error("Can't get data from Quandl: $qdset")
  end
end

qdf1 = quandl("WIKI/AAPL.csv", QAPIKEY);
df = float.(qdf[:Close]);
sf = 1.0/df[1];
aapl = sf .* df;
qdf2 = quandl("WIKI/MSFT.csv", QAPIKEY);
df = float.(qdf[:Close]);
sf = 1.0/df[1];
julia > msft = sf .* df ;
n = minimum([length(aapl),length(msft)]);
t = collect(1:n);.

conn = RedisConnection() julia> for i = n-1:-1:0
rpush(conn,'APPL~Close',aapl[end-i])
rpush(conn,'MSFT~Close',msft[end-i])
end

# Redis lists are zero-based
aapl_rdata = float.(lrange(conn,"AAPL~Close",0,n-1);
msft_rdata = float.(lrange(conn,"MSFT~Close",0,n-1);

plot(t, aapl_rdata, Color="red", label="AAPL")
plot(t, msft_rdata, Color="blue", label="MSFT")
xlabel("Time")
ylabel("Scaled Price")
legend()

# Other Key-Value systems
# Memcache
# LMDB
# LevelDB

# Document databases

# MongoDB

$> mongo
$> use test
$> db.getCollectionNames()

$> cd /Users/malcolm/PacktPub/Chp09/Quotes
$> mongo < quotes.js

$> mongo
> qc = db.quotes
> qc.count()

> qc.find({author:"Oscar Wilde"})

(v1.1) pkg> add LibBSON#master
(v1.1) pkg> add https://github.com/ScottPJones/Mongo.jl#master

import Mongoc
const mc = Mongoc
client = mc.Client(); # Create a connection to Mongo
db = client["test"];
quotes = db["quotes"];

# Check that we are connected and can see the quotes collection
length(quotes)
```
# Note the need for triple quoting the BSON string
# Escaping in single "s (i.e. \") which work equally well

mc.find_one(quotes, mc.BSON("""{"author" : "Oscar Wilde"}"""))

# CRUD

doc1 = mc.BSON()
doc1["author"] = "Orson Welles";
doc1["quote"] = "If you want a happy ending, it all depends on where you stop your story";
doc1["category"] = "Words of Wisdom";

doc2 = mc.BSON()
doc2["author"] = "Bo Bennett";
doc2["quote"] = "Visualization is daydreaming with a purpose";
doc2["category" ] = "Computing";

mc.as_dict(doc2)

res1 = push!(quotes, doc1)
oid1 = res1.inserted_oid
length(quotes)

selr = mc.BSON();
selr["_id"] = oid1;
mc.delete_one(quotes, selr)
length(quotes)
```
# RESTful interfacing

rts = chomp(read(`curl –s http://amisllp.com/now.php`, String))

println(rts);

# JSON and BSON

company = """{
"founder" : "Malcolm Sherrington",
"company" : "AMIS Consulting LLP",
"website" : "amisllp.com",
"partners" : 5,
"incorporated" : "1981-06-01",
"areas_of_work" :
["Aerospace","Healthcare","Finance","Web Development"],
"experience_years : [31, 26, 15, 21],
"programming_languages" :
["Fortran","C","Ada","Perl","PHP","Julia","Python","R"]
}"""

import JSON
amis = JSON.parse(company)

using Dates
Date(amis["incorporated"])
Integer.(amis["experience_years"])

# BSON
# Web Databases (CouchDB)

# HTTP package

testpage = "http://amisllp.com/testpage.html";
r = HTTP.get(url, retry=false, readtimeout=30);
r.status
r.headers

String(r.body)

HTTP.open("GET", testpage) do http
i = 0
while (!eof(http) && (i <= 8))
i = i + 1
s = read(http)
(i > 1) && println(s)
end
end

# CouchDB

$> curl http://localhost:5984

$> curl -X PUT http://localhost:5984/quotes
$> curl http://localhost:5984/_all_dbs

$> curl -H 'Content-Type: application/json' \

using HTTP, JSON
cdp ="http://localhost:5984"
dbs = String(HTTP.get(cdp*("/_all_dbs")).body)
JSON.parse(dbs)

json = JSON.get("http://127.0.0.1:5984/quotes/_all_docs");

rec = JSON.parse(json)
rec["rows"][1]["key"]

db = "quotes";
key = "ecc520bf0083e5a48907e52f1f0013bb";
json = String(HTTP.get("$cdp/$db/$key").body);
doc = JSON.parse(json)

using Printf
@printf "%s [%s]" doc["quote"] doc["author"]

$> npm install -g couchimport

$> cat quotes.tsv | \
couchimport --url http://localhost:5984 --db quotes

doc = JSON.parse(String(HTTP.get(cdp*("//quotes")).body))
doc["doc_count"]


# JuliaDB (defunct)

julia > using JuliaDB, IndexedTables
julia > path = joinpath(homedir(),

# Indicate that date and ticker fields should be indexed.
stockdata = loadndsparse(path,

# Single values may be shown in the usual way
using Dates
stockdata[Date("2010-06-01"), "GOOGL"]

# Or we can define a date range and select a couple of stocks
stockdata[Date("2012-01"):Dates.Month(1):Date("2014-12"), ["GOOGL", "KO"]]

# Values of Goldman Sachs with closing prices in [100.0,140.0]
filter(x -> x.close >= 100.0 && x.close <= 140.0, stockdata[:, "GS"])

# Get records for Xerox where the first tracing day (of the month) is a Friday

filter((1=>Dates.isfriday, googl = stockdata[:, ["GOOGL"]];
spread = map(x -> x.high - x.low, googl)
round(reduce(+,(mean.(spread)))/length(spread), digits=4)
gain = map(x -> x.open - x.close, googl)
round(reduce(+,(mean.(gain)))/length(gain), digits=4)

using OnlineStats
groupreduce(Mean(),stockdata,:ticker; select=:close)

using StatsPlots
@df stockdata plot(:date, :close, group=:ticker, layout = 4, legend = :topleft)
