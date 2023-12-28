#! /usr/bin/env julia
#
# Check on the number of arguments, print usage unless 1 or 2.
#
using DelimitedFiles

nargs = length(ARGS);
if nargs == 1 || nargs > 2 
    tsvfile = ARGS[1];
else
    println("usage: etl.jl tsvfile");
    exit(); 
end

# One liner to double up single quotes
escticks(s) = replace(s, "'" => "''");

# Read all file into a matrix, first dimensions is number of lines
qq = readdlm(tsvfile, '\t');
n = size(qq)[1];

# Going to store all categories in a dictionary
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

# Now dump the categories:

for cat = keys(cats)
  jd = cats[cat];
    println("insert into categories values($jd,'$cat');");  
end

