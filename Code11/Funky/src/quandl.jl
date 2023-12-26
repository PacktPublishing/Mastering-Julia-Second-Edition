const QURL = "https://www.quandl.com/api/v3/datasets/";

function quandl(qdset::String; apikey="", order="A", type="TB") 
  url  = string(QURL,qdset,".csv")
  (length(apikey) > 0) && (url = string(url,"?apikey=",apikey))        
  if !((order=="A") || (order=="D"))
    println("Sort order must be A or D")
    return nothing  
  end
  if !((type=="TA") || (type=="TB") || (type=="DF"))
    println("Output type must be TA or TB or DF")
    return nothing  
  end
  resp = HTTP.get(url);
  if resp.status == 200
    df = CSV.read(IOBuffer(resp.body));
  else
    println("Can't get data from Quandl for dataset: ", qdset)
  end
  if order == "A"
    sort!(df) 
  end
  if type == "TA"
    return df2ta(df)
  elseif type == "TB"
    return df2tb(df)  
  else
    return df
  end
end

function df2ta(df::DataFrame)
# Only check that df has a Date field
  nm = names(df);
  if reduce(|, map(x -> x == :Date, nm ))
    ts = Date.(df[:Date]);
    vals = float.(hcat(df[:Open],df[:High],df[:Low],df[:Close],df[:Volume]));
    cols = [:Open,:High,:Low,:Close,:Volume];
    return TimeArray(ts,vals,cols,nothing)
  end
  nothing
end

function df2tb(df::DataFrame)
  return table(df)
end

timearray(df::DataFrame) = df2ta(df)
idxtable(df::DataFrame)  = df2tb(df)
