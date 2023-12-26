module Funky

using HTTP,CSV,DataFrames,TimeSeries,IndexedTables
using Printf, Dates, Statistics

include("ftop.jl")
include("queens.jl")
include("panto.jl"
include("quandl.jl")

const PUNCTS =  [' ','\n','\t','-','.',',',':',';','!','?','\'','"'];

function basel(N::Integer)
  @assert N > 0
  s = 0.0
  for i = 1:N
    s += 1.0/float(i)^2
  end
  return s
end

fac(N::Integer) =
  (N < 1) ? throw(ArgumentError("N must be positive")) : reduce(*,collect(big.(1:N)))

function fib(n::Integer)
  @assert n > 0
  (a, b) = (big(0), big(1))
  while n > 0
    (a, b) = (b, a+b)
    n -= 1
  end
  return a
end

function hailstone(n::Integer)
   k = 1
   a = [n]
   while n > 1
      n = (n % 2 == 0) ? n >> 1 : 3n + 1
      push!(a,n)
      k += 1
   end
   a
end

function kempner(n::Integer)
  @assert n > 0
  s = 0.0
  r9  = r"9"
  r9x = r"9{2}"
  for i in 1:n
    if (match(r9,string(i)) == nothing) || (match(r9x,string(i)) != nothing)
      s += 1.0/float(i)
    end
  end
  return s
end

function isdate(s::String)
  b = false
  try    
    Date(s)
    b = true 
  finally
    return b
  end    
end

function wordcount(text)
  wds = split(lowercase(text), PUNCTS; keepempty = false)
  d = Dict()
  for w = wds
     d[w] = get(d,w,0)+1
  end
  return d
end

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

function filter(pat::Regex, dir=".")
  a = Any[]
  for f in readdir(dir)
    occursin(pat,f) && push!(a, f)
  end
  return a
end

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

macro until(condition, block)
    quote
        while true
            $(esc(block))
            if $(esc(condition))
                break
            end
        end
    end
end

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
end

export basel, fac, fib, hailstone, kempner,  isdate, ftop, quandl
export @traprun, @bmk, @until 

end
