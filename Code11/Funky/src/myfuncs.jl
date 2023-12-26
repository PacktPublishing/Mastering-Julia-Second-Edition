function basel(N::Integer)
  @assert N > 0
  s = 0.0
  for i = 1:N
    s += 1.0/float(i)^2
  end
  return s
end

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
