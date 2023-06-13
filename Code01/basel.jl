using Printf

function basel(M::Integer)
  @assert M > 0
  s = 0.0
  for i = 1:M
    s += 1.0/float(i)^2
  end
  return s
end

N = 10^9
BS = basel(N)
@Printf.printf "Estimate of Basel function for %d iterations is %.6f\n" N BS

##############################################################################

@time sum(collect(1:N)) do x
  1/(x*x)
end

0.663908 seconds (29.54 k allocations: 764.890 MiB, 9.09% gc time, 5.05% compilation time)
1.6449340568482282