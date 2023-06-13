using Printf

mutable struct Params
  N::Int
  K::Int
end

P = Params(10000, 0) 

sumsq(x,y) = x*x + y*y;

for i in 1:P.N
 if sumsq(rand(), rand()) < 1.0
    P.K += 1
 end
end

@Printf.printf "Estimate of PI for %d trials is %8.5f\n" P.N 4.0*(P.K / P.N);
