using Printf

P =  [100000 0];

sumsq(x,y) = x*x + y*y;

for i in 1:P[1]
 if sumsq(rand(), rand()) < 1.0
    P[2] += 1
 end
end

@Printf.printf "Estimate of PI for %d trials is %8.5f\n" P[1] 4.0*(P[2] / P[1]);
