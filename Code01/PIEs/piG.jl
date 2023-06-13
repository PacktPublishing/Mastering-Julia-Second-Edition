using Printf

N =  100000
K = 0

sumsq(x,y) = x*x + y*y;

for i in 1:N
 if sumsq(rand(), rand()) < 1.0
    global K += 1
 end
end

@Printf.printf "Estimate of PI for %d trials is %8.5f\n" N 4.0*(K / N);

