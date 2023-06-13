sumsq(x,y) = x*x + y*y;

N = 1000000;
begin
K = 0
for i in 1:N
 if sumsq(rand(), rand()) < 1.0
    K += 1
 end
end
end

P = 4.0*(K / N);
println("Estimate of PI for $N trials is $P")

