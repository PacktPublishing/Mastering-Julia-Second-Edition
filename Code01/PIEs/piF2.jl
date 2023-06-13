using Printf

function compute_pi()
  sumsq(x,y) = x*x + y*y;
  N = 100000
  K = 0
  for i in 1:N
    if sumsq(rand(), rand()) < 1.0
      K += 1
    end
  end
  @Printf.printf "Estimate of PI for %d trials is %8.5f\n" N 4.0*(K / N);
end

compute_pi()
