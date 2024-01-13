)using BenchmarkTools, UnicodePlots, PythonPlot

# Define a function to compute the sum of the squares of (x,y)

sumsq(x,y) = x*x + y*y;

# Set the number of trials and initialise the counter
# Using zero(Integer) will ensure the counter is an integer

N = 10^8;
K = zero(Integer);

using Printf

# Sum the number of pairs which lie within the inscribed circle

for i = 1:N
  if sumsq(rand(), rand()) < 1.0
    global K += 1
  end
end

# ... and output the value for PI

@Printf.printf "Estimate of PI for %d trials is %8.5f\n" N 4.0*(K / N);

# We can time this using the @time macro
# Note that his takes some time, will be different in the REPL

@time for i = 1:N
  if sumsq(rand(), rand()) < 1.0
    global K += 1
  end
end

# Now need the global K to access top level varables in the REPL
# The error message is less than informative
# i.e. UndefVarError: K not defined

# Try v0.7 to get a warning and understand the problem
# ┌ Warning: Deprecated syntax `implicit assignment to global variable `K``.
# │ Use `global K` instead.
# └ @ none:0

# This ONLY applies to the top level loops
# So wrapping in let/end, begin/end block is OK, since K is in local scope

# Has is by using using a function and passing K,N as parameters
#
# Here is a contrived example which demonstrates this

while (true)
    N = 1000000
    K = zero(Integer)
    for i = 1:N
        if sumsq(rand(), rand()) < 1.0
          K += 1
        end
    end
    println(4.0*(K/N))
    break
end


# We can rewrite as a generator expressions
# This does not require the variable 'K'
#
# Comprehensions can also be written without the enclosing square brackets,
# producing an object known as a generator.

((rand()^2 + rand()^2) < 1.0 for i = 1:N)

# The generator can be "evaluated" as the argument to an aggregate function

count((rand()^2 + rand()^2) < 1.0 for i = 1:N) * 4 / N

# And this runs much more quickly with a much better memory allocation

@time count((rand()^2 + rand()^2) < 1.0 for i = 1:N) * 4 / N

# JULIA_HOME is not longer defined in Version 1
# See: https://docs.julialang.org/en/v1/manual/environment-variables/
#
# JULIA_BINDIR sets Sys.BINDIR

Sys.BINDIR

# Rather than just using the @time macro we can get a more accuracy timing using the BenchmarkTools package

using BenchmarkTools

# Define the function to sum the series
# It converges very slowly

function basel(N::Integer)
  @assert N > 0
  s = 0.0
  for i = 1:N
    s += 1.0/float(i)^2
  end
  return s
end

basel(10^8)   # Evaluate it over 10^8 terms

π*π/6.0    # Compare this to the analytic solution

# Use the BenchmarkTools macro @benchmark to get more unbiased timings

@benchmark basel(10^8) samples=10

using UnicodePlots

# Generate an array of the numbers from 1 to 100
# The ordinate value is create using a list comprehensive

x = collect(1:100);
y = [x[i]*sin(0.3*x[i])*exp(-0.03*x[i]) for i = 1:length(x)];

myPlot = lineplot(x, y, title = "My Plot", name = "chp-01")

# Alternatively this can be done using a map() construction

t = collect(0.0:0.1:10.0);
y =  map(x -> x*sin(3.0*x)*exp(-0.3*x), t);
myPlot = lineplot(t, y, title = "Second Plot", name = "chp-01")

using PyPlot
plt = PyPlot   # Will get function name clashes with UnicodePlots.

S0  = 100;      # Spot price
K   = 102;      # Strike price
r   = 0.05;     # Risk free rate
q   = 0.0;      # Dividend yield
v   = 0.2;      # Volatility
tma = 0.25;     # Time to maturity
T   = 90;       # Number of time steps

dt  = tma/T;
N = T + 1;
x = collect(0:T);

# Qualify the function calls since UnicodePlot has exports same functions into Main.

plt.title("Asian Option trajectories");
plt.xlabel("Time");
plt.ylabel("Stock price");

for k = 1:5
  S = zeros(Float64,N)
  S[1] = S0;
  dW = randn(N)*sqrt(dt);
  [ S[t] = S[t-1] * (1 + (r - q - 0.5*v*v)*dt + v*dW[t] + 0.5*v*v*dW[t]*dW[t]) for t=2:N ]
  plt.plot(x,S)
end

