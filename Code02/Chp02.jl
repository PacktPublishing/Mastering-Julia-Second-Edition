# Get some packages we will need later.

using Statistics, Images, ImgView

x = 2;   typeof(x)
x = 2.0;  typeof(x)

x0 = 2^65
x1 = big(2)^65
@assert x0 == x1

for T = Any[Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128]
    println("$(lpad(T,7)): [$(typemin(T)),$(typemax(T))]")
end

x = 0xbb31; y = 0xaa5f;  xor(x,y)
x = 0xbb31;  x << 8
x = 0xbb31; p = (2 < 3); x + p

# Mean of 15 random numbers in range 0:100
#
A = rand(0:100,15)
mean(A)

# Broadcasting -- HERE ???
# Define a 2x3 matrix of rational numbers and convert to floats
X = convert.(Float64, [11/17 2//9 3//7; 4//13 5//11 6//23])

# Round to 4 decimal places
round.(X,digits=4)

# This does NOT work
f(x) = x*sin(3.0x)*exp(-0.03x)
Y = f(X)

# But this does
Y = f.(X)

# Of course this could also be done without the temporary function
Y = X .* sin.(3.0 .* X) .* exp.(- 0.03 .* X)

# First 15 Fibonnaci series

A = Array{Int64}(undef,15);
A[1]=1
A[2]=1
[A[i] = A[i-1] + A[i-2] for i = 3:length(A)]


# The 'recursive' definition of factorial function
# A simple loop is much quicker
#
function fac(n::Integer)
  @assert n > 0
  (n == 1) ? 1 : n*fac(n-1)
end

# This has difficulties with integer overflow
#
using Printf; pf=Printf;

for i = 1:30
  @pf.printf "%3d : %d \n" i fac(i)
end

# But since a BigInt <: Integer if we pass a BigInt the reoutine returns one
#
fac(big(30))

# We can check this using the gamma function
#
gamma(31)     # Γ(n+1)  <=>  n!


# This non-recursive one liner works!
# Note that this returns a BigInt regardless of the input
#
fac(N::Integer) =   (N < 1) ? throw(ArgumentError("N must be positive")) : reduce(*,collect(big.(1:N)))

@time(fac(402))

gamma(big(403.0))


# The 'standard' recursive definition

function fib(k::Integer)
  @assert k > 0
  (k < 3) ? 1 : fib(k-1) + fib(k-2)
end

@time fib(15)


# A better version

function fib(n::Integer)
  @assert n > 0
  a = Array{typeof(n)}(undef,n)
  a[1] = 1
  a[2] = 1
  for i = 3:n
    a[i] = a[i-1] + a[i-2]
  end
  return a[n]
end

@time(fib(big(402)))

# A still better version

function fib(n::Integer)
  @assert n > 0
  (a, b) = (big(0), big(1))
  while n > 0
    (a, b) = (b, a+b)
    n -= 1
  end
  return a
end

# THis converts to the Golden ratio
@printf "%.15f" fib(101)/fib(100)


# Is a built-in constant in Julia
Base.MathConstants.golden


# Reseed the random number generator: 
# rng will give a reproducible sequence of numbers
# if and only if a seed is provided

using Unicode
function bacs()
  bulls = cows = turns = 0
  a = Any[]
  while length(unique(a)) < 4 
    push!(a,rand('0':'9'))
    end
  my_guess = unique(a)
  println("Bulls and Cows")
  while (bulls != 4)
    print("Guess?> ")
    readline()
	if eof(STDIN)
      s = "q"
    else
      print("Guess?> ")
      s = chomp(readline())
    end 
    if (s == "q")
      print("My guess was "); [print(my_guess[i]) for i=1:4]
      return
    end
    guess = collect(s)
    if !(length(unique(guess)) == length(guess) == 4 && all(Unicode.isdigit,guess))
      print("\nEnter four distinct digits or q to quit: ")
      continue
    end
    bulls = sum(map(==, guess, my_guess))
    cows = length(intersect(guess,my_guess)) - bulls
    println("$bulls bulls and $cows cows!")
    turns += 1
  end
  println("\nYou guessed my number in $turns turns.")
end

bacs()


# Look at different definitions of the norm function
# For a Gaussian distribution of size N we should expect the answer ~= √N
# The first call f1(1) is to run in the function and not affext the timing
# This version uses the function in Base

import LinearAlgebra:norm 
f1(n) = norm(randn(n))

f1(10);

@time f1(100_000_000)


# We can get the same result using a mapreduce procedure
# Note that it is a new set of random number, so the answer is slightly different
# The time is about the same

f2(n) = sqrt(mapreduce(x -> x*x, +, randn(n)))
f2(10);
@time f2(100_000_000)


# Using a conventional mapping we need to pipe the result to sum it 
# and then take the square root
# This takes a little longer than the previous.

f3(n) = map(x -> x*x,randn(n)) |> sum |> sqrt
f3(10);
@time f3(100_000_000)

# Finally we can non-vectorize the code, which is much quicker,
# In Julia non-vectorized (i.e loopy) code is invariably faster
# han the vectorized equivalent.

function f4(n)
    t = 0.0
    for i = 1:n
        t += randn()^2
    end
    return sqrt(t)
end

f4(10);
@time f4(100_000_000)

# Generate primes with Sieve of Eratoshenes
#
# Define a helper function and then the actual sieve.
        
cop(X, i) = any(j -> i % j == 0, X)        

function erato(N::Integer)
  @assert N > 0
  P = Int[]
  for i in 2:N
    if !cop(P, i)
      push!(P, i)
    end
  end
  return P
end
        
# Run it for the number of primes upto 1 million        

tm = @elapsed A = erato(1_000_000);
print("Computed $(length(A)) primes in $(round(tm, digits=4)) sec.")

        
# Julia Sets
# Define high level function to create a Julia set

function juliaset(z, z0, nmax::Int64)
    for n = 1:nmax
        if abs(z) > 2 (return n-1) end
        z = z^2 + z0
    end
    return nmax
end

# And how to write a PGN file to disk.

function create_pgmfile(img, outf::String)
    s = open(outf, "w")
    write(s, "P5\n")    
    n, m = size(img)
    write(s, "$m $n 255\n")
    for i=1:n, j=1:m
        p = min(img[i,j],255)
        write(s, UInt8(p))
    end
    close(s)
end

# Run function and create a 800x400 image
# Redefine the c0 parameter to create a different set
            
h = 400; 
w = 800; 
m = Array{Int64,2}(undef,h,w);

c0 = -0.8 + 0.16im;
pgm_name = "jset.pgm";

t0 = time();
for y=1:h, x=1:w
    c = complex((x-w/2)/(w/2), (y-h/2)/(w/2))
    m[y,x] = juliaset(c, c0, 256)
end
t1 = time();

# You should find the file in the same chapter as the notebook

create_pgmfile(m, pgm_name);
@printf "Written %s\nFinished in %.4f seconds.\n" pgm_name (t1-t0);


# Display the image using Images package

img = load("jset.pgm");
imgshow(img)

