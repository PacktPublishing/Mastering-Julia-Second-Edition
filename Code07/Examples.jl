
using LinearAlgebra
la = LinearAlgebra

A = [1 -2  2; 1 -1  2; -1  1  1];
det(A)

b = [5, 7, 5];
v = A\b

transpose(v)

A1 = A[:, 2:3] 

(A1\b)' 

A2 = A[1:2,:]; b2 = b[1:2];
(A2\b2)'

Alu = lu(A)

Alu.U # similar for the lower, i.e. Alu.L

Alu.L*Alu.U  # => A (original matrix)

Alu.P

#=
Suppose we have equations:

 x - 2y + 2z = 5
 x -  y + 2z = 3
-x +  y + z  = 6

=#

b = [5; 3; 6];

# Ax = b => LUx = b :  x = inv(U)*inv(L)*b

(x,y,z) = inv(Alu.U)*inv(Alu.L)*b

# Compute the eigenvalues of A
# (These are complex numbers)
U = eigvals(A)

# ... and the eigenvectors
V = eigvecs(A)

#=
The eigenvectors are the columns of the V matrix.
=#

A*V[:,1] - U[1]*V[:,1]

#=
That is, all the real and imaginary parts are of the e-16 order, 
so this is in effect a zero matrix of complex numbers.
=#

A = [1 -2  2; 1 -1  2; -1  1  1];
Diagonal(diag(A))

diag(A)

n = 2000;
B = randn(n,n);
B1 = B + B';
B2 = copy(B1);
B2[1,2] += 1eps();
B2[2,1] += 2eps();

issymmetric(B1)'
issymmetric(B2)'

@time eigvals(B1);
@time eigvals(B2);
@time eigvals(Symmetric(B2));

#---------------------------------------------------------------------------------------------------

# Signal Processing

using Plots

fq = 500.0;
N = 512;
T = 6 / fq;
t = collect(range(0, stop=T, length=N));

x1 = sin.(2π * fq * t);
x2 = cos.(8π * fq * t);
x3 = cos.(16π * fq * t);
x = x1 + 0.4*x2 + 0.2*x3;

Plots.plot(t,x)

## using Pkg; Pkg.add("FFTW")
using FFTW

X = rfft(x)

sr = N / T

fd = collect(range(0, stop = sr/2, length = div(N,2) + 1))

yy = abs.(X);

Plots.plot(fd[1:64], yy[1:64])

using DSP

ns = 0.1*randn(length(x));
xn = x + ns;
M = 16;
xm = ones(Float64, M) / M;
xf = conv(xn, xm)

Plots.plot(1:length(xf), xf)

responsetype = Lowpass(0.2)

prototype = Elliptic(4, 0.5, 30)


tf = tffilt(digitalfilter(responsetype, prototype))

tf = convert(TFFilter, digitalfilter(responsetype, prototype))



numerator_coefs   = coefb(tf)
denominator_coefs = coefa(tf)

responsetype = Bandpass(10, 40; fs=1000)
prototype = Butterworth(4)

xb = filt(digitalfilter(responsetype, prototype), x)›
plot(1:length(xb), xb)

# Image Processing

img = open("Files/lena.pgm");
magic  = chomp(readline(img));
params = chomp(readline(img));
pm = split(params)

# Remember the GSD 

try
  global wd = parse(Int64,pm[1]);
  global ht = parse(Int64,pm[2]);
catch
  error("Can't figure out the image dimensions")
end

# Version 1.0 way of defining a byte array
# readbytes!() will read in place

data = Array{UInt8,2}(undef,wd,ht)
readbytes!(img, data, wd*ht);

data = reshape(data,wd,ht);
close(img);


# Define a convolution mask

Gx = [1 2 1; 0 0 0; -1 -2 -1];
Gy = [1 0 -1; 2 0 -2; 1 0 -1];

dout = copy(data);
for i = 2:wd-1
  for j = 2:ht-1
    temp = data[i-1:i+1, j-1:j+1];
    x = sum(Gx.*temp)
    y = sum(Gy.*temp)
    p = Int64(floor(sqrt(x*x + y*y)))
    dout[i,j] = (p < 256) ? UInt8(p) : 0xff
  end
end

# ... and output the result
out = open("lenaX.pgm","w");
println(out,magic);
println(out,params);
write(out,dout);
close(out);

# This only works if you have Imagemagick (or similar) installed
run(`display lenaX.pgm`);

using Plots; gr()

# http://docs.juliadiffeq.org/latest/

#=
OrdinaryDiffEq.jl is part of the JuliaDiffEq common interface, 
but can be used independently of DifferentialEquations.jl. 

User passes to OrdinaryDiffEq.jl an algorithm to solve
=#

using OrdinaryDiffEq

function ff(d,u,p,t)
  d[1] =  u[1] - u[1]*u[2]
  d[2] = -u[2] + u[1]*u[2] - u[2]*u[3]
  d[3] = -u[3] + u[2]*u[3]
end

           
u0 = [0.5; 1.0; 2.0];    # Setup the initial conditions
tspan = (0.0,10.0);       # and the time range


#=
In OrdinaryDiffEq.jl, some good "go-to" choices for ODEs are:

AutoTsit5(Rosenbrock23()) handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation.
BS3() for fast low accuracy non-stiff.
Tsit5() for standard non-stiff. This is the first algorithm to try in most cases.
Vern7() for high accuracy non-stiff.
Rodas4() for stiff equations with Julia-defined types, events, etc.
radau() for really high accuracy stiff equations (requires installing ODEInterfaceDiffEq.jl)

=#

prob = ODEProblem(ff,u0,tspan)


u = solve(prob, Tsit5());

# Plot API will plot the array
styles = [:solid; :dash; :dot]
N = length(styles)
styles = reshape(styles, 1, N)  # styles is now a 1xN Vector

Plots.plot(u, line = (2,styles))

using Sundials
function exotherm(t, x, dx; n=1, a=1, b=1)
  p = x[2]^n * exp(x[1])
  dx[1] = p - a*x[1]
  dx[2] = -b*p
  return(dx)
end

t = collect(range(0.0; stop=5.0,length=1001))
fexo(t,x,dx) = exotherm(t, x, dx, a=0.6, b=0.1)
x1 = Sundials.cvode(fexo, [0.0, 1.0], t)

PyPlot.plot(x1[:,1])

#=
using Roots, Printf

f(x,a) = exp(x) - a*x

for p = 2.8:-0.02:2.6
  try
    ff(x) = f(x,p) 
    @printf "%.2f : %.5f\n" p find_zero(ff,1.0)
  catch
    error("No convergence for parameter value: $p")
  end
end

// There is no solution for a <= exp(1) => 2.7182 .....
=#



function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob,CVODE_Adams())
Plots.plot(sol,vars=(1,2,3))

using DifferentialEquations

f(du,u,p,t) = (du .= u)
g(du,u,p,t) = (du .= u)

u0 = rand(4,2);
W = WienerProcess(0.0,0.0,0.0);

prob = SDEProblem(f,g,u0,(0.0,1.0),noise=W)
sol = solve(prob,SRIW1());

# using Plots
# gr()
Plots.plot(sol)

function f(du,u,p,t)
  du[1] = u[1]
end

prob = ODEProblem(f,[0.2],(0.0,10.0))

rate(u,p,t) = 2
affect!(integrator) = (integrator.u[1] = integrator.u[1]/2)
jump = ConstantRateJump(rate,affect!)
jump_prob = JumpProblem(prob,Direct(),jump)

sol = solve(jump_prob,Tsit5())
Plots.plot(sol)

using Calculus

f(x)=sin(x)*cos(x)
derivative(f,1.0)

# Check since d(f) => cos*cos - sin*sin
cos(1.0)^2 - sin(1.0)^2

# Possible to curry the function
df = derivative(f)
df(1.0)

# Also defined is the 2nd derivative
d2f = second_derivative(f)
d2f(1.0)

# Also can use tick notation forhigher derivatives
# e.g. Value of the 3rd derivative
f'''(1.0)

#
# There are 2D functions, argument is a N-vector
# Be careful of name clashes
#
h(x) = (1+x[1])*exp(x[1])*sin(x[2])*cos(x[2])
gd=Calculus.gradient(h)
gd([1.0,1.0])

hs = Calculus.hessian(h)
hs([1.0,1.0])

# It is possible to output the symbolic version of the derivate
differentiate("sin(x)*cos(x)", :x)

# Not that clear but can be simplfied somewhat
# ... although not entirely perfect.
simplify(differentiate("sin(x)*cos(x)", :x)) 

# These techniques work more than a single variable too
# We get a 2-D array of partial derivatives
# Clearly the terms with a '0' multiplier can be ignored
#
simplify(differentiate("x*exp(-x)*sin(y)", [:x, :y]))

#---------------------------------------------------------------------------

struct D <: Number
  d1::Float64
  d2::Float64
end

import Base: +, /, convert, promote_rule
+(x::D, y::D) = D(x.d1+y.d1, x.d2+y.d2)
/(x::D, y::D) = D(x.d1/y.d1, (y.d1*x.d2 - x.d1*y.d2)/y.d1^2)
convert(::Type{D}, x::Real) = D(x,zero(x))
promote_rule(::Type{D}, ::Type{<:Number}) = D

da = D(17,1)

function hero(x; k::Integer = 10)
  @assert k > 0
  t = (1+x)/2
  for i = 2:k
    t = (t + x/t)/2
  end
  return t
end

Call the function with the D() structure.
db = hero(da)

db.d1^2   # Confirm the solution

How does this work?
Decode it without the 'dual' number

using Printf
function dhero(x; k = 10) 
    t  = (1+x)/2
    dt = 1
    @printf "%3d : %.5f : %.5f\n" 1 t dt
    for i = 2:k;  
        t  = (t+x/t)/2; 
        dt = (dt + (t - x*dt)/t^2)/2; 
        @printf "%3d : %.5f : %.5f\n" i t dt
    end    
    (t,dt)
end

dhero(17, k=5)


using SymPy

Display Hero's first 4 terms
Gets a little complex after that, try upping M!

M = 4
x = symbols("x")
for i = 1:M
  display(simplify(hero(x, k = i )))
end


#-----------------------------------------------------------------------------

using QuadGK
f(x) = sin(x)*(1.0 + cos(x))
quadgk(f,0.0,1.0)

#=
Addition function gives points and weights over the interval [-1,1]
Pick a sigmodial function : 1 - x exp(-|x|)
=#
g(u) = 1 - u*exp(-abs(u))

using PyPlot
x = collect(-1.0:0.1:1.0)
y = g.(x)
PyPlot.plot(x,y)

# The gauss() function creates a tuple of N array of points and weights
(x,w) = gauss(20)
s = sum([w[i]*g(x[i]) for i = 1:20])

# Alternate package which will do multidimension quadratutes
# Written by Steven Johnson (of PyCall, IJulia, etc. )
# It also does 1-D integration

using HCubature
hquadrature(f,0.0,1.0)

# But the power is the n-D quadratures
# Notice that the function expects an array argument

h(x) = 2.0*x[1]*exp(-x[1])*sin(x[2])*cos(x[2])
hcubature(h, [0,0], [1,1])



# Add JuMP#master to get v0.19

using JuMP, Clp

#=
Maximize the 5x+3y function subject to the constraint that: 3x+5y < 7
=#
m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 0 <= x <= 5 );
@variable(m, 0 <= y <= 10 );

@objective(m, Max, 5x + 3y );
@constraint(m, 2x + 5y <= 7.0 );

JuMP.optimize!(m)

println("x = ", JuMP.value(x), " y = ", JuMP.value(y))

JuMP.objective_value(m)

using JuMP, LinearAlgebra, Printf

N = 6;
m = Model()
@variable(m, x[1:N], Bin);      # Define array variable to hold results
profit = [ 5, 3, 2, 7, 4, 4 ];  # Profit vector of size N
weight = [ 2, 8, 4, 2, 5, 6 ];  # Weights vector of size 

maxcap = 15;

@objective(m, Max, dot(profit, x));
@constraint(m, dot(weight, x) <= maxcap);


using GLPK
JuMP.optimize!(m, with_optimizer(GLPK.Optimizer))

println("Objective is : ", JuMP.objective_value(m))
println("\nSolution is :")

for i = 1:N
    print("\tx[$i] = ", JuMP.value(x[i]))
    println(", p[$i]/w[$i] = ", profit[i]/weight[i])
end

using Optim

rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
result = Optim.optimize(rosenbrock, zeros(2), BFGS())

# Look at some values across a diagonal (x[1] == x[2]) cut
#
using Printf
for x in 0.0:0.1:1.3
  @printf "%3.1f : %f\n" x rosenbrock([x,x])
end



using SimJulia, ResumableFunctions
using Distributions
using Printf, Random

NUM_CUSTOMERS = 16     # total number of customers generated
NUM_TELLERS   = 2      # number of servers
QUEUE_MAX     = 2      # Maximum size of queue

μ = 0.4                # service rate
λ = 0.9                # arrival rate

arrival_dist = Exponential(1/λ)  # interarrival time distriubtion
service_dist = Exponential(1/μ)  # service time distribution

iseed = ccall((:clock,"libc"),Int32,())
Random.seed!(iseed);

queue_length = 0;
queue_stack  = Array{Integer,1}(undef,0);

@resumable function visit(env::Environment, 
                        teller::Resource, 
                        id::Integer, 
                        time_arrvl::Float64, 
                        dist_serve::Distribution)

# customer arrives
# recall the crazy scoping rules, we can see the array but not the scalar
#
    global queue_length
    @yield timeout(env, time_arrvl)
    @printf "Customer %2d %15s : %.3f\n" id "arrives" now(env)
    if queue_length > 0
        push!(queue_stack,id)
        println("CHECK: Length of the queue is $queue_length")
    end
    queue_length += 1
# customer starts to be served
    @yield request(teller)
    queue_length -= 1
    @printf "Customer %2d %15s : %.3f\n" id "being served" now(env) 
#  teller is busy
    @yield timeout(env, rand(dist_serve)) 
# customer leaves
    @yield release(teller) 
    @printf "Customer %2d %15s : %.3f\n" id "leaves" now(env)
end

# initialize simulation <: environment
sim     = Simulation()  

# initialize service resources
service = Resource(sim, NUM_TELLERS) 

# initialize customers and set arrival time
# customers arrive randomly baed on Poisson distribution
arrival_time = 0.0
for i = 1:NUM_CUSTOMERS 
    arrival_time += rand(arrival_dist)
    @process visit(sim, service, i, arrival_time, service_dist)
end

# run the simulation
run(sim)

# Check on which customers had to wait
queue_stack
















