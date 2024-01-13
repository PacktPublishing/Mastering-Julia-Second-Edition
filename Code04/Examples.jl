#
# We can define these here or as required
using BenchmarkTools, MacroTools, Lazy, PythonPlot

# We need these as some functions have moved from Base to Stdlib
using Printf, Statistics, SpecialFunctions

recip(x::Number) =  (x == zero(typeof(x))) ? error("Invalid reciprocal") : one(typeof(x)) / x
recip(2)
recip(recip(2))
recip(11//17)
recip(11 + 17im)
aa = rand(3)

recip(aa)

## recip(a::Array) = [a[i] = recip(a[i]) for i = 1:length(a)]
## OR

recip(a::Array) = map(recip,a)
recip(aa)

map(sin,recip(aa)) # and can map functions to the array
bb = [2.1 3.2 4.3; 9.8 8.7 7.6]
recip(bb) # our definition works but returns a 2-D array

cc = recip(aa)'.*recip(bb)  # however we can do still do this.

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

incr(x) = x + 1
@code_native incr(2)
@code_native incr(2.7)
@code_native incr(2//7)
@code_native incr(2.0 + 7.0im)

# Look at some intermediate stages
dump(:(incr(2.7)))
@code_lowered incr(2.7)
@code_typed incr(2.7)
@code_llvm(incr(2.7))

ex1 = :((x^2 + y^2 - 2*x*y)^0.5)
ex2 = quote
    (x^2 + y^2 - 2*x*y)^0.5
end

# If we instantiate the variables x,y,z then we can evaluate the expression
# This is NOT a function call, so values need to be known before hand
#
x = 1.1; y = 2.5;
eval(ex1)

# Note that the @eval macro does NOT behave as the function call
@eval ex1

# Unless the expression is prefixed with '$'
@eval $ex1

# i.e.  @eval $(ex) === eval(ex)

# Although ex1 and ex2 are equivalent the AST of the expression 
# is slightly different
# Here it is for ex1, we leave ex2 to the reader
#
dump(ex1)

# Dump is very verbose, so using show_expr creates an S-expression version
#
Meta.show_sexpr(ex1)

# Traversing a tree
#
function traverse!(ex, symbols) end

function traverse!(ex::Symbol, symbols) 
    push!(symbols, ex) 
end

function traverse!(ex::Expr, symbols)
    if ex.head == :call  # function call
        for arg in ex.args[2:end]
            traverse!(arg, symbols)  # recursive
        end
    else
        for arg in ex.args
            traverse!(arg, symbols)  # recursive
        end
    end
end

function traverse(ex::Expr) 
    symbols = Symbol[]
    traverse!(ex, symbols) 
    return unique(symbols)  # Don't output duplicate
end

traverse(ex1)

macro pout(n)
    if typeof(n) == Expr 
       println(n.args)
    end
    return n
end
@pout x
@pout (x^2 + y^2 - 2*x*y)^0.5

macro dotimes(n, body)
    quote
        for i = 1:$(esc(n))
            $(esc(body))
        end
    end
end
@dotimes 3 print("Hi")
i = 0; @dotimes 3 [global i += 1; println(i*i)]

# Expand the @assert function
macroexpand(Main,:(@assert n > 0))

# Expand our @dotimes function, whic is somewhat simpler
#
macroexpand(Main,:(@dotimes 3 [global i += 1; println(i*i)]))

# Not all macro expansions produce short boiler plate
#
using Printf, Statistics
aa = [rand() for i = 1:100000];
@printf "The average value is %f.3 over %d trials"  mean(aa) length(aa)

# A more useful version of dotimes macro is until
# Creates a loop and breaks out when a condition fails
#
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
i = 0; @until (i >= 3) [global i += 1; println(i*i)]

# The until macro can be used to implement a simple if-then-else
#
macro iif(cond, body1, body2)
    :(if !$cond
        $(esc(body1))
    else
        $(esc(body2))
    end)
end


fac(n::Integer) = (n == 1) ? 1 : n*fac(n-1)
n = 10
@iif (n < 1) fac(n) ArgumentError("$n not positive")

n = -1; @iif (n < 1) fac(n) ArgumentError("$n not positive")

function fib(n::Integer)
  @assert n > 0
  (n == 1 || n == 2) ? 1 : fib(n-1) + fib(n-2)
end

macroexpand(Main,:(@time(fac(big(402)))))

# A timing macro

# First defined a modified Kempner function
# This converts very slowly

function kempner(n::Integer)
  @assert n > 0
  s = 0.0
  r9  = r"9"
  r9x = r"9{2}"
  for i in 1:n
    if (match(r9,string(i)) == nothing) || (match(r9x,string(i)) != nothing)
      s += 1.0/float(i)
    end
  end
  return s
end

[kempner(10^i) for i in 1:7]

macro bmk(fex, n::Integer)
    quote 
      let s = 0.0
      if $(esc(n)) > 0
        val = $(esc(fex))
        for i = 1:$(esc(n))
          local t0 = Base.time_ns()
          local val = $(esc(fex))
          s += Base.time_ns() - t0
        end
        return s/($(esc(n)) * 10e9)
      else
        Base.error("Number of trials must be positive")
      end
    end
  end
end

@bmk kempner(10^7) 10

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Evaluate a polynomial with a function
# This is NOT a macro version

function poly_native(x, a...)
  p=zero(x)
  for i = 1:length(a)
    p = p + a[i] *  x^(i-1)
  end
  return p
end

f_native(x) = poly_native(x,1,2,3,4,5)
f_native(2.1)

# Neither IS this!

function poly_horner(x, a...)
  b=zero(x)
  for i = length(a):-1:1
    b = a[i] + b * x
  end
  return b
end

# f -> (((5*x + 4)*x + 3)*x + 2)*x + 1
#
f_horner(x) = poly_horner(x,1,2,3,4,5)
f_horner(2.1)

# for the macro version we need a help function
# Define mad(x,a,b)
# [Julia has this function too :- muladd(x,a,b)]

mad(x,a,b) = a*x + b
mad(2.1,5,4)

# And NOW use Horner's method in a macro
# [ p... is a variable list of arguments passed as an array]

macro horner(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(mad(t, $ex, $(esc(p[i]))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end

@horner 2.1 1 2 3 4 5

# The saving of eliminating the loop for a 5th order polynomial is not great
# However for larger arrays and complex calculations this can be substantial
#
macroexpand(Main,:(@horner 2.1 1 2 3 4 5))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

macro reverse(ex)
    if isa(ex, Expr) && ex.head == :call
        return Expr(:call, ex.args[1], reverse(ex.args[2:end])...)
    else
        return ex
    end
end

@reverse  8//11 - 5//11 * 3//11


# i.e evaluates as :
3//11 * 5//11 - 8//11

macro q(s)     
  s0 = eval(s)
  try
    if length(s0) > 0
      ss = split(ss,'\n')
      for i = 1:length(ss)
         println(reverse(ss[i]))
      end
    end
  catch
    return  
  end
end

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The postwalk function splits an expression into symbols and then 
# reconstructs it, so we can apply different operations to each symbol

using MacroTools:postwalk

ex = :(1 + (2 + 3) + 4)
p = postwalk(ex) do x
    x isa Integer ? fac(x) : x
end

# Evaluate the expression
eval(p)

map(x -> @show(x), [1,2,3,4]);

postwalk(ex) do x
         @show x
end

@capture(ex, a_ + b_ + c_)
b

a*eval(b) + c    # => 1*5 + 4
reduce(+, 1:10)
plus(a, b) = :($a + $b)
p = reduce(plus, 1:10)

eval(p)

# Do something useful
# This is the series expansion for the SINE function
# We need the factorial function and will use the version from the STDLIB
#
using SpecialFunctions
k = 2
pp = [:($((-1)^k) * x^$(1+2k) / $(factorial(1+2k))) for k = 0:5]

# We can reduce this to a single expression
reduce(plus,pp)

# ... and evaluate it for a specific value of x
x = 2.1; eval(reduce(plus,pp))

using Lazy
import Lazy: cycle, range, drop, take

# Create a list of Fibonacci numbers and pick off first 20, using the @lazy macro.
#
fibs = @lazy 0:1:(fibs + drop(1, fibs));
take(20, fibs)

# Function style, pass the argument to a function
#
@> π/6 sin exp   # ==> exp(sin(π/6))

# Currying: The @> can also have functional arguments
#
f(x,μ) = -(x - μ)^2
@> π/6 f(1.6) exp

# The @>> macro reverse the order of the arguments
# Use this toutput the first 15 EVEN squares
#
esquares = @>> range() map(x -> x^2) filter(iseven);
take(15, esquares)

# Next create a list of primes
# The takewhile function is defined in Lazy
#
isprime(n) =
  @>> primes begin
    takewhile(x -> x<=sqrt(n))
    map(x -> n % x == 0)
    any; !
  end;

# We need to initialise the primes list
#
primes = filter(isprime, range(2));
take(20, primes)

isprime(113)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# A simple generated function to execute a multiply + add
# The result returned is a symbol
#
@generated function mad(a,b,c)
    Core.println("Calculating: a*b + c")
    return :(a * b + c)
end

# Call the function
mad(2.3,1.7,1.1)

# ... and again
# ... this time it is not re-evaluated
mad(2.3,1.7,1.1)

# With different types of arguments - the function is reevaluated
mad(2.3,1.7,1)

# To illustrate the use consider the following function ...
# ... which multiplies the size of dimensions of a n-D array
#
# First a 'normal' function definition

function pdims(x::Array{T,N}) where {T,N}
  s = 1
  for i = 1:N
    s = s * size(x, i)
  end
  return s
end

# And then a 'generated function' version
@generated function gpdims(x::Array{T,N}) where {T,N}
  ex = :(1)
  for i = 1:N
     ex = :(size(x, $i) * $ex)
  end
  return ex
end

# We need an array to test the function
aa = reshape([rand() for i = 1:1000],10,5,5,4); size(aa)

# Both functions provide the same result
#
pdims(aa) == gpdims(aa)

# But if we look at the lowered code
@code_lowered pdims(aa)

# But if we look at the lowered code
@code_lowered gpdims(aa)

# Which therefore results in very different generated native code
@code_native pdims(aa)

# This version does not have the looping
@code_native gpdims(aa)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Previous versions of Julia had an example section, absent in version v1.0
# An interesting one was modular integer, which required a little re-editing
#
module ModInts
export ModInt
import Base: +, -, *, /, inv

struct ModInt{n} <: Integer
  k::Int
  ModInt{n}(k) where n = new(mod(k,n))
end

Base.show(io::IO, k::ModInt{n}) where n =
    print(io, get(io, :compact, false) ? k.k : "$(k.k) mod $n")

(a::ModInt{n} + b::ModInt{n}) where n = ModInt{n}(a.k+b.k)
(a::ModInt{n} - b::ModInt{n}) where n = ModInt{n}(a.k-b.k)
(a::ModInt{n} * b::ModInt{n}) where n = ModInt{n}(a.k*b.k)
-(a::ModInt{n}) where n = ModInt{n}(-a.k)

inv(a::ModInt{n}) where n = ModInt{n}(invmod(a.k, n))
(a::ModInt{n} / b::ModInt{n}) where n = a*inv(b) 

Base.promote_rule(::Type{ModInt{n}}, ::Type{Int}) where n = ModInt{n}
Base.convert(::Type{ModInt{n}}, i::Int) where n = ModInt{n}(i)

end


# Test the module
#
using Main.ModInts
m1 = ModInt{11}(2)
m2 = ModInt{11}(7)
m3 = 3*m1 + m2  # => mod(13,11) => 2 

# Because of multiple dispatch we can do the following
#
mm = reshape([ModInt{11}(rand(0:10)) for i = 1:100],10,10)

ma = [ModInt{11}(rand(0:10)) for i = 1:10]

mm.*ma'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## An Ordered Pair Module
#
module OrdPairs
import Base: +,-,*,/,==,!=,>,<,>=,<=
import Base: abs,conj,inv,zero,one,show
import LinearAlgebra: transpose,adjoint,norm
export OrdPair
struct OrdPair{T<:Number}
    a::T
    b::T
end
OrdPair(x::Number) = OrdPair(x, zero(T))
value(u::OrdPair)   = u.a
epsilon(u::OrdPair) = u.b
zero(::Type{OrdPairs.OrdPair}) = OrdPair(zero(T),zero(T))
one(::Type{OrdPairs.OrdPair})  = OrdPair(one(T),zero(T))
abs(u::OrdPair)  = abs(value(u))
norm(u::OrdPair) = norm(value(u))
+(u::OrdPair, v::OrdPair) = OrdPair(value(u) + value(v), epsilon(u) + epsilon(v))
-(u::OrdPair, v::OrdPair) = OrdPair(value(u) - value(v), epsilon(u) - epsilon(v))
*(u::OrdPair, v::OrdPair) = OrdPair(value(u)*value(v), epsilon(u)*value(v) + value(u)*epsilon(v))
/(u::OrdPair, v::OrdPair) = OrdPair(value(u)/value(v),(epsilon(u)*value(v) - value(u)*epsilon(v))/(value(v)*value(v)))
==(u::OrdPair, v::OrdPair) = norm(u) == norm(v)
!=(u::OrdPair, v::OrdPair) = norm(u) != norm(v)
>(u::OrdPair, v::OrdPair)  = norm(u) > norm(v)
>=(u::OrdPair, v::OrdPair) = norm(u) >= norm(v)
<(u::OrdPair, v::OrdPair)  = norm(u) < norm(v)
<=(u::OrdPair, v::OrdPair) = norm(u) <= norm(v)
+(x::Number, u::OrdPair) = OrdPair(value(u) + x, epsilon(u))
+(u::OrdPair, x::Number) = x + u

-(x::Number, u::OrdPair) = OrdPair(x - value(u), epsilon(u))
-(u::OrdPair, x::Number) = OrdPair(value(u) - x, epsilon(u))

*(x::Number, u::OrdPair) = OrdPair(x*value(u), x*epsilon(u))
*(u::OrdPair, x::Number) = x*u

/(u::OrdPair, x::Number) = (1.0/x)*u

conj(u::OrdPair)  = OrdPair(value(u),-epsilon(u))
inv(u::OrdPair)   = one(OrdPair)/u

transpose(u::OrdPair) = u
transpose(uu::Array{OrdPairs.OrdPair,2}) = [uu[j,i] for i=1:size(uu)[1],j=1:size(uu)[2]]
adjoint(u::OrdPair) = u

convert(::Type{OrdPair}, x::Number) = OrdPair(x,zero(x))
promote_rule(::Type{OrdPair}, ::Type{<:Number}) = OrdPair

## show(io::IO,u::OrdPair) = print(io,value(u)," + (",epsilon(u),")ϵ")

function show(io::IO,u::OrdPair) 
 op::String = (epsilon(u) < 0.0) ? " - " : " + ";
 print(io,value(u),op,abs(epsilon(u)),"ϵ")
end

end

# Exercise the module
# Notice how the show() routine provides pretty-print output

using Main.OrdPairs
p1 = OrdPair(2.3,-1.7)
p2 = OrdPair(4.4,0.9)
p1 * p2
p2/p1

# Again we can operate on array of Ordered Pairs

using Statistics
pp = [OrdPair(rand(),rand()) for i in 1:100]

mean(pp)

# And promote a rational in a mixed OP
p3 = OrdPair(2.3, 11/7)

# The module is not a full implementation
p4 = OrdPair(2.3, 11.0 + 7.2im)

# Also a number of other arithmetic functions need to be imported from Base
std(p1)


