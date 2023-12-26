"""
Calculate the sum of the series i/(i+1)^2 using the genie
function for an integer i in the range [1:n].
"""
function aladdin(n::Integer)
  @assert n > 0
  s = 0.0
  for i in 1:n
    s += genie(i,2)
  end
  return s
end

"""
Compute the value of the expression x/(x+1)^k" 
where x is a numeric and k a (non-complex) number.
""" 
genie(x,k) = x/(x+1)^k

const N_THIEVES = 40; # Ali Baba has 40 thieves in the fairy story.

"""
Compute and store the items using Aladdin's genie storing each partial sum
of the series i/(i+1)^2 in an array upto a value of 40.

So the preferable count to choose is a multiple of 40 in order to
capture the final value sum of the series in the array.
"""
function alibaba(n::Integer)
   @assert n >= N_THIEVES
   s = 0.0
   k = n/N_THIEVES
#=
We could define a fixed array but for only 41 items
i.e., Ali Baba and his 40 thieves s  push values instead 
=#
   a = Float64[]
   push!(a,s)
   for i in 1:n
     s += genie(i,2)
     (mod(i, k) == 0) && push!(a,s)
   end
   return a
end

"Try running panto() for list of functions"
function panto()
  helptxt = "Get help on individual pantomime characters: alibabi, aladdin, genie."
  println(helptxt);
end
