#! /usr/bin/env julia
#
using Getopt
short_list = "ho:q"
long_list  = ["help", "quiet", "output="]

println("\nSize of original ARGS array: ", size(ARGS))
println("Short list :: ", short_list)
println("Long list  :: ", long_list, "\n")

for (opt, arg) in getopt(ARGS, short_list, long_list)
   @show (opt, arg)
end
println("Length of modified ARGS array: ", size(ARGS))

