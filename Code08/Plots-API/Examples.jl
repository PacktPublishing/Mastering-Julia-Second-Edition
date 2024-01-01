using Plots 
 
# Use GR as backend 
julia> gr() 
Plots.GRBackend()   

# [1]. Plot some simple sinusoidal curves

# There will be 25 data points in 3 separate series  
xs = 0 : 2π/25 : 2π;

# Define a sine, cosine and a scaled function(by x)  
data = [sin.(xs) cos.(xs) 0.5.*xs.*sin.(xs).*cos.(xs)];

# We put labels in a vector: one applies to each series  
labels = ["Sine" "Cosine" "Scaled function"]

# Marker shapes are to be applied to the data points 
markershapes = [:diamond :circle :star5]

# Marker colors in a matrix: applies to series and data 
markercolors = [:orange :red :blue]

# Now create the display 
julia> p = plot(xs, data, label = labels,  shape = markershapes,  
                markercolor = markercolors, markersize = 5)


# [2]. Multiple plots using layouts.

# The following is an example of a simple recipe:

mutable struct MyRecipe end  
@recipe  function f(::MyRecipe, n::Integer = 10; add_marker = false)  
  linecolor --> :blue 
  seriestype := :path  
  markershape --> (add_marker ? :circle : :none)  
  delete!(plotattributes, :add_marker)  
  rand(n)  
end

??????

plot(plot(mt, 25, linecolor = :black),  
     plot(mt, 100, linecolor = :red),  
     plot(mt, marker = (:star5,5)),  
     plot(mt, add_marker = true) )


# [3].Saving plots with HDF5

import Plots;  
hdf5();

p = plot(plot(mt, 25, linecolor = :black), plot(mt, 100, linecolor = :red)); 
Plots.hdf5plot_write(p, "plotsave.hdf5")

# Switch the backend to GR 
gr();
p = Plots.hdf5plot_read("plotsave.hdf5")
