### [ Basic Makie ] ###

# Define two target functions ...
 
using Makie
GLMakie.activate!(title="Plot title", fxaa=false);
 
 f1(u) = sin.(u) ./ (1 .+ u)
 f2(u) = u.*exp.(-0.5u) .* cos.(u)
 
#= ... and create some values for a specific range = 0:4pi, length = 80)f
   or both functions  =#

h = pi/20;
x = collect(0:h:4pi)
y1 = f1(x); 
y2 = f2(x);

#= Create the scene and the line plot of y1, overlaying it
   with the data points =#
   
scene = lines(x, y1, color = :blue)
scatter!(scene,x,y1,color = :red, markersize = 0.1)

#= Do the same for y2, using a different colour and marker,
   then display the finished "scene". =#

lines!(scene,x,y2, color = :black)
scatter!(scene, x, y2, color = :green, marker = :utriangle, markersize = 0.1)
scene

using RDatasets, DataFrames
iris = dataset("datasets", "iris")
scene = Scene()
colors = [:red, :green, :blue]

i = 1    # A color incrementOr
for sp in unique(iris[!,:Species])
  idx = iris[!,:Species] .== sp
  sel = iris[idx, [:SepalWidth, :SepalLength]]
  scatter!(scene, sel[:,1], sel[:,2], 
           color = colors[i],
           limits = Rect(1.5, 4.0, 3.0, 4.0))
   global i = i + 1
end

# Add the axes and label then show the scene
axis = scene[Axis]
axis[:names][:axisnames] = ("Sepal width","Sepal length")
scene

### [ Test images ] ###

using TestImages 
 
jet = testimage("jetplane"); 
size(jet) 
jetC = jet[121:320, :];

# To display the resulting image, using the image() function.

using GLMakie 
image(rotr90(jetC), axis=(aspect=DataAspect(), title="Leaving on a Jetplane ?"))

# Resizing

using Colors, TestImages
a17 = testimage("earth_apollo17"); 
typeof(a17[1,1]) 

a17g = Gray.(a17); 
typeof(a17g[1,1]) 

k = 6; kr = 1/k*k; 
sz = size(a17g); 
(m,n) = [convert(Integer,floor(sz[i]/k)) for i = 1:2];

# The resizing will be done by averaging over 36 pixels
# Define a zeroed array to hold the final values

a17r = zeros(Gray{Float32},m,n);
for i in 1:(m-1) 
  for j in 1:(n-1) 
    for h in 1:k 
       a17r[i,j] += kr*a17g[6*i + h, 6*j + h] 
     end 
  end 
end

