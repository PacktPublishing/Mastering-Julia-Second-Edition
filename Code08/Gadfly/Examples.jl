# Create scatter diagram

using Gadfly, Cairo, Fontconfig  
gd = Gadfly;
X = Float64[]; Y = Float64[];
N = 1000; 
for i in 1:N 
  x = rand() 
  y = rand() 
  if (x*x + y*y) < 1 
    push!(X,x) 
    push!(Y,y) 
  end 
end

dd = gd.plot(x = X, y = Y) 
draw(PNG("random-pts.png", 15cm, 12cm) , dd)

K = lastindex(X); # => 782 
println("Estimate of PI is ", round(4*K/N, digits=3)) 

#-------------------------------------------------------

# Load the GCSE results from RDatasets

using Gadfly, RDatasets, DataFrames;  
mlmf = dataset("mlmRev","Gcsemv")  
df = mlmf[completecases(mlmf), : ]  

# Data values for the exam and coursework differentiating between boys and girls

plot(df, x = "Course", y = "Written", color = "Gender")

#-------------------------------------------------------
#=  
Take note of the extensive use of Julia's broadcasting style. 
=# 
plot((x,y) ->  x .* exp.(-(x - floor.(x))).^2 .- y.^2, -8.0, 8.0, -2.0, 2.0)

#= 
  Display the following function: f(x) -> x*(5-x)*sin(5x) 
  A line plot of together with a scatter plot of the dithered points.
=#

import Gadfly 
const gd = Gadfly

x = collect(0.0:0.02:5.0); 
n = length(x);  # => Will be 251 points
y1 = [3*y*(5-y)*sin(5*y) for y in x]; 
y2 = Vector{Float64}(undef,n); 
[y2[i] = y1[i] + randn() for i in 1:n];

gd.plot(gd.layer(x=x,y=y1,gd.Geom.line, gd.Theme(default_color=gd.colorant"red")), 
gd.layer(x=x, y=y2, gd.Geom.point, gd.Theme(default_color=gd.colorant"blue")))


#-------------------------------------------------------

# Build a complex drawing based on the Sierpinski fractal

using Compose

function sierpinski(n)
  if n == 0
    compose(context(), polygon([(1,1), (0,1), (1/2, 0)]));
  else t = sierpinski(n - 1);
    compose( context(),(context( 1/4, 0, 1/2, 1/2), t),
                       (context( 0, 1/2, 1/2, 1/2), t),
                       (context( 1/2, 1/2, 1/2, 1/2), t));
  end
end

cxt = compose(sierpinski(1), linewidth(0.2mm),fill(nothing), stroke("black"));
draw(SVG("sierp1.svg", 10cm, 8.66cm), cxt);

cxt = compose(sierpinski(3), linewidth(0.2mm),fill(nothing), stroke("black"));
draw(SVG("sierp3.svg", 10cm, 8.66cm), cxt);

cxt = compose(sierpinski(5), linewidth(0.2mm),fill(nothing), stroke("black"));
draw(SVG("sierp5.svg", 10cm, 8.66cm), cxt);
