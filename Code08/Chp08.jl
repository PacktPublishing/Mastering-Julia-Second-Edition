### [  Unicode Plots ] ###

#= 
Horizontal bar plot of the populations of the inner London Boroughs according 
to the 2021 census.  The data is contained in the file inner-london-boros.txt 
as comma-separated pairs of values. 
=#
import UnicodePlots
const ui = UnicodePlots 

X = Vector{String}(); 
Y = Vector{Int64}();

for s in eachline("inner-london-boros.txt")
  (x,y) = split(s,",") 
  push!(X,x)
  push!(Y, parse(Int64,y)) 
end

ui.barplot(X,Y, xlabel="Inner London Boroughs Population")

### [ Luxor ] ###
#= Use Luxor to display a drawing based on Sierpinski triangles =#

using Random, Colors
import Luxor
const lx = Luxor

#= 
Let's now define a function to create a triangle and another which uses this 
to build up the Sierpinski triangle 
=#

function triangle(points, degree)
  lx.sethue(cols[degree])
  lx.poly(points, :fill) 
 end

function sierpinski(points, degree)
  triangle(points, degree)
  if degree > 1
    p1, p2, p3 = points
    sierpinski([p1, lx.midpoint(p1, p2), lx.midpoint(p1, p3)], degree-1)
    sierpinski([p2, lx.midpoint(p1, p2), lx.midpoint(p2, p3)], degree-1)
    sierpinski([p3, lx.midpoint(p3, p2), lx.midpoint(p1, p3)], degree-1)
  end 
 end

# Now a function to draw the graphic and preview the result
function draw(n)
  lx.circle(lx.O, 100, :clip)
  points = lx.ngon(lx.O, 150, 3, -pi/2, vertices=true)
  sierpinski(points, n)
 end

# Set some parameters in Luxor
lx.Drawing(400, 250)
lx.background("white")
lx.origin()
depth = 8
cols = distinguishable_colors(depth)

# And call draw to display the result
draw(depth)
lx.finish()
lx.preview()

### [ Turtle graphics ] ###
#= 
Create simple drawing using “turtle” graphics - routines to control the turtle 
using commands such as: Forward, Turn, Circle, Orientation, Rectangle, Pendown, 
Penup, Pencolor, Penwidth, Reposition  
=#
import Luxor
lx = Luxor
lx.Drawing(600, 400, "turtles.png")
lx.origin(); lx.background("midnightblue");

tur = lx.Turtle();
lx.Pencolor(tur, "cyan");
lx.Penwidth(tur, 1.5);

n = 5;
for i in 1:400
  global n 
  lx.Forward(tur, n)
  lx.Turn(tur, 89.5)
  lx.HueShift(tur)
  n += 0.75
end 

lx.fontsize(20)
lx.finish()

### [ PyPlot ] ###
import PyPlot  
const py = PyPlot  

# Plot a modulated sinusiod
x = collect(range(0.0, stop=2pi, length=1000))  
y = sin.(3*x + 4*cos.(2*x));  
 
py.title("A sinusoidally modulated sinusoid");  
py.plot(x, y, color="red", linewidth=2.0, linestyle="--");  
py.savefig("sinusoid.svg");

# Plot a simple 3D surface
y = collect(range(0, stop=3π,length=250))  
py.surf(y, y, y.*sin.(y).*cos.(y)'.*exp.(-0.4y))
 
#=
  Create a plot using the XKCD "comic" mode from Python's Matplotlib.
  SO the xkcd() routine must be available in the (default) Python insatllation
=#
py.xkcd()  
x = collect(range(1, length=101, stop=10));  
y = sin.(3x + cos.(5x))  

py.title("XKCD fun") 
py.plot(x,y)

### [ Winston ] ###

# Plot simple functions over range [0, 3pi]
import Winston;  
const wn = Winston 
t = collect(range(0, stop=3pi, length=1000));

# Define 3 functions and create arrays based on the t variate and display them:
f(x::Array) = 8x .* exp.(-0.3x) .* sin.(3x);  
g(x::Array) = 0.1x.*(2pi .- x).*(3pi .- x);  
h(x::Array) = 10.0 ./ (1 .+ x.*x).^0.5;
y1 = f(t); y2 = g(t); y3 = h(t); 
wn.plot(t,y1,"b",t,y2,"r",t,y3,"k")
 
# Plot y1 against log(t)  
wn.semilogx(t,y1)  
wn.title("log(t) vs 10x * exp(-0.3x) * sin(3x)")

#=
  Use a "Framed Plot" for a more complex graph.
  Create and a linear relationship between two variables dithering 
  each by applying a random Gaussian variate
=#
p = wn.FramedPlot(aspect_ratio=1, xrange=(-10,110), yrange=(-10,110)); 

n = 21;  
x = collect(range(0.0, length=n, stop=100.0));  
yA = 10.0*randn(n) .+ 40.0;  
yB = x .+ 5.0*randn(n);

# Set labels and symbol styles for the plot 
a = wn.Points(x, yA, kind="circle");  
b = wn.Points(x, yB);  
wn.setattr(a,label="'a' points");  
wn.setattr(b,label="'b' points");  
wn.style(b, kind="filled circle");

# Plotine that 'fits' through the yB points ...
# ...and add a legend in the top LHS part of the graph
s = wn.Slope(1, (0,0), kind="dotted");  
wn.setattr(s, label="slope"); 
lg = wn.Legend(.1, .9, Any[a,b,s] );  
wn.add(p, s, a, b, lg);

# Display the completed graph
wn.display(p)

 


