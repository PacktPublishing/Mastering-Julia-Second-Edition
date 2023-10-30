using Sundials, Plots

function exotherm(t, x, dx; n=1, a=1, b=1)
  p = x[2]^n * exp(x[1])
  dx[1] = p - a*x[1]
  dx[2] = -b*p
  return(dx)
end

t = collect(range(0.0; stop=5.0,length=1001))

fexo(t,x,dx) = exotherm(t, x, dx, a=0.6, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y11 = x[:,1];
y12 = x[:,2];

fexo(t,x,dx) = exotherm(t, x, dx, a=1.0, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y21 = x[:,1];
y22 = x[:,2];

fexo(t,x,dx) = exotherm(t, x, dx, a=1.2, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y31 = x[:,1];

y32 = x[:,2];

fexo(t,x,dx) = exotherm(t, x, dx, a=1.4, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y41 = x[:,1];
y42 = x[:,2];

fexo(t,x,dx) = exotherm(t, x, dx, a=1.6, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y51 = x[:,1];
y52 = x[:,2];

fexo(t,x,dx) = exotherm(t, x, dx, a=2.0, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y61 = x[:,1];
y62 = x[:,2];
plot!(t,y21)
plot!(t,y31)
plot!(t,y41)
plot!(t,y51)
plot!(t,y61)

#=
plot(t,y11,y12)
plot!(t,y21,y22)
plot!(t,y31,y32)
plot!(t,y41,y42)
plot!(t,y51,y52)
plot!(t,y61,y62)
=#

#-----------------------------------------------------------------

using Sundials, Plots

function exotherm(t, x, dx; n=1, a=1, b=1)
  p = x[2]^n * exp(x[1])
  dx[1] = p - a*x[1]
  dx[2] = -b*p
  return(dx)
end

fexo(t,x,dx) = exotherm(t, x, dx, a=0.6, b=0.1);
x = Sundials.cvode(fexo, [0.0, 1.0], t);
y = x[:,1];
plot(t,y)

#-----------------------------------------------------------------

using Sundials, PythonPlot

function exotherm(t, x, dx, n, a, b)
  p = x[2]^n * exp(x[1])
  dx[1] = p - a*x[1]
  dx[2] = -b*p
  return(dx)
end

t = collect(range(0.0; stop=5.0,length=1001));
bb = 0.1;

for aa in [0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
  fexo(t,x,dx) = exotherm(t, x, dx, 1, aa, bb);
  t = collect(range(0.0; stop=5.0,length=1001));
  x = Sundials.cvode(fexo, [0.0, 1.0], t);
  y = x[:,1];
  plot(t,y)
end















