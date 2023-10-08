
using SymPy

u = symbols("u")
x = symbols("x", real=true)
y1, y2 = symbols("y1, y2", positive=true)
alpha = symbols("alpha", integer=true, positive=true)

typeof(x)

solve(x^2 + 1)

solve(u^2 + 1)

@vars x y
ex = x^2 + 2x + 1
subs(ex, x, y)

ex |> subs(x, 1)

p = x^2 + 3x + 2

factor(p)

expand(prod([(x-i) for i in 1:5]))

p = x^2 + 6x + 9; factor(p)

q = x*y + x*y^2 + x^2*y + x

collect(q, x)

collect(q,y)

r = 1/x + 1/x^2

together(r)

apart( (4x^3 + 21x^2 + 10x + 12) /  (x^4 + 5x^3 + 5x^2 + 4x))

theta = symbols("theta", real=true)

p = cos(theta)^2 + sin(theta)^2

simplify(sin(2theta) - 2sin(theta)*cos(theta))

# SymPy can solve polynomial equations
p = (x-3)^2*(x-2)*(x-1)*x*(x+1)*(x^2 + x + 1)

real_roots(p)

solve(p)

# We can use the Plot API to display graphs
x = symbols("x")

using Plots
plotly()

plot(x*sin(3.0x)*exp(-0.3x), 0, 10)

# This is the Basel series
#
i, n = symbols("i, n")
summation(i^2, (i, 1, n))

sn = Sum(1/i^2, (i, 1, n))
doit(sn)

limit(doit(sn), n, oo)

pp = 0.3 + x*sin(3.0x)*exp(-0.3x)
plot(pp, 0, 10)

# An initial valued problem
#
y = SymFunction("y")
a, x = symbols("a,x")
eqn = y'(x) - 3*x*y(x) - 1

x0, y0 = 0, 4
out = dsolve(eqn, x, (y, x0, y0))

x0, y0 = 0, a
out = dsolve(eqn, x, (y, x0, y0))

as = -2:0.6:2
ex = rhs(out)
p = plot(ex(a=>as[1]), -1.8, 1.8, ylims=(-4, 4))
[plot!(p, ex(a=>i), -1.8, 1.8, ylims=(-4, 4)) for i in as[2:end]]
p  


