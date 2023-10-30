
https://www.youtube.com/watch?v=hjmTax2DOHA

alpha = 0.75
beta = 0.1

du = f(u, p, t)dt + g(u, p, t)
f(u, p, t) = alpha * u
g(u, p, t) = beta * u

u0 = 1.0

t0 = 0.0
t1 = 1.0

tspan = (t0, t1)
step = 0.0001

prob_ODE = ODEProblem(f, u0, tspan)
prob_SDE = SDEProblem(f, g, u0, tspan)

sol_ODE = solve(prob_ODE)
sol_SDE = solve(prob_SDE, EM(), dt = step)

#--------------------------------------------------------------------------------

https://www.juliabloggers.com/summary-of-julia-plotting-packages
https://docs.juliahub.com/UnitfulRecipes/KPSlU/1.1.0/examples/2_Plots/

using DifferentialEquations

const g = 9.81                   # gravitational acceleration [m/s²]

function pendulum(du, u, t, p)
    du[1] = u[2]
    du[2] = -g*sin(u[1])
end

u0 = [3.0, 0.0];                  # initial state vector
tt = (0.0, 10.0);                 # time interval
ps = [1.0];					 	  # parameters

prob = ODEProblem(pendulum, u0, tt, ps)
sol = solve(prob)

# plot(sol)

θ = [sol.u[i][1] for i = 1:length(sol.u)]
ϕ = [sol.u[i][2] for i = 1:length(sol.u)]
τ = [sol.t[i] for i = 1:length(sol.u)]

using PythonPlot
plot(τ,θ)

#--------------------------------------------------------------------------------
