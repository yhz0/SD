using TwoSD

using Distributions, JuMP, CPLEX
using LinearAlgebra

# Data
N = 7
# holding cost
h = 1
# backorder cost
p = 4
# transshipment cost
c = 0.5

# Distribution:
mu = [100.0, 200, 150, 170, 180, 170, 170]
sigma = [20.0, 50, 30, 50, 40, 30, 50]

# IMPORTANT: We truncate the demand at plus or minus 3 sigmas to avoid negative demands
d_dist = product_distribution(
    truncated.(Normal.(mu, sigma), mu - 3*sigma, mu + 3*sigma)
)

# A direct model is needed to use TwoSD with Julia.
model = direct_model(CPLEX.Optimizer())

# Stage 1
@variable(model, s[1:N] >= 0)

# Stage 2
@variables(model, begin
    e[1:N] >= 0
    f[1:N] >= 0
    q[1:N] >= 0
    r[1:N] >= 0
    T[1:N,1:N] >= 0
end)

# @constraint(model, DIAGZERO[i=1:N], T[i,i] == 0)
@constraint(model, SUPPLY[i=1:N], f[i]+sum(T[i,j] for j = 1:N if j != i) + e[i] - s[i] == 0)
# RHS: d_i
@constraint(model, DEM[i=1:N], f[i]+sum(T[j,i] for j = 1:N if j != i) + r[i] == 0)
# RHS: sum(d)
@constraint(model, BAL, sum(r) + sum(q) == 0)
@constraint(model, END_INVENTORY[i=1:N], e[i] + q[i] - s[i] == 0)

# It is recommended that you formulate the problem in Min
@objective(model, Min, sum(h*e) + sum(c*T) + sum(p*r))

split_position = Position(SUPPLY[1], e[1])

using Random
rng = Random.MersenneTwister(1234)

function mystoc()
    d = rand(rng, d_dist)
    binding = [Position(DEM[i], "RHS") => d[i] for i in 1:N]
    push!(binding, Position(BAL, "RHS") => sum(d))
    return OneRealization(binding)
end

# We also need to construct the mean value of the random positions
# Note the binding is [d1, ..., dn, sum(d)]
user_mean = copy(mu)
push!(user_mean, sum(mu))

solve_sd(model, split_position, user_mean, mystoc)
