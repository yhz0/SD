using CPLEX, JuMP
using Distributions

# Demand distributions
dem_high = DiscreteNonParametric([0.0, 1, 2, 3, 5], [0.05, 0.2, 0.5, 0.2, 0.05])
dem_low = DiscreteNonParametric([0.0, 1, 2, 3], [0.1, 0.4, 0.4, 0.1])
link_up = Bernoulli(0.7)

# scenario generation
xi_tilde = [dem_high, dem_high, dem_low, link_up, link_up]

xi = vec(collect(Base.product(support.(xi_tilde)...)))
prob_list = vec(collect(Base.product(probs.(xi_tilde)...)))
prob = prod.(prob_list)

# Deterministic
# xi = [mean.(xi_tilde)]
# prob = [1.0]

# Check all the scenarios sum to one
@assert(sum(prob) ≈ 1)

# Total number of scenarios
N = length(prob)

# Build model
model = direct_model(CPLEX.Optimizer())

@variables(model, begin
# Stage 1 Variables
    x[1:3] >= 0
# Stage 2 Variables
    s[1:N, 1:3] >= 0
    f[1:N, 1:3, 1:2] >= 0
end)

@constraint(model, budget, x[1] + x[2] + 4*x[3] <= 10.0)


# Note: xi[i] = (omega_1[i], omega_2[i], omega_3[i], theta_1[i], theta_2[i]),
# where i = 1..N is the scenario number
@constraint(model, link1[i = 1:N], s[i, 1] + f[i,1,1] + f[i,1,2] == xi[i][1])
@constraint(model, link2[i = 1:N], s[i, 2] + f[i,2,1] + f[i,2,2] == xi[i][2])
@constraint(model, link3[i = 1:N], s[i, 3] + f[i,3,1] + f[i,3,2] == xi[i][3])

@constraint(model, capx1[i = 1:N], f[i,1,1]+f[1,2,2]+f[i,3,2] <= xi[i][4] * x[1])
@constraint(model, capx2[i = 1:N], f[i,1,2]+f[i,2,1]+f[i,3,2] <= xi[i][5] * x[2])
@constraint(model, capx3[i = 1:N], f[i,1,2]+f[i,2,2]+f[i,3,1] <= x[3])

@objective(model, Min, sum(prob[i] * (s[i,1] + s[i,2] + s[i,3]) for i in 1:N))

# Optimization
optimize!(model)

# Assert problem is solved
@assert(termination_status(model) == OPTIMAL)

@show value.(x)
@show objective_value(model)

# fix.(x, [5.0, 5.0, 0.0]; force=true)
# optimize!(model)
# objective_value(model)