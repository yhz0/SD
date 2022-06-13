using CPLEX, JuMP
using TwoSD
using Distributions

# First we create a model to store the LP template
model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

# Stage 1 Variables
@variable(model, x[1:3] >= 0)

# Stage 2 Variables
@variables(model, begin
    s[1:3] >= 0
    f[1:3, 1:2] >= 0
end)

# Stage 1 Constraints
@constraint(model, budget, x[1] + x[2] + 4*x[3] <= 10)

# Stage 2 Constraints
# RHS randomness will be declared in the stoc part.
# RHS: omega[1:3] in link[1], link[2], link[3]
@constraint(model, link, s + f * [1, 1] .== 0)

# Random coefficient x1: -theta1
@constraint(model, capx1, -x[1] +f[1,1]+f[2,2]+f[3,2] <= 0)
# Random coefficient x2: -theta2
@constraint(model, capx2, -x[2] +f[1,2]+f[2,1]+f[3,2] <= 0)
@constraint(model, capx3, -x[3] +f[1,2]+f[2,2]+f[3,1] <= 0)

# The objective
@objective(model, Min, s[1] + s[2] + s[3])

# Specify where the first stage and second stage separates in the template.
split_position = Position(link[1], s[1])

# Demand distributions
dem_high = DiscreteNonParametric([0.0, 1, 2, 3, 5], [0.05, 0.2, 0.5, 0.2, 0.05])
dem_low = DiscreteNonParametric([0.0, 1, 2, 3], [0.1, 0.4, 0.4, 0.1])
link_up = Bernoulli(0.7)

function mystoc()
    return OneRealization(
        [
            # Demands are generated
            Position(link[1], "RHS") => rand(dem_high),
            Position(link[2], "RHS") => rand(dem_high),
            Position(link[3], "RHS") => rand(dem_low),

            # Is the link working? yes = 1 no = 0
            # Note negative sign in theta because it is on LHS.
            Position(capx1, x[1]) => -rand(link_up),
            Position(capx2, x[2]) => -rand(link_up)
        ]
    )
end

# By order of how the randomness are specified in mystoc(), we calculate the mean.
user_mean = [mean(dem_high), mean(dem_high), mean(dem_low), mean(link_up), mean(link_up)]

# Call SD
solution = solve_sd(model, split_position, user_mean, mystoc)

# Print solution if available
@show decision.(x, Ref(solution), AverageSolution)
@show decision.(x, Ref(solution), CompromiseSolution)
