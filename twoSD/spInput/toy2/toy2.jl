using Distributions, JuMP, CPLEX
using TwoSD

a = 0.1

t1_dist = Normal(0, 1)
t2_dist = Normal(a, 1)

model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

@variable(model, x >= 0)
@variable(model, y)

@constraint(model, c1, x <= 1)
# RHS: t1, x: -t2
@constraint(model, c2, -x + y >= 0)

@objective(model, Min, y)

# TIM section
split_position = Position(c2, y)

# STOC section
function mystoc()::OneRealization
    t1 = rand(t1_dist)
    t2 = rand(t2_dist)

    return OneRealization([
        Position(c2, "RHS") => t1,
        Position(c2, x) => -t2
    ])
end

user_mean = [0.0, -a]

bmodel = backend(model)
CPLEX.CPXwriteprob(bmodel.env, bmodel.lp, "toy2.lp", "LP")
CPLEX.CPXwriteprob(bmodel.env, bmodel.lp, "toy2.mps", "MPS")

solution = solve_sd(model, split_position, user_mean, mystoc)
