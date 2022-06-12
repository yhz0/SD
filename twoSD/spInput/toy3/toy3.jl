# A toy problem for debugging T matrix randomness.

using JuMP, CPLEX


model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

@variable(model, x1 >= 0)
@variable(model, x2 >= 0)
@variable(model, y >= 0)

@constraint(model, c1, x1 <= 1)
@constraint(model, c2, x2 <= 1)

# x1: -xi1, x2:-xi2, RHS:xi0
@constraint(model, c3, x1+x2+y >= 0)

@objective(model, Min, y)

bm = backend(model)
CPLEX.CPXwriteprob(bm.env, bm.lp, "toy3.cor", "MPS")
CPLEX.CPXwriteprob(bm.env, bm.lp, "toy3.lp", "LP")