import itertools

from pyomo.core import *

import baa99_base

model, d1_rhs_table, d2_rhs_table = baa99_base.define_stochastic_model()

model.x1 = Var(bounds=(0,217))
model.x2 = Var(bounds=(0,217))
model.v1 = Var(within=NonNegativeReals)
model.v2 = Var(within=NonNegativeReals)
model.u1 = Var(within=NonNegativeReals)
model.u2 = Var(within=NonNegativeReals)
model.w11 = Var(within=NonNegativeReals)
model.w12 = Var(within=NonNegativeReals)
model.w22 = Var(within=NonNegativeReals)

model.FirstStageCost = \
    Expression(initialize=(4*model.x1 + 2*model.x2))
model.SecondStageCost = \
    Expression(initialize=(-8*model.w11 - 4*model.w12 - 4*model.w22 +\
                           0.2*model.v1 + 0.2*model.v2 + 10*model.u1 + 10*model.u2))

model.obj = Objective(expr=model.FirstStageCost + model.SecondStageCost)

model.d1 = Constraint(expr=model.w11 + model.u1 == model.d1_rhs)
model.PySP_StochasticRHS[model.d1] = True

model.d2 = Constraint(expr=model.w12 + model.w22 + model.u2 == model.d2_rhs)
model.PySP_StochasticRHS[model.d2] = True

model.s1 = Constraint(expr=-model.x1 + model.w11 + model.w12 + model.v1 == 0)

model.s2 = Constraint(expr=-model.x2 + model.w22 + model.v2 == 0)

num_scenarios = len(d1_rhs_table) * len(d2_rhs_table)
def pysp_scenario_tree_model_callback():

    scenario_tree_model = \
        baa99_base.generate_scenario_tree_model(num_scenarios)

    return scenario_tree_model

scenario_data = itertools.product(d1_rhs_table,
                                  d2_rhs_table)
def pysp_instance_creation_callback(scenario_name, node_names):

    #
    # Clone a new instance and update the stochastic
    # parameters from the sampled scenario
    #

    instance = model.clone()

    d1_rhs_val, d2_rhs_val = next(scenario_data)
    instance.d1_rhs.value = d1_rhs_val
    instance.d2_rhs.value = d2_rhs_val

    return instance
