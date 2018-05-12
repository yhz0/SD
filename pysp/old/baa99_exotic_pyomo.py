import random

from pyomo.core import *

import baa99_base

#
# A more exotic form of the baa99 model used for
# testing edge cases:
#  - Objective has a constant term
#  - Uses range constraints
#  - Variable bounds are moved to constraints
#    in order to create "first-stage" constraints

model = baa99_base.define_stochastic_model()

model.x1 = Var()
model.x2 = Var()
model.v1 = Var()
model.v2 = Var()
model.u1 = Var()
model.u2 = Var()
model.w11 = Var()
model.w12 = Var()
model.w22 = Var()

model.FirstStageCost = \
    Expression(initialize=(4*model.x1 + 2*model.x2 + 1.0))
model.SecondStageCost = \
    Expression(initialize=(-8*model.w11 - 4*model.w12 - 4*model.w22 +\
                           0.2*model.v1 + 0.2*model.v2 + 10*model.u1 + 10*model.u2))

model.obj = Objective(expr=model.FirstStageCost + model.SecondStageCost)

model.s2 = Constraint(expr=0 <= -model.x2 + model.w22 + model.v2 <= 0)

model.d2 = Constraint(expr=model.d2_rhs <= model.w12 + model.w22 + model.u2 <= model.d2_rhs)
model.PySP_StochasticRHS[model.d2] = model.d2_rhs

model.s1 = Constraint(expr=0 <= -model.x1 + model.w11 + model.w12 + model.v1 <= 0)

model.d1 = Constraint(expr=model.d1_rhs <= model.w11 + model.u1 <= model.d1_rhs)
model.PySP_StochasticRHS[model.d1] = model.d1_rhs

model.bounds = ConstraintList(noruleinit=True)
model.bounds.add(0 <= model.x1 <= 217)
model.bounds.add(0 <= model.x2 <= 217)
model.bounds.add(0 <= model.v1)
model.bounds.add(0 <= model.v2)
model.bounds.add(0 <= model.u1)
model.bounds.add(0 <= model.u2)
model.bounds.add(0 <= model.w11)
model.bounds.add(0 <= model.w12)
model.bounds.add(0 <= model.w22)

num_scenarios = 100
sample_data = None
def pysp_scenario_tree_model_callback():
    global sample_data

    scenario_tree_model = \
        baa99_base.generate_scenario_tree_model(num_scenarios)

    sample_data = baa99_base.sample_into_scenario_tree_model(model, scenario_tree_model)

    return scenario_tree_model

def pysp_instance_creation_callback(scenario_name, node_names):
    assert sample_data is not None

    #
    # Clone a new instance and update the stochastic
    # parameters from the sampled scenario
    #

    instance = model.clone()

    for param_cuid, param_value in sample_data[scenario_name]:
        param_data = param_cuid.find_component(instance)
        param_data.value = param_value

    return instance
