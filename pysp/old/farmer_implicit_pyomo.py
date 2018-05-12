import numpy as np

from pyomo.core import *

import farmer_base

model = farmer_base.define_stochastic_model()

#
# Generate a random table for each parameter
#

marginals = {}
marginals['WHEAT'] = np.random.normal(2.5, 0.5, 1000)
marginals['CORN'] = np.random.normal(3.0, 0.4, 1000)
marginals['SUGAR_BEETS'] = np.random.normal(20.0, 4.0, 1000)

model.PySP_StochasticParameters = Suffix()
for i in model.CROPS:
    weights, points = np.histogram(marginals[i], bins=25, density=True)
    model.PySP_StochasticParameters[model.Yield[i]] = zip(points, weights)

#
# Define 
#

sample_data = None
def pysp_scenario_tree_model_callback():
    global sample_data

    scenario_tree_model = farmer_base.generate_scenario_tree_model(100)

    sample_data = farmer_base.sample_into_scenario_tree_model(model, scenario_tree_model)

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
