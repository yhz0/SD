from pyomo.core import *

import farmer_base

model = farmer_base.define_stochastic_model()

#
# Define the per-scenario random data realizations
#

ScenarioYield = {}
ScenarioYield['Scenario1'] = \
    {'WHEAT':2.0,'CORN':2.4,'SUGAR_BEETS':16.0}
ScenarioYield['Scenario2'] = \
    {'WHEAT':2.5,'CORN':3.0,'SUGAR_BEETS':20.0}
ScenarioYield['Scenario3'] = \
    {'WHEAT':3.0,'CORN':3.6,'SUGAR_BEETS':24.0}

def pysp_scenario_tree_model_callback():
    return farmer_base.generate_scenario_tree_model(len(ScenarioYield))

def pysp_instance_creation_callback(scenario_name, node_names):
    global ScenarioYield

    instance = model.clone()
    instance.Yield.store_values(ScenarioYield[scenario_name])

    return instance
