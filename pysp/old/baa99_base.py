import random
import bisect
import collections

from pyomo.core import *

def define_stochastic_model():

    model = ConcreteModel()

    model.d1_rhs = Param(mutable=True, initialize=100)
    model.d2_rhs = Param(mutable=True, initialize=100)

    model.PySP_StochasticParameters = Suffix()
    model.PySP_StochasticParameters[model.d1_rhs] = \
        [(17.75731865,0.04),
         (32.96224832, 0.04),
         (43.68044355, 0.04),
         (52.29173734, 0.04),
         (59.67893765, 0.04),
         (66.27551249, 0.04),
         (72.33076402, 0.04),
         (78.00434172, 0.04),
         (83.40733268, 0.04),
         (88.62275117, 0.04),
         (93.71693266, 0.04),
         (98.74655459, 0.04),
         (103.7634931, 0.04),
         (108.8187082, 0.04),
         (113.9659517, 0.04),
         (119.2660233, 0.04),
         (124.7925174, 0.04),
         (130.6406496, 0.04),
         (136.9423425, 0.04),
         (143.8948148, 0.04),
         (151.8216695, 0.04),
         (161.326406, 0.04),
         (173.7895514, 0.04),
         (194.0396804, 0.04),
         (216.3173937, 0.04)]

    model.PySP_StochasticParameters[model.d2_rhs] = \
        [(5.960319592, 0.04),
         (26.21044859, 0.04),
         (38.673594, 0.04),
         (48.17833053, 0.04),
         (56.10518525, 0.04),
         (63.05765754, 0.04),
         (69.35935045, 0.04),
         (75.20748263, 0.04),
         (80.73397668, 0.04),
         (86.03404828, 0.04),
         (91.18129176, 0.04),
         (96.2365069, 0.04),
         (101.2534454, 0.04),
         (106.2830673, 0.04),
         (111.3772488, 0.04),
         (116.5926673, 0.04),
         (121.9956583, 0.04),
         (127.669236, 0.04),
         (133.7244875, 0.04),
         (140.3210624, 0.04),
         (147.7082627, 0.04),
         (156.3195565, 0.04),
         (167.0377517, 0.04),
         (182.2426813, 0.04),
         (216.3173937, 0.04)]

    # define on the model after constraints are declared
    model.PySP_StochasticRHS = Suffix()

    return model

def generate_scenario_tree_model(scenario_count):

    from pyomo.pysp.util.scenariomodels import generate_simple_twostage

    st_model = generate_simple_twostage(scenario_count)

    first_stage = st_model.Stages.first()
    second_stage = st_model.Stages.last()

    # First Stage
    st_model.StageCostVariable[first_stage] = 'FirstStageCost'
    st_model.StageVariables[first_stage].add('x1')
    st_model.StageVariables[first_stage].add('x2')

    # Second Stage
    st_model.StageCostVariable[second_stage] = 'SecondStageCost'
    st_model.StageVariables[second_stage].add('v1')
    st_model.StageVariables[second_stage].add('v2')
    st_model.StageVariables[second_stage].add('u1')
    st_model.StageVariables[second_stage].add('u2')
    st_model.StageVariables[second_stage].add('w11')
    st_model.StageVariables[second_stage].add('w12')
    st_model.StageVariables[second_stage].add('w22')

    return st_model

def sample_into_scenario_tree_model(reference_model, scenario_tree_model):

    sample_data = {}
    for scenario_name in scenario_tree_model.Scenarios:

        # Draw random data
        sample_data[scenario_name] = []
        for param, table in reference_model.PySP_StochasticParameters.items():
            sample_data[scenario_name].append((ComponentUID(param), _sample(table)))

    return sample_data
