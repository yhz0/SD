from pyomo.core import *

def define_stochastic_model():

    model = ConcreteModel()

    model.d1_rhs = Param(mutable=True, initialize=0.0)
    model.d2_rhs = Param(mutable=True, initialize=0.0)

    d1_rhs_table = \
        [17.75731865,
         32.96224832,
         43.68044355,
         52.29173734,
         59.67893765,
         66.27551249,
         72.33076402,
         78.00434172,
         83.40733268,
         88.62275117,
         93.71693266,
         98.74655459,
         103.7634931,
         108.8187082,
         113.9659517,
         119.2660233,
         124.7925174,
         130.6406496,
         136.9423425,
         143.8948148,
         151.8216695,
         161.326406,
         173.7895514,
         194.0396804,
         216.3173937]

    d2_rhs_table = \
        [5.960319592,
         26.21044859,
         38.673594,
         48.17833053,
         56.10518525,
         63.05765754,
         69.35935045,
         75.20748263,
         80.73397668,
         86.03404828,
         91.18129176,
         96.2365069,
         101.2534454,
         106.2830673,
         111.3772488,
         116.5926673,
         121.9956583,
         127.669236,
         133.7244875,
         140.3210624,
         147.7082627,
         156.3195565,
         167.0377517,
         182.2426813,
         216.3173937]

    # define on the model after constraints are declared
    model.PySP_StochasticRHS = Suffix()

    return model, d1_rhs_table, d2_rhs_table

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
