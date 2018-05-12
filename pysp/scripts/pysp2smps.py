#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

import os
import gc
import sys
import time
import contextlib
import copy
from optparse import OptionParser, OptionGroup
try:
    import pstats
    pstats_available=True
except ImportError:
    pstats_available=False
# for profiling
try:
    import cProfile as profile
except ImportError:
    import profile
try:
    from guppy import hpy
    guppy_available = True
except ImportError:
    guppy_available = False

from pyutilib.pyro import shutdown_pyro_components
from pyutilib.services import TempfileManager

from pyomo.util import pyomo_command
from pyomo.core.base import maximize, minimize
from pyomo.opt.parallel import SolverManagerFactory

from pyomo.pysp.phutils import  _OLD_OUTPUT
from pyomo.pysp.scenariotree import ScenarioTreeInstanceFactory

import scenariotreeserverutils

#
# utility method to construct an option parser for ph arguments,
# to be supplied as an argument to the runph method.
#

def construct_smps_options_parser(usage_string):

    parser = OptionParser()
    parser.usage = usage_string

    # NOTE: these groups should eventually be queried from the PH,
    # scenario tree, etc. classes (to facilitate re-use).
    inputOpts        = OptionGroup(parser, 'Input Options')
    solverOpts       = OptionGroup(parser, 'Solver Options')
    outputOpts       = OptionGroup(parser, 'Output Options')
    otherOpts        = OptionGroup(parser, 'Other Options')

    parser.add_option_group( inputOpts )
    parser.add_option_group( solverOpts )
    parser.add_option_group( outputOpts )
    parser.add_option_group( otherOpts )

    inputOpts.add_option('-m','--model-directory',
      help='The directory in which all model (reference and scenario) definitions are stored. Default is ".".',
      action="store",
      dest="model_directory",
      type="string",
      default=".")
    inputOpts.add_option('-i','--instance-directory',
      help='The directory in which all instance (reference and scenario) definitions are stored. This option is required if no callback is found in the model file.',
      action="store",
      dest="instance_directory",
      type="string",
      default=None)
    def objective_sense_callback(option, opt_str, value, parser):
        if value in ('min','minimize',minimize):
            parser.values.objective_sense = minimize
        elif value in ('max','maximize',maximize):
            parser.values.objective_sense = maximize
        else:
            parser.values.objective_sense = None
    inputOpts.add_option('-o','--objective-sense-stage-based',
      help='The objective sense to use for the auto-generated scenario instance objective, which is equal to the '
           'sum of the scenario-tree stage costs. Default is None, indicating an Objective has been declared on the '
           'reference model.',
      action="callback",
      dest="objective_sense",
      type="choice",
      choices=[maximize,'max','maximize',minimize,'min','minimize',None],
      default=None,
      callback=objective_sense_callback)
    inputOpts.add_option('--bounds-cfgfile',
      help="The name of python script containing a pysp_boundsetter_callback function to compute and update scenario variable bounds. Default is None.",
      action="store",
      dest="bounds_cfgfile",
      default=None)
    inputOpts.add_option('--aggregate-cfgfile',
      help="The name of python script containing a pysp_aggregategetter_callback function to collect and store aggregate scenario data on the main scenario tree manager. Default is None.",
      action="store",
      dest="aggregate_cfgfile",
      type="string",
      default=None)
    inputOpts.add_option('--scenario-tree-seed',
      help="The random seed associated with manipulation operations on the scenario tree (e.g., down-sampling or bundle creation). Default is 0, indicating unassigned.",
      action="store",
      dest="scenario_tree_random_seed",
      type="int",
      default=0)
    inputOpts.add_option('--scenario-tree-downsample-fraction',
      help="The proportion of the scenarios in the scenario tree that are actually used. Specific scenarios are selected at random. Default is 1.0, indicating no down-sampling.",
      action="store",
      dest="scenario_tree_downsample_fraction",
      type="float",
      default=1.0)
    inputOpts.add_option('--scenario-bundle-specification',
      help="The name of the scenario bundling specification to be used when executing Progressive Hedging. Default is None, indicating no bundling is employed. If the specified name ends with a .dat suffix, the argument is interpreted as a filename. Otherwise, the name is interpreted as a file in the instance directory, constructed by adding the .dat suffix automatically",
      action="store",
      dest="scenario_bundle_specification",
      default=None)
    inputOpts.add_option('--create-random-bundles',
      help="Specification to create the indicated number of random, equally-sized (to the degree possible) scenario bundles. Default is 0, indicating disabled.",
      action="store",
      dest="create_random_bundles",
      type="int",
      default=None)
    inputOpts.add_option("--flatten-expressions", "--linearize-expressions",
      help="EXPERIMENTAL: An option intended for use on linear or mixed-integer models " \
           "in which expression trees in a model (constraints or objectives) are compacted " \
           "into a more memory-efficient and concise form. The trees themselves are eliminated. ",
      action="store_true",
      dest="flatten_expressions",
      default=False)

    solverOpts.add_option('--solver-manager',
      help="The type of solver manager used to coordinate scenario sub-problem solves. Default is serial.",
      action="store",
      dest="solver_manager_type",
      type="string",
      default="serial")
    solverOpts.add_option('--pyro-hostname',
      help="The hostname to bind on. By default, the first dispatcher found will be used. This option can also help speed up initialization time if the hostname is known (e.g., localhost)",
      action="store",
      dest="pyro_manager_hostname",
      default=None)
    solverOpts.add_option('--handshake-with-phpyro',
      help="When updating weights, xbars, and rhos across the PHPyro solver manager, it is often expedient to ignore the simple acknowledgement results returned by PH solver servers. Enabling this option instead enables hand-shaking, to ensure message receipt. Clearly only makes sense if the PHPyro solver manager is selected",
      action="store_true",
      dest="handshake_with_phpyro",
      default=False)
    solverOpts.add_option('--phpyro-required-workers',
      help="Set the number of idle phsolverserver worker processes expected to be available when the PHPyro solver manager is selected. This option should be used when the number of worker threads is less than the total number of scenarios (or bundles). When this option is not used, PH will attempt to assign each scenario (or bundle) to a single phsolverserver until the timeout indicated by the --phpyro-workers-timeout option occurs.",
      action="store",
      type=int,
      dest="phpyro_required_workers",
      default=None)
    solverOpts.add_option('--phpyro-workers-timeout',
     help="Set the time limit (seconds) for finding idle phsolverserver worker processes to be used when the PHPyro solver manager is selected. This option is ignored when --phpyro-required-workers is set manually. Default is 30.",
      action="store",
      type=float,
      dest="phpyro_workers_timeout",
      default=30)
    solverOpts.add_option('--shutdown-pyro',
      help="Shut down all Pyro-related components associated with the Pyro and PH Pyro solver managers (if specified), including the dispatch server, name server, and any solver servers. Default is False.",
      action="store_true",
      dest="shutdown_pyro",
      default=False)

    outputOpts.add_option('--symbolic-solver-labels',
      help='When interfacing with the solver, use symbol names derived from the model. For example, \"my_special_variable[1_2_3]\" instead of \"v1\". Useful for debugging. When using the ASL interface (--solver-io=nl), generates corresponding .row (constraints) and .col (variables) files. The ordering in these files provides a mapping from ASL index to symbolic model names.',
      action='store_true',
      dest='symbolic_solver_labels',
      default=False)
    outputOpts.add_option('--output-times',
      help="Output timing statistics during various runtime stages",
      action="store_true",
      dest="output_times",
      default=False)
    outputOpts.add_option('--output-instance-construction-times',
      help="Output timing statistics for instance construction (client-side only when using PHPyro",
      action="store_true",
      dest="output_instance_construction_times",
      default=False)
    outputOpts.add_option('--verbose',
      help="Generate verbose output for both initialization and execution. Default is False.",
      action="store_true",
      dest="verbose",
      default=False)
    otherOpts.add_option('--disable-gc',
      help="Disable the python garbage collecter. Default is False.",
      action="store_true",
      dest="disable_gc",
      default=False)
    if guppy_available:
        otherOpts.add_option("--profile-memory",
                             help="If Pympler is available (installed), report memory usage statistics for objects created after each PH iteration. A value of 0 indicates disabled. A value of 1 forces summary output after each PH iteration >= 1. Values greater than 2 are currently not supported.",
                             action="store",
                             dest="profile_memory",
                             type=int,
                             default=0)
    otherOpts.add_option('--profile',
      help="Enable profiling of Python code.  The value of this option is the number of functions that are summarized.",
      action="store",
      dest="profile",
      type="int",
      default=0)
    otherOpts.add_option('--traceback',
      help="When an exception is thrown, show the entire call stack. Ignored if profiling is enabled. Default is False.",
      action="store_true",
      dest="traceback",
      default=False)

    outputOpts.add_option('-d','--output-directory',
      help='The directory in which all SMPS related output files will be stored. ** required **',
      action="store",
      dest="output_directory",
      type="string",
      default=None)
    outputOpts.add_option('-b','--output-basename',
      help='The basename to use for all SMPS related output files. ** required **',
      action="store",
      dest="output_name",
      type="string",
      default=None)
    outputOpts.add_option('--implicit',
      help='Generate SMPS files using implicit parameter distributions.',
      action="store_true",
      dest="implicit",
      default=False)

    return parser


def DefaultOptions():
    parser = construct_smps_options_parser("")
    options, _ = parser.parse_args([''])
    return options

def GenerateScenarioTree(options,
                         scenario_instance_factory,
                         include_scenarios=None):

    scenario_tree = scenario_instance_factory.generate_scenario_tree(
        include_scenarios=include_scenarios,
        downsample_fraction=options.scenario_tree_downsample_fraction,
        bundles_file=options.scenario_bundle_specification,
        random_bundles=options.create_random_bundles,
        random_seed=options.scenario_tree_random_seed)

    #
    # print the input tree for validation/information purposes.
    #
    if options.verbose:
        scenario_tree.pprint()

    #
    # validate the tree prior to doing anything serious
    #
    if not scenario_tree.validate():
        raise RuntimeError("Scenario tree is invalid")
    else:
        if options.verbose:
            print("Scenario tree is valid!")

    return scenario_tree

class ScenarioTreeManager(object):

    def __init__(self, options):

        self._initialized = False
        self._options = copy.deepcopy(options)

        ### Updated During initialize()
        self._scenario_tree = None
        self._solver_manager = None
        self._instances = None
        self._bundle_binding_instance_map = {}
        self._bundle_scenario_instance_map = {}
        self._objective_sense = None
        # distributed worker information
        self._phpyro_worker_jobs_map = {}
        self._phpyro_job_worker_map = {}
        ###

        # callback information
        self._callback_function = {}
        self._callback_module_key = {}
        self._callback_module_name = {}
        # For the users to modify as they please in the aggregate
        # callback as long as the data placed on it can be serialized
        # by Pyro
        self._aggregate_user_data = {}

        # validate that callback functions exist in specified modules
        for options_attr_file, callback_name in (("aggregate_cfgfile",
                                                  "pysp_aggregategetter_callback"),
                                                 ("bounds_cfgfile",
                                                  "pysp_boundsetter_callback")):
            module_name = getattr(self._options, options_attr_file)
            if module_name is not None:
                sys_modules_key, module = load_external_module(module_name)
                callback = None
                for oname, obj in inspect.getmembers(module):
                    if oname == callback_name:
                        callback = obj
                        break
                if callback is None:
                    raise ImportError("PySP callback with name '%s' could "
                                      "not be found in module file: %s"
                                      % (callback_name, module_name))
                self._callback_function[callback_name] = callback
                self._callback_module_key[callback_name] = sys_modules_key
                self._callback_module_name[callback_name] = module_name

    def deactivate(self):
        import pyomo.solvers.plugins.smanager.phpyro

        if not self._initialized:
            return

        if self._solver_manager is not None:
            if isinstance(self._solver_manager,
                          pyomo.solvers.plugins.smanager.\
                          phpyro.SolverManager_PHPyro):
                scenariotreeserverutils.release_scenariotreeservers(self)
                self._solver_manager.release_workers()

                self._solver_manager.deactivate()

        self._solver_manager = None
        self._scenario_tree = None
        self._instances = None
        self._phpyro_worker_jobs_map = {}
        self._phpyro_job_worker_map = {}
        self._initialized = False

    def release(self):
        ret = (self._scenario_tree,
               self._solver_manager,
               self._phpyro_worker_jobs_map,
               self._phpyro_job_worker_map)
        self._solver_manager = None
        self._scenario_tree = None
        self._instances = None
        self._phpyro_worker_jobs_map = {}
        self._phpyro_job_worker_map = {}
        self._initialized = False
        return ret

    def initialize(self, scenario_tree):
        import pyomo.environ
        import pyomo.solvers.plugins.smanager.phpyro

        init_start_time = time.time()

        print("Initializing Scenario Tree Manager")
        print("")

        if scenario_tree is None:
            raise ValueError("A scenario tree must be supplied to the "
                             "ScenarioTreeManager initialize() method")
        self._scenario_tree = scenario_tree

        # construct the solver manager.
        if self._options.verbose:
            print("Constructing solver manager of type="
                  +self._options.solver_manager_type)
        self._solver_manager = \
            SolverManagerFactory(self._options.solver_manager_type,
                                 host=self._options.pyro_manager_hostname)
        if self._solver_manager is None:
            raise ValueError("Failed to create solver manager of "
                             "type="+self._options.solver_manager_type)

        isPHPyro =  isinstance(self._solver_manager,
                               pyomo.solvers.plugins.\
                               smanager.phpyro.SolverManager_PHPyro)

        if isPHPyro:

            if self._scenario_tree.contains_bundles():
                num_jobs = len(self._scenario_tree._scenario_bundles)
                if not _OLD_OUTPUT:
                    print("Bundle solver jobs available: "+str(num_jobs))
            else:
                num_jobs = len(self._scenario_tree._scenarios)
                if not _OLD_OUTPUT:
                    print("Scenario solver jobs available: "+str(num_jobs))

            workers_expected = self._options.phpyro_required_workers
            if (workers_expected is None):
                workers_expected = num_jobs

            timeout = self._options.phpyro_workers_timeout if \
                      (self._options.phpyro_required_workers is None) else \
                      None

            self._solver_manager.acquire_workers(workers_expected,
                                                 timeout)

        initialization_action_handles = []
        if isPHPyro:

            if self._options.verbose:
                print("Broadcasting requests to initialize "
                      "distributed scenario tree workers")

            initialization_action_handles.extend(
                scenariotreeserverutils.\
                initialize_scenariotree_workers(self))

            if self._options.verbose:
                print("Distributed scenario tree initialization "
                      "requests successfully transmitted")

        else:

            build_start_time = time.time()

            if not _OLD_OUTPUT:
                print("Constructing scenario tree instances")

            self._instances = \
                self._scenario_tree._scenario_instance_factory.\
                construct_instances_for_scenario_tree(
                    scenario_tree,
                    flatten_expressions=self._options.flatten_expressions,
                    report_timing=self._options.output_times,
                    preprocess=False)

            if self._options.verbose or self._options.output_times:
                print("Time to construct scenario instances=%.2f seconds"
                      % (time.time() - build_start_time))

            if not _OLD_OUTPUT:
                print("Linking instances into scenario tree")

            build_start_time = time.time()

            # with the scenario instances now available, link the
            # referenced objects directly into the scenario tree.
            self._scenario_tree.linkInInstances(
                self._instances,
                objective_sense=self._options.objective_sense,
                create_variable_ids=True)

            if self._options.verbose or self._options.output_times:
                print("Time link scenario tree with instances=%.2f seconds"
                      % (time.time() - build_start_time))

            if self._scenario_tree.contains_bundles():
                build_start_time = time.time()

                if self._options.verbose:
                    print("Forming binding instances for all scenario bundles")

                self._bundle_binding_instance_map.clear()
                self._bundle_scenario_instance_map.clear()

                if not self._scenario_tree.contains_bundles():
                    raise RuntimeError("Failed to create binding instances for scenario "
                                       "bundles - no scenario bundles are defined!")

                for scenario_bundle in self._scenario_tree._scenario_bundles:

                    if self._options.verbose:
                        print("Creating binding instance for scenario bundle=%s"
                              % (scenario_bundle._name))

                    self._bundle_scenario_instance_map[scenario_bundle._name] = {}
                    for scenario_name in scenario_bundle._scenario_names:
                        self._bundle_scenario_instance_map[scenario_bundle._name]\
                            [scenario_name] = self._instances[scenario_name]

                    # IMPORTANT: The bundle variable IDs must be idential to
                    #            those in the parent scenario tree - this is
                    #            critical for storing results, which occurs at
                    #            the full-scale scenario tree.

                    scenario_bundle._scenario_tree.linkInInstances(
                        self._instances,
                        create_variable_ids=False,
                        master_scenario_tree=self._scenario_tree,
                        initialize_solution_data=False)

                    bundle_ef_instance = create_ef_instance(
                        scenario_bundle._scenario_tree,
                        ef_instance_name=scenario_bundle._name,
                        verbose_output=self._verbose)

                    self._bundle_binding_instance_map[scenario_bundle._name] = \
                        bundle_ef_instance

                if self._output_times:
                    print("Scenario bundle construction time=%.2f seconds"
                          % (time.time() - build_start_time))

        # If specified, run the user script to collect aggregate
        # scenario data. This can slow down PH initialization as
        # syncronization across all phsolverservers is required
        if self._options.aggregate_cfgfile is not None:

            callback_name = "pysp_aggregategetter_callback"

            if isPHPyro:

                # Transmit invocation to phsolverservers
                print("Transmitting user aggregate callback invocations "
                      "to phsolverservers")
                if self._scenario_tree.contains_bundles():
                    for scenario_bundle in self._scenario_tree._scenario_bundles:
                        ah = scenariotreeserverutils.\
                             transmit_external_function_invocation_to_worker(
                                 self,
                                 scenario_bundle._name,
                                 self._callback_module_name[callback_name],
                                 callback_name,
                                 invocation_type=(scenariotreeserverutils.InvocationType.\
                                                  PerScenarioChainedInvocation),
                                 return_action_handle=True,
                                 function_args=(self._aggregate_user_data,))
                        while(1):
                            action_handle = self._solver_manager.wait_any()
                            if action_handle in initialization_action_handles:
                                initialization_action_handles.remove(action_handle)
                                self._solver_manager.get_results(action_handle)
                            elif action_handle == ah:
                                result = self._solver_manager.get_results(action_handle)
                                break
                        assert len(result) == 1
                        self._aggregate_user_data = result[0]

                else:
                    for scenario in self._scenario_tree._scenarios:
                        ah = scenariotreeserverutils.\
                             transmit_external_function_invocation_to_worker(
                                 self,
                                 scenario._name,
                                 self._callback_module_name[callback_name],
                                 callback_name,
                                 invocation_type=(scenariotreeserverutils.InvocationType.\
                                                  SingleInvocation),
                                 return_action_handle=True,
                                 function_args=(self._aggregate_user_data,))
                        while(1):
                            action_handle = self._solver_manager.wait_any()
                            if action_handle in initialization_action_handles:
                                initialization_action_handles.remove(action_handle)
                                self._solver_manager.get_results(action_handle)
                            elif action_handle == ah:
                                result = self._solver_manager.get_results(action_handle)
                                break
                        assert len(result) == 1
                        self._aggregate_user_data = result[0]

                # Transmit final aggregate state to phsolverservers
                print("Broadcasting final aggregate data to phsolverservers")
                initialization_action_handles.extend(
                    scenariotreeserverutils.transmit_external_function_invocation(
                        self,
                        "pyomo.pysp.ph",
                        "assign_aggregate_data",
                        invocation_type=(scenariotreeserverutils.InvocationType.\
                                         SingleInvocation),
                        return_action_handles=True,
                        function_args=(self._aggregate_user_data,)))

            else:

                print("Executing user aggregate getter callback function")
                for scenario in self._scenario_tree._scenarios:
                    result = self._callback_function[callback_name](
                        self,
                        self._scenario_tree,
                        scenario,
                        self._aggregate_user_data)
                    assert len(result) == 1
                    self._aggregate_user_data = result[0]

        # if specified, run the user script to initialize variable
        # bounds at their whim.
        if self._options.bounds_cfgfile is not None:

            callback_name = "pysp_boundsetter_callback"

            if isPHPyro:

                # Transmit invocation to phsolverservers
                print("Transmitting user bound callback invocations to "
                      "phsolverservers")
                if self._scenario_tree.contains_bundles():
                    for scenario_bundle in self._scenario_tree._scenario_bundles:
                        initialization_action_handles.append(
                            scenariotreeserverutils.\
                            transmit_external_function_invocation_to_worker(
                                self,
                                scenario_bundle._name,
                                self._callback_module_name[callback_name],
                                callback_name,
                                invocation_type=(scenariotreeserverutils.InvocationType.\
                                                 PerScenarioInvocation),
                                return_action_handle=True))
                else:
                    for scenario in self._scenario_tree._scenarios:
                        initialization_action_handles.append(
                            scenariotreeserverutils.\
                            transmit_external_function_invocation_to_worker(
                                self,
                                scenario._name,
                                self._callback_module_name[callback_name],
                                callback_name,
                                invocation_type=(scenariotreeserverutils.InvocationType.\
                                                 SingleInvocation),
                                return_action_handle=True))

            else:

                print("Executing user bound setter callback function")
                for scenario in self._scenario_tree._scenarios:
                    self._callback_function[callback_name](
                        self,
                        self._scenario_tree,
                        scenario)

        # gather scenario tree data if not local
        if isPHPyro:

            if self._options.verbose:
                print("Broadcasting requests to collect scenario tree "
                      "instance data from PH solver servers")

            scenariotreeserverutils.\
                gather_scenario_tree_data(self,
                                          initialization_action_handles)
            assert len(initialization_action_handles) == 0

            if self._options.verbose:
                print("Scenario tree instance data successfully collected")

            if self._options.verbose:
                print("Broadcasting scenario tree id mapping"
                      "to PH solver servers")

            scenariotreeserverutils.transmit_scenario_tree_ids(self)

            if self._options.verbose:
                print("Scenario tree ids successfully sent")

        self._objective_sense = \
            self._scenario_tree._scenarios[0]._objective_sense

        if self._options.verbose:
            print("Scenario tree manager is successfully initialized")

        if self._options.output_times:
            print("Overall initialization time=%.2f seconds"
                  % (time.time() - init_start_time))

        # gather and report memory statistics (for leak
        # detection purposes) if specified.
        if (guppy_available) and (self._options.profile_memory >= 1):
            print(hpy().heap())

        # indicate that we're ready to run.
        self._initialized = True

#
# Create a ScenarioTreeManager object from scratch using
# the options object.
#

def ScenarioTreeManagerFromScratch(options):

    start_time = time.time()
    if options.verbose:
        print("Importing model and scenario tree files")

    scenario_instance_factory = \
        ScenarioTreeInstanceFactory(options.model_directory,
                                    options.instance_directory,
                                    options.verbose)

    if options.verbose or options.output_times:
        print("Time to import model and scenario tree "
              "structure files=%.2f seconds"
              %(time.time() - start_time))

    try:

        scenario_tree = \
            GenerateScenarioTree(options,
                                 scenario_instance_factory)

    except:

        print("Failed to generate scenario tree")
        scenario_instance_factory.close()
        raise

    scenario_tree_manager = None
    try:

        scenario_tree_manager = ScenarioTreeManager(options)

        scenario_tree_manager.initialize(scenario_tree)

    except:

        print("A failure occurred while generating the "
              "scenario tree manager. Cleaning up...")
        if scenario_tree_manager is not None:
            scenario_tree_manager.deactivate()
        scenario_instance_factory.close()
        raise

    return scenario_tree_manager

#
# There is alot of cleanup that should be be done before a
# ScenarioTreeManager object goes out of scope (e.g.  releasing PHPyro
# workers, closing file archives, etc.). However, many of the objects
# requiring context management serve a purpose beyond the lifetime of
# the ProgressiveHedging object that references them. This function
# assumes the user does not care about this and performs all necessary
# cleanup when we exit the scope of the 'with' block.
# Example:
#
# with PHFromScratchManagedContext(options) as ph:
#    ph.run()
#

@contextlib.contextmanager
def ScenarioTreeManagerFromScratchManagedContext(options):

    scenario_tree_manager = None
    try:

        scenario_tree_manager = ScenarioTreeManagerFromScratch(options)
        yield scenario_tree_manager

    except:
        ScenarioTreeManagerCleanup(scenario_tree_manager)
        raise
    else:
        ScenarioTreeManagerCleanup(scenario_tree_manager)

def ScenarioTreeManagerCleanup(scenario_tree_manager):

    if scenario_tree_manager is None:

        return

    scenario_tree_manager.deactivate()

    if scenario_tree_manager._scenario_tree is not None:

        if scenario_tree_manager._scenario_tree._scenario_instance_factory is not None:

            scenario_tree_manager._scenario_tree._scenario_instance_factory.close()

#
# Convert a PySP scenario tree formulation to SMPS input files
#

def run_pysp2smps(options):
    import smpsutils

    if (options.output_directory is None):
        raise ValueError("Output directory name is required. "
                         "Use --output-directory command-line option")

    if (options.output_name is None):
        raise ValueError("Output base name is required. "
                         "Use --output-basename command-line option")

    if not os.path.exists(options.output_directory):
        os.makedirs(options.output_directory)

    start_time = time.time()

    io_options = {'symbolic_solver_labels':
                  options.symbolic_solver_labels}

    if options.implicit:

        if options.verbose:
            print("Importing model and scenario tree files")

        with ScenarioTreeInstanceFactory(
                options.model_directory,
                options.instance_directory,
                options.verbose) as scenario_instance_factory:

            if options.verbose or options.output_times:
                print("Time to import model and scenario tree "
                      "structure files=%.2f seconds"
                      %(time.time() - start_time))

            smpsutils.convert_implicit(options.output_directory,
                                       options.output_name,
                                       scenario_instance_factory,
                                       io_options=io_options)

    else:

        # This context manages releasing pyro workers and
        # closing file archives
        with ScenarioTreeManagerFromScratchManagedContext(options) \
             as scenario_tree_manager:
            smpsutils.convert_explicit(options.output_directory,
                                       options.output_name,
                                       scenario_tree_manager,
                                       io_options=io_options)

    end_time = time.time()

    print("")
    print("Total conversion execution time=%.2f seconds" %(end_time - start_time))

    scenario_tree_manager.deactivate()

#
# The main PH initialization / runner routine. Really only branches
# based on the construction source - a checkpoint or from scratch.
#

def exec_pysp2smps(options):
    import pyomo.environ

    start_time = time.time()

    try:

        run_pysp2smps(options)

    # This context will shutdown the pyro nameserver if requested.
    # Ideally, pyro workers can be reused without restarting the
    # nameserver
    finally:
        # if an exception is triggered, and we're running with
        # pyro, shut down everything - not doing so is
        # annoying, and leads to a lot of wasted compute
        # time. but don't do this if the shutdown-pyro option
        # is disabled => the user wanted
        if ((options.solver_manager_type == "pyro") or \
            (options.solver_manager_type == "phpyro")) and \
            options.shutdown_pyro:
            print("\n")
            print("Shutting down Pyro solver components.")
            shutdown_pyro_components(num_retries=0)

    print("")
    print("Total execution time=%.2f seconds"
          % (time.time() - start_time))

#
# the main driver routine for the pysp2smps script.
#

def main(args=None):
    #
    # Top-level command that executes the extensive form writer.
    # This is segregated from run_ef_writer to enable profiling.
    #

    #
    # Import plugins
    #
    import pyomo.environ
    #
    # Parse command-line options.
    #
    try:
        options_parser = construct_smps_options_parser("runph [options]")
        (options, args) = options_parser.parse_args(args=args)
    except SystemExit as _exc:
        # the parser throws a system exit if "-h" is specified - catch
        # it to exit gracefully.
        return _exc.code

    #
    # Control the garbage collector - more critical than I would like
    # at the moment.
    #

    if options.disable_gc:
        gc.disable()
    else:
        gc.enable()

    #
    # Run PH - precise invocation depends on whether we want profiling
    # output.
    #

    # if an exception is triggered and traceback is enabled, 'ans'
    # won't have a value and the return statement from this function
    # will flag an error, masking the stack trace that you really want
    # to see.
    ans = None

    if pstats_available and options.profile > 0:
        #
        # Call the main routine with profiling.
        #
        tfile = TempfileManager.create_tempfile(suffix=".profile")
        tmp = profile.runctx('exec_pysp2smps(options)',globals(),locals(),tfile)
        p = pstats.Stats(tfile).strip_dirs()
        p.sort_stats('time', 'cumulative')
        p = p.print_stats(options.profile)
        p.print_callers(options.profile)
        p.print_callees(options.profile)
        p = p.sort_stats('cumulative','calls')
        p.print_stats(options.profile)
        p.print_callers(options.profile)
        p.print_callees(options.profile)
        p = p.sort_stats('calls')
        p.print_stats(options.profile)
        p.print_callers(options.profile)
        p.print_callees(options.profile)
        TempfileManager.clear_tempfiles()
        ans = [tmp, None]
    else:
        #
        # Call the main routine without profiling.
        #

        if options.traceback:
            ans = exec_pysp2smps(options)
        else:
            try:
                try:
                    ans = exec_pysp2smps(options)
                except ValueError:
                    str = sys.exc_info()[1]
                    print("VALUE ERROR:")
                    print(str)
                    raise
                except KeyError:
                    str = sys.exc_info()[1]
                    print("KEY ERROR:")
                    print(str)
                    raise
                except TypeError:
                    str = sys.exc_info()[1]
                    print("TYPE ERROR:")
                    print(str)
                    raise
                except NameError:
                    str = sys.exc_info()[1]
                    print("NAME ERROR:")
                    print(str)
                    raise
                except IOError:
                    str = sys.exc_info()[1]
                    print("IO ERROR:")
                    print(str)
                    raise
                except pyutilib.common.ApplicationError:
                    str = sys.exc_info()[1]
                    print("APPLICATION ERROR:")
                    print(str)
                    raise
                except RuntimeError:
                    str = sys.exc_info()[1]
                    print("RUN-TIME ERROR:")
                    print(str)
                    raise
                except:
                    print("Encountered unhandled exception")
                    traceback.print_exc()
                    raise
            except:
                print("\n")
                print("To obtain further information regarding the "
                      "source of the exception, use the --traceback option")

    gc.enable()

    return ans

@pyomo_command('pysp2smps',
               'Convert a PySP Scenario Tree Formulation to SMPS input format')
def pysp2smps_main(args=None):
    return main(args=args)

if __name__ == "__main__":
    main(args=sys.argv[1:])
