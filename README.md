# Stochastic Decomposition 

Stochastic decomposition (SD) is a sequential sampling-based algorithm for two-stage stochastic linear programs (2-SLP). The algorithmic details can be found in:

1. Higle, J. L. and Sen, S. (1994). Finite master programs in regularized stochastic decomposition. Mathematical Programming, 67(1):143–168.
2. Sen, S. and Liu, Y. (2016). Mitigating uncertainty via compromise decisions in two-stage stochastic linear programming: Variance reduction. Operations Research, 64(6):1422–1437.
3. Gangammanavar, H., Liu, Y. and Sen, S. (2018) Stochastic Decomposition for Two-stage Stochastic Linear Programs with Random Cost Coefficients, availaable on Optimization Online.

This software is developed by Jason Mai, Lei Zhao, Yifan Liu, and Harsha Gangammanavar and Suvrajeet Sen.

### Support
Please report bugs [here on GitHub](https://github.com/USC3DLAB/SD/issues).

## Version information: 
`v1.0`: Supports 2-SLPs with randomness in right-hand side and the technology matrix.
`v2.0`: Supports 2-SLPs with randomness in right-hand side vector, the technology matrix and cost coefficient vector.

## Input file format
The implementation of accepts stochastic programs described using the SMPS file format:

* Core file: Describes the optimization problem corresponding to a single scenario
* Time file: Provides information regarding decomposition of the problem into master and subproblem
* Stoch file: Describes the stochastic information associated with the problem. (Currently, INDEP and BLOCK types are supported)

For more information about the SMPS file format [here](https://doi.org/10.1137/1.9780898718799.ch2)

## Installation

Prerequisite: CPLEX should be available on your machine. (If not, follow this link to get a trial version: [here](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/).

1). Download and compile SD with the following command. Simply copy and paste the following to your terminal (only support unix-like systems such as OS X and Ubuntu):

`git clone https://github.com/USC3DLAB/SD && cd sd/src && sudo make && cd ../instance && ln -s ../src/sd`

(note: In order to locate the path to cplex header and library file, "sudo make" is required as you can see. You can check the content of the makefile to verify this.)

2). Then excute sd by typing the following into your terminal:

`./sd`

At the prompt, enter an instance name for example: `pgp2`

(note: all instances are stored in the sdinput folder.)

Results will be stored in instance/sdoutput/pgp2. Please check pgp2.detailed_soln.out for solutions and time_sample.out for CPU time and number of samples used.
