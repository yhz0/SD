# TwoSD Julia Interface

## What type of problems can SD solve?

TwoSD solves two-stage stochastic linear programs in the following form.

$$
\begin{align*}
\min \quad & c^\top x + E_\xi[h(x, \xi)]\\
s.t.\quad & Ax = b, x \geq 0
\end{align*}
$$

$$
h(x,\xi) = \min \{d(\xi)^\top y: T(\xi)x + Wy = r(\xi), \quad y\geq 0\}
$$

### Important Assumptions on the Problem

1. The problem should satisfy the __relatively complete recourse property__.
2. As a minimization problem, the second stage objective function must be uniformly bounded below for $\xi$ almost surely and for all $x$.
3. The user need to supply a subroutine/function that generates independent realizations (samples) of the random LP coefficients each time it is called. Although the random variables do not have to be finite supported, the structure should not change. (see the Sampling Subroutine)

## Installation

0. Clone or download the repository.

1. Compile SD

    You will need the following to compile SD, which is written in C.

    - CPLEX

        CPLEX is a proprietary MIP/LP solver. Get an academic version.

        Windows users: It is strongly recommended that you __do NOT install CPLEX into a path that contains SPACE or special non-ASCII symbols__, specifically do not use the default path "Program Files".

    - make
    - gcc

        Windows users: only tested under MSYS2.

2. Specify in the Makefile where the CPLEX installation is.

    ```makefile
    CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio201/cplex
    ```

    If you installed CPLEX version other than 20.10, change the path accordingly.

    (The following is for Windows build only.)

    ```makefile
    LIBS := -lcplex2010
    ```

    Change the library name as necessary, if you are using a CPLEX version other than 20.10.

3. Run "make all" under "twoSD" __subfolder__.

4. Install the package using Julia's package manager. To avoid CPLEX issues, make sure you have set up CPLEX and JuMP correctly in Julia before you install this package. Open a terminal or shell, change directory to where you placed the repository. Run julia in that terminal, in Julia's package manager prompt (by pressing "]" key in REPL), run `dev "."`. 
This will add Julia as a development package.

## How to write the input?

To formulate an SP for this interface, we need three components:
1. The LP "template" that determines the fixed part of the SP coefficients.
2. The splitting position, i.e. the row and column that separates the first and second stages.
3. A subroutine that generates realization of the random part of LP coefficients.

In many SP problems, only a small portion of the coefficients are truly random. For perfomance reasons, we partition the second stage LP coefficients into two parts: the fixed part (which stay constant for all realizations), and the random part. The fixed part goes into a "template", and the random part is sampled from the subroutine.

When `solve_sd()` is called, it uses the splitting position to separate the template into first stage and second stage. The template and the splitting point are passed into SD.

At each iteration or whenever SD needs to generate a new scenario, the user supplied subroutine is invoked to obtain the random part of the second stage LP coefficients. The template is then correspondingly modified in SD to form a scenario second stage problem.

### The LP Template

A template is a deterministic LP with the first stage and second stage combined, as shown below.

$$
\begin{align*}
\min \quad & c^\top x + \bar d^\top y \\
s.t. \quad & Ax = b, x \geq 0\\
& \bar Tx + Wy = \bar r, y\geq 0
\end{align*}
$$

When writing the template, only the fixed coefficient matters. The random coefficients in the second stage cost, the matrix $T$, and right-hand-side $r$ will be overwritten latter in the algorithm, so it can be taken the mean as a placeholder, or left out altogether.

The template is created using `direct_model` in JuMP, with CPLEX as the optimizer.

```julia
using TwoSD

model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

# OR for older CPLEX version
# model = direct_model(CPLEX.Optimizer())
```

The variables and constraints are added using the usual `@constraint` and `@variable`.

In addition to the usual requirements imposed by JuMP, we require the following conventions when defining the template:

1. Each constraint and variable must have a name, or belongs to a named container (e.g. array of variables/constraints). The name of the variable/constraint/container must not conflict with existing identifiers.
2. __Both stages are required to be Minimization problem.__ A max-max problem can be converted into min-min problem by negating the signs of the cost coefficients in both stages.
3. The first stage variables (resp. constraints) must be defined before the second stage variables (resp. constraints).

Note: It is also recommended that the constraints are written in the form of $a^\top x = b, a^\top x \leq b, a^\top x \geq b$. The variables are put on the left side of the equations, and the constant is on the right. This will make sure the signs will be correct when defining the random part.

### The Splitting Position

We need to let SD know the stageness of the constraints and variables listed in the template. Since we have listed all the first stage variables before the second stage variables, we can specify the first "second stage variable. We do the same for the constraints.

If the splitting constraint or variable happens to be inside a container, e.g. in an array, we specify the first element in that container.

We specify the splitting constraint and variable via the `Position` structure. The `Position` structure is a pair containing a reference to a constraint and a variable, and it can constructed using the following code.

```julia
split_position = Position(some_constraint, some_variable)
```

### Sampling Subroutine

The sampling subroutine should generate the random coefficients each time it is called. It should return an instance of `OneRealization`, which contains the instructions of which second stage LP coefficient to change (by specifying a `Position`), and what the value of that coefficient should be in that specific realization.

We begin by creating a function that generate some random numbers. These numbers should be used to construct a list of Pairs consisting of `Position` => values. The value must be a floating point value. The list is then used to construct the `OneRealization` and returned.

To signal a change in the right-hand-side, we put the string literal `"RHS"` as the variable when specifying the `Position`. If the change is in the second stage cost coefficient, put string literal `"OBJ"` as the constraint in the `Position`. The `Position` cannot be both "RHS" and "OBJ".

Note that the Position in the list should not change across scenarios. Also, if a bound on a second stage variable is random, a separate constraint should be written, since the bound constraint cannot be referenced if it is declared along with the variable.

```julia
function mystoc()
    # Generate some random numbers according to some distribution, or load a random entry from a dataset etc
    r1 = some_random_number()
    r2 = some_random_number()
    binding = [
        Position(second_stage_constraint, "RHS") => r1,
        Position(another_constraint, x1) => r2
    ]
    return OneRealization(binding)
end
```

__Important__: The random number generator should be seeded outside the function. Seeding the RNG inside would cause the function to return the same number.

Finally we need to specify the mean of the coefficients, __in the order of how the binding is constructed__. For the above example, that would be:

```julia
user_mean = [r1_mean, r2_mean]
```

## Running SD

After we gather the three parts above, we call SD by running:

```julia
result = solve_sd(model, split_position, user_mean, mystoc)
```

After the process completes, the output will be written into several files on the working directory.

In addition, the result will be populated with the incumbent solutions, and the compromise solution (if multiple replications are run, by default). To get the solutions, we call

```julia
# Get incumbent solution for one first stage variable x
decision(x, result, IncumbentSolution)

# Get compromise solution for x
decision(x, result, CompromiseSolution)

# Broadcasting the function to get incumbent solution for a container of variables s
decision.(s, Ref(result), CompromiseSolution)
```

Calling decision on second stage variables will cause an out-of-bound error. 

## Full example

Here is a demonstration of using TwoSD to solve a small Network Planning problem.

```julia
using CPLEX, JuMP
using TwoSD
using Distributions

# First we create a model to store the LP template
model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

# Stage 1 Variables
@variable(model, x[1:3] >= 0)

# Stage 2 Variables
@variables(model, begin
    s[1:3] >= 0
    f[1:3, 1:2] >= 0
end)

# Stage 1 Constraints
@constraint(model, budget, x[1] + x[2] + 4*x[3] <= 10)

# Stage 2 Constraints
# RHS randomness will be declared in the stoc part.
# RHS: omega[1:3] in link[1], link[2], link[3]
@constraint(model, link, s + f * [1, 1] .== 0)

# Random coefficient x1: -theta1
@constraint(model, capx1, -x[1] +f[1,1]+f[2,2]+f[3,2] <= 0)
# Random coefficient x2: -theta2
@constraint(model, capx2, -x[2] +f[1,2]+f[2,1]+f[3,2] <= 0)
@constraint(model, capx3, -x[3] +f[1,2]+f[2,2]+f[3,1] <= 0)

# The objective
@objective(model, Min, s[1] + s[2] + s[3])

# Specify where the first stage and second stage separates in the template.
split_position = Position(link[1], s[1])

# Demand distributions
dem_high = DiscreteNonParametric([0.0, 1, 2, 3, 5], [0.05, 0.2, 0.5, 0.2, 0.05])
dem_low = DiscreteNonParametric([0.0, 1, 2, 3], [0.1, 0.4, 0.4, 0.1])
link_up = Bernoulli(0.7)

function mystoc()
    return OneRealization(
        [
            # Demands are generated
            Position(link[1], "RHS") => rand(dem_high),
            Position(link[2], "RHS") => rand(dem_high),
            Position(link[3], "RHS") => rand(dem_low),

            # Is the link working? yes = 1 no = 0
            # Note negative sign in theta because it is on LHS.
            Position(capx1, x[1]) => -rand(link_up),
            Position(capx2, x[2]) => -rand(link_up)
        ]
    )
end

# By order of how the randomness are specified in mystoc(), we calculate the mean.
user_mean = [mean(dem_high), mean(dem_high), mean(dem_low), mean(link_up), mean(link_up)]

# Call SD
solution = solve_sd(model, split_position, user_mean, mystoc)

# Print solution if available
@show decision.(x, Ref(solution), AverageSolution)
@show decision.(x, Ref(solution), CompromiseSolution)
```

## Troubleshooting / Current Issues
    
1. Compilation issue: missing cplex header ilcplex.h
    
    Make sure you change Makefile to point to the correct CPLEX installation path. Make sure there is no space or non-ascii character in the path.

2. SD runs 5000 iterations without stopping and zero objective value.
    
    SD probably has received a feasibility LP. Make sure the objective is in the model.

3. SD crashes on start with ReadOnlyMemoryError.

    Please make sure the binding and split_position are correct. You should not try to modify the coefficients in the first stage.

4. SD crashes on start with exit code 1.

    Please form the problem in Min form.
