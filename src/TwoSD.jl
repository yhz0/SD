module TwoSD

# Includes
include("SDTypes.jl")
include("SDPosition.jl")
include("SDAPI.jl")

## Imports
using CPLEX, MathOptInterface, JuMP

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using .SDPosition
using .SDAPI

# SD Structure

function build_constraint_mapping(model::Model)::Dict{ConstraintRef, Int}
    constr = Vector{ConstraintRef}()
    for (f, s) in list_of_constraint_types(model)
        if f == VariableRef
            continue
        end

        for c in all_constraints(model, f, s)
            push!(constr, c)
        end
    end
    sort!(constr, by= con -> optimizer_index(con).value)
    return Dict(zip(constr, range(1, length(constr))))
end

function build_variable_mapping(model)::Dict{VariableRef, Int}
    v = all_variables(model)
    sort!(v, by= v -> optimizer_index(v).value)
    return Dict(zip(v, range(1, length(v))))
end

function filler_template(observ_p::Ptr{Cdouble}, gen::Function, pattern_dict)
    or::OneRealization = gen()
    observ = unsafe_wrap(Array, observ_p, length(pattern_dict))
    for (pos, val) in or.value
        ind = pattern_dict[pos]
        observ[ind] = convert(Cdouble, val)
    end
end

function initialize_sd()
    SDAPI.openSolver()
    SDAPI.readConfig_jl()
end

# Parse a comma delimited string as numbers.
# There might be some precision loss during the conversion.
function parse_solution(s::String)
    v = split(s, ",")
    return parse.(Float64, v)
end

# Load all solutions from a file, one solution on each line
function load_solution_file(path::String)
    solutions = Vector{Vector{Float64}}()
    if !isfile(path)
        @warn "load_solution_file: $(path) does not exist."
        return solutions
    end
    open(path, "r") do io
        for line in eachline(io)
            push!(solutions, parse_solution(line))
        end
    end
    return solutions
end

# Structure for storing the sd output
mutable struct SDResult
    row_map # constraint mapping from jump to index
    col_map # variable mapping from jump to index

    row_map_r # mapping from index to constraints
    col_map_r # mapping from index to variables

    incumbentX # incumbent solution
    compromiseX # compromise solution
    avgX # average solution
end

# Reverse key/value pair in a dict
function reverse_dict(my_dict::Dict)
    return Dict(value => key for (key, value) in my_dict)
end

# Load SD output files from current directory.
function build_result(row_map=nothing, col_map=nothing, basepath = ".")
    incumbX = load_solution_file(joinpath(basepath, "spOutputincumb.dat"))
    compromiseX = load_solution_file(joinpath(basepath, "spOutputcompromiseX.dat"))
    avgX = load_solution_file(joinpath(basepath, "spOutputavgX.dat"))

    row_map_r = row_map === nothing ? nothing : reverse_dict(row_map)
    col_map_r = col_map === nothing ? nothing : reverse_dict(row_map)

    return SDResult(row_map, col_map, row_map_r, col_map_r, incumbX, compromiseX, avgX)
end

@enum SDSolutionType IncumbentSolution=1 CompromiseSolution=2 AverageSolution=3

# Get the decision from a given result, with selected type (Incumbent, Compromise, Average)
function decision(var::VariableRef, result::SDResult, type::SDSolutionType=Incumbent)
    if type == IncumbentSolution
        v = result.incumbentX[1]
    elseif type == AverageSolution
        v = result.avgX[1]
    elseif type == CompromiseSolution
        v = result.compromiseX[1]
    end
    index = result.col_map[var]
    return v[index]
end

# Check if stoc_row and stoc_col are in the supported place,
# relatively to split_index. The randomness should only appear in
# T matrix, d second stage coefficients, r second stage RHS.
# If not throw out a warning.
function check_stoc_pos_validity(stoc_row::Int, stoc_col::Int, split_index_row::Int, split_index_col::Int)
    if stoc_row == -1 && stoc_col >= split_index_col
        # Second stage objective
        return "RANDOM_OBJ"
    elseif stoc_col == -1 && stoc_row >= split_index_row
        # Second stage RHS
        return "RANDOM_RHS"
    elseif stoc_row >= split_index_row && stoc_col < split_index_col
        # T Matrix
        return "RANDOM_T"
    else
        @warn "Unsupported Random Position at ($stoc_row, $stoc_col). Check stageness."
    end
end

function solve_sd(model::Model, split_position::Position,
    user_mean::AbstractVector{<:AbstractFloat}, mystoc::Function)
    # This is what build should do
    cpx_model = backend(model)
    initialize_sd()

    # Build COR
    cor = SDAPI.populateCore(cpx_model.lp, "JuliaExportProblem")

    # Use these to define the position translator
    row_map = build_constraint_mapping(model)
    col_map = build_variable_mapping(model)
    index(p::Position) = SDPosition.pos_to_index(p, row_map, col_map)

    # Build TIM
    split_index = index(split_position)
    tim = SDAPI.buildTwoStageTime(split_index.row, split_index.col)

    # Build Stoc
    sample_realization = mystoc()
    stoc_numOmega = length(sample_realization.value)
    pattern = SDPosition.get_pattern(sample_realization)
    pattern_dict = SDPosition.bind_pattern(pattern) # save this for later index lookups

    stoc_row = [p.row for p in index.(pattern)]
    stoc_col = [p.col for p in index.(pattern)]
    stoc_mean = convert(Vector{Float64}, user_mean) # have the user define this

    # Validity check: to ensure the randomness only appears in
    # T matrix, d second stage coefficients, r second stage RHS
    check_stoc_pos_validity.(stoc_row, stoc_col, split_index.row, split_index.col)

    # Make the user generation function to a c-func pointer
    mystoc_filler(observ_p::Ptr{Cdouble}) = filler_template(observ_p, mystoc, pattern_dict)
    mystoc_filler_p = @cfunction($mystoc_filler, Cvoid, (Ptr{Cdouble}, )) 

    stoc = SDAPI.buildExtStocType(stoc_numOmega, stoc_row, stoc_col, stoc_mean, mystoc_filler_p)

    SDAPI.dumpCore(cor) # print
    SDAPI.dumpTime(tim) # print
    SDAPI.dumpStoc(stoc)

    SDAPI.algo(cor, tim, stoc, ".", "JuliaExportProblem")

    result = build_result(row_map, col_map)
    return result
end

export Position, OneRealization, solve_sd, decision
export IncumbentSolution, CompromiseSolution, AverageSolution
end