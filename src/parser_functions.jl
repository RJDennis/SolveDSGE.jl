function open_model_file(path::Q) where {Q <: AbstractString}

    #= Opens the model file and reads the contents, which are stored
       in a String-vector. =#

    model_file = open(path)
    model_array = readlines(model_file)
    close(model_file)

    # Remove lines or parts of lines that have been commented out.

    for i in eachindex(model_array)
        if occursin("#",model_array[i]) == true
            components = split(model_array[i],"#")
            keep = prod(components[isodd.(eachindex(components))])
            model_array[i] = keep
        end
    end

    # Remove blank lines

    model_array = model_array[model_array .!= ""]

    return model_array

end

function find_term(model_array::Array{Q,1}, term::Q) where {Q <: AbstractString}

    # Finds the position in the String-vector where model-terms (states, jumps,
    # parameters, equations, etc) are located.

    locations = findall(y->y == term, model_array)
    if length(locations) == 1
        return locations[1]
    elseif length(locations) == 0
        error("The $term-designation does not appear in the model file")
    else
        error("The $term-designation appears multiple times in the model file")
    end

end

function find_end(model_array::Array{Q,1}, startfrom::S) where {Q <: AbstractString,S <: Integer}

    # Each of the key model terms must be followed by an 'end', whose location
    # in the String-array this function finds.
    # White space at the start of a line is stripped out.

    for end_location = startfrom:length(model_array)
        if startswith(strip(model_array[end_location]), "end") == true
            return end_location
        end
    end

    error("Model file is missing an 'end' term")

end

function get_variables(model_array::Array{Q,1}, term::Q) where {Q <: AbstractString}

    #= This function extracts the variable names and ensures that no names
       are repeated =#

    term_begin = find_term(model_array, term) + 1
    term_end   = find_end(model_array, term_begin) - 1

    if term_begin > term_end
        if term == "shocks:"
            return [""]
        else
            error("The model file contains no $(term[1:end-1])")
        end
    end

    term_block = model_array[term_begin:term_end]

    # Extract the variable names

    # First remove any trailing variable separators: "," or ";".

    for i = 1:length(term_block)
        if endswith(term_block[i],union(",",";"))
            term_block[i] = term_block[i][1:end-1]
        end
    end

    variables = String.(strip.(split(term_block[1], union(",",";"))))
    for i = 2:length(term_block)
        variables = [variables;String.(strip.(split(term_block[i], union(",",";"))))]
    end

    if length(variables) != length(unique(variables))
        error("Some $(term[1:end-1]) are repreated")
    end

    return variables

end

function count_variables(x::Array{Q,1}) where {Q <: AbstractString}

    # Counts the number of variables.

    if x[1] .== ""
        return 0
    else
        return length(x)
    end

end

function combine_states_and_jumps(x::Array{Q,1}, y::Array{Q,1}) where {Q <: AbstractString}

    # This function combines the states, "x", with the jump variables, "y", to
    # genetate the model's variables.

    if x[1] .== ""
        return y
    else
        variables = union(x, y)
        if length(variables) != length(x) + length(y)
            error("Some states and jumps have the same name.")
        end
        return variables
    end

end

function create_lag_variables(x::Array{Q,1}) where {Q <: AbstractString}

    # Creates an array containing the names of all possible lagged variables.
    # None, or only a small set, of these lagged variables may appear in any
    # given model.

    # x is the model variables (states + jumps)

    lag_variables = copy(x)
    for i = 1:length(lag_variables)
        lag_variables[i] = string(x[i],"(-1)")
    end

    return lag_variables

end

function create_lead_variables(x::Array{Q,1}) where {Q <: AbstractString}

    # Creates an array containing the names of all possible future variables.
    # Not all of these future variables may appear in any given model.

    # x is the model variables (states + jumps)

    lead_variables = copy(x)
    for i = 1:length(lead_variables)
        lead_variables[i] = string(x[i],"(+1)")
    end

    return lead_variables

end

function get_parameters_and_values(model_array::Array{Q,1}, term::Q) where {Q <: AbstractString}

    # This function extracts the names and associated values for each of the
    # model's parameters.  The parameter names are sorted so that larger names
    # come first.

    parametersbegin = find_term(model_array, term) + 1
    parametersend   = find_end(model_array, parametersbegin) - 1
    if parametersbegin > parametersend
        error("The model file contains no $(term[1:end-1])")
    end

    parameterblock  = model_array[parametersbegin:parametersend]

    # Extract the parameter names and values

    # First remove any trailing separators: "," or ";".

    for i = 1:length(parameterblock)
        if endswith(parameterblock[i],union(",",";"))
            parameterblock[i] = parameterblock[i][1:end-1]
        end
    end

    revised_parameterblock = String.(strip.(split(parameterblock[1], union(",",";"))))
    for i = 2:length(parameterblock)
        revised_parameterblock = [revised_parameterblock;String.(strip.(split(parameterblock[i], union(",",";"))))]
    end

    parameters = Array{Q}(undef,length(revised_parameterblock))
    values     = Array{Float64}(undef,length(revised_parameterblock))
    for i = 1:length(revised_parameterblock)
        if occursin("=",revised_parameterblock[i]) == false
            error("Parameter line $i does not contain an '=' sign")
        else
            pair = strip.(split(revised_parameterblock[i],"="))
            parameters[i] = pair[1]
            values[i] = Meta.parse(pair[2])
        end
    end

    parameter_order   = sortperm(length.(parameters),rev = true)
    sorted_parameters = parameters[parameter_order]
    sorted_values     = values[parameter_order]

    return sorted_parameters, sorted_values

end

function get_equations(model_array::Array{Q,1}, term::Q) where {Q <: AbstractString}

    # Extract the model's equations.

    equationsbegin  = find_term(model_array, term) + 1
    equationsend    = find_end(model_array, equationsbegin) - 1
    if equationsbegin > equationsend
        error("The model file contains no $(term[1:end-1])")
    end

    equation_block  = model_array[equationsbegin:equationsend]

    # Extract the equations

    # First remove any trailing separators: "," or ";".

    for i = 1:length(equation_block)
        if endswith(equation_block[i],union(",",";"))
            equation_block[i] = equation_block[i][1:end-1]
        end
    end

    revised_equations = String.(strip.(split(equation_block[1], union(",",";"))))
    for i = 2:length(equation_block)
        revised_equations = [revised_equations;String.(strip.(split(equation_block[i], union(",",";"))))]
    end

    for i = 1:length(equation_block)
        if occursin("=",revised_equations[i]) == false
            error("Equation line $i does not contain an '=' sign")
        end
    end

    return revised_equations

end

function get_re_model_primatives(model_array::Array{Q,1}) where {Q <: AbstractString}

    # This function takes the model-array read from a model file, extracts the
    # critical model information, does some basic error checking, and returns
    # it in a structure.

    states                        = get_variables(model_array, "states:")
    jumps                         = get_variables(model_array, "jumps:")
    shocks                        = get_variables(model_array, "shocks:")
    variables                     = combine_states_and_jumps(states, jumps)
    lag_variables                 = create_lag_variables(variables)
    lead_variables                = create_lead_variables(variables)
    equations                     = get_equations(model_array, "equations:")
    (parameters, parametervalues) = get_parameters_and_values(model_array, "parameters:")

    for i in [variables;parameters]
        if sum(occursin.(i,equations)) == false
            println("Warning: $i is not in any equation")
        end
    end

    combined_names = [parameters;variables;shocks]
    if length(unique(combined_names)) != length(combined_names)
        error("Some parameters, variables, or shocks have the same name")
    end

    if sum(combined_names .== "exp") != 0
        error("'exp' cannot be a name for a variable, a shock, or a parameter.")
    end

    re_model_primatives = REModelPrimatives(states,jumps,shocks,variables,lag_variables,lead_variables,parameters,parametervalues,equations)

    return re_model_primatives

end

function reorder_equations(model::DSGEModel)

    # This function reorders the model's equations so that the equations
    # containing the shock processes appear first.

    equations           = model.equations
    reordered_equations = copy(equations)

    shocks = model.shocks
    states = model.states
    jumps  = model.jumps
    reordered_states = copy(states)

    states_that_have_been_reordered = Int64[]

    # Construct summary information about each equation

    shocks_number = zeros(Int64,length(equations))
    states_number = zeros(Int64,length(equations))
    jumps_number  = zeros(Int64,length(equations))
    for i = 1:length(equations)
        if length(shocks) != 0
            shocks_number[i] = sum(occursin.(shocks,equations[i]))
        end
        states_number[i] = sum(occursin.(states,equations[i]))
        jumps_number[i] = sum(occursin.(jumps,equations[i]))
    end
    number_eqns_with_shocks = sum(shocks_number .!= 0)

    # Put the equations with no jumps in them at the top

    pos = 1
    for i = 1:length(equations)
        if jumps_number[i] == 0 && i != pos
            reordered_equations[pos], reordered_equations[i] = reordered_equations[i], reordered_equations[pos]
            shocks_number[pos], shocks_number[i] = shocks_number[i], shocks_number[pos]
            states_number[pos], states_number[i] = states_number[i], states_number[pos]
            jumps_number[pos], jumps_number[i]   = jumps_number[i], jumps_number[pos]
            pos += 1
        end
    end

    # Put the shock equations with the fewest states at the top

    reordered_equations[1:number_eqns_with_shocks] .= reordered_equations[sortperm(states_number[1:number_eqns_with_shocks])]
    shocks_number[1:number_eqns_with_shocks] .= shocks_number[sortperm(states_number[1:number_eqns_with_shocks])]
    states_number[1:number_eqns_with_shocks] .= states_number[sortperm(states_number[1:number_eqns_with_shocks])]
    jumps_number[1:number_eqns_with_shocks]  .= jumps_number[sortperm(states_number[1:number_eqns_with_shocks])]

    # Now sort out the order of the states in the system

    pos = 1
    for i = 1:number_eqns_with_shocks
        for j = pos:length(states)
            if occursin(states[j],reordered_equations[i]) == true && sum(j.==states_that_have_been_reordered) == 0
                reordered_states[pos], reordered_states[j] = reordered_states[j], reordered_states[pos]
                states_that_have_been_reordered = [states_that_have_been_reordered;j;]
                pos += 1
            end
        end
    end

    return reordered_equations, reordered_states

end

function reorganize_equations(equations::Array{Q,1},model::DSGEModel) where {Q <: AbstractString}

    # This function replaces lagged variables with pseudo current variables,
    # augmenting the state vector and the model's equations accordingly.

    states          = model.states
    jumps           = model.jumps
    variables       = model.variables
    lag_variables   = model.lag_variables

    reorganized_equations = copy(equations)
    reorganized_states    = copy(states)

    # First we want to determine if any equation contains a lagged variable.

    model_has_lags = false
    for i = 1:length(equations)
        if sum(occursin.(lag_variables,equations[i])) != 0
            model_has_lags = true
            break
        end
    end

    # If a model contains lagged variables, then we introduce a pseudo variable
    # in its place and augment the list of state variables and the set of model
    # equations.

    if model_has_lags == true
        for j = 1:length(lag_variables)
            if sum(occursin.(lag_variables[j],equations)) != 0
                for i = 1:length(equations)
                    reorganized_equations[i] = replace(reorganized_equations[i], lag_variables[j] => string(variables[j],"lag"))
                end
                reorganized_states    = [reorganized_states;string(variables[j],"lag")]
                reorganized_equations = [reorganized_equations;string(variables[j],"lag(+1) = ",variables[j])]
            end
        end
    end

    reorganized_variables = union(reorganized_states,jumps)

    return reorganized_equations, reorganized_states, reorganized_variables

end

function repackage_equations(equations::Array{Q,1},variables::Array{Q,1},model::DSGEModel) where {Q <: AbstractString}

    # This function is critical for repackaging the model's equations, replacing
    # parameter names with values, and numbering variables.

    shocks          = model.shocks
    lead_variables  = model.lead_variables
    parameters      = model.parameters
    parametervalues = model.parametervalues

    repackaged_equations = copy(equations)

    if count_variables(shocks) != 0
        combined_names = [variables;parameters;shocks]
    else
        combined_names = [variables;parameters]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]

    # Now we go through every equation and replace future variables, variables, and
    # shocks with a numbered element of a vector, "x".  We also replace parameter
    # names with parameter values.

    for j in sorted_combined_names
        if sum(j .== variables) == 1
            variable_index = findfirst(isequal(j),variables)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],"$j(+1)" => "x[$(length(variables) + variable_index)]")
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(variable_index)]")
            end
        elseif sum(j .== parameters) == 1
            parameter_index = findfirst(isequal(j),parameters)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => parametervalues[parameter_index])
            end
        elseif count_variables(shocks) != 0 && sum(j .== shocks) == 1
            shock_index = findfirst(isequal(j),shocks)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(length(variables) + length(lead_variables) + shock_index)]")
            end
        end
    end

    return repackaged_equations

end

function create_steady_state_equations(equations::Array{Q,1}, model::DSGEModel) where {Q <: AbstractString}

    # Make the model static by replacing leads and lags with current variables
    # and setting shocks equal to zero.  Also, replace parameter names with
    # their associated value.

    steady_state_equations = copy(equations)

    variables       = model.variables
    shocks          = model.shocks
    parameters      = model.parameters
    parametervalues = model.parametervalues

    if count_variables(shocks) != 0
        combined_names = [variables;parameters;shocks]
    else
        combined_names = [variables;parameters]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]

    # Now we go through every equation and replace future variables, variables, and
    # shocks with a numbered element of a vector, "x".  We also replace parameter
    # names with parameter values

    for j in sorted_combined_names
        if j in variables
            variable_index = findfirst(isequal(j),variables)
            for i = 1:length(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],"$j(-1)" => "x[$(variable_index)]")
                steady_state_equations[i] = replace(steady_state_equations[i],"$j(+1)" => "x[$(variable_index)]")
                steady_state_equations[i] = replace(steady_state_equations[i],j => "x[$(variable_index)]")
            end
        elseif j in parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i = 1:length(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => parametervalues[parameter_index])
            end
        elseif count_variables(shocks) != 0 && j in shocks
            shock_index = findfirst(isequal(j),shocks)
            for i = 1:length(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => 0.0)
            end
        end
    end

    return steady_state_equations

end

function make_equations_equal_zero(equations::Array{Q,1}) where {Q <: AbstractString}

    # Reexpress all of the model's equations such that they equal zero.

    zeroed_equations = similar(equations)

    for i = 1:length(equations)
        pair = strip.(split(equations[i], "="))
        zeroed_equations[i] = string(pair[1], " - (", pair[2], ")")
    end

    return zeroed_equations

end

function create_projection_equations(equations::Array{Q,1},model::DSGEModel) where {Q <: AbstractString}

    projection_equations = copy(equations)

    nx = length(model.states)
    ny = length(model.jumps)
    ns = length(model.shocks)
    ne = length(equations)
    nv = nx + ny

    for j = 1:nx
        for i = 1:ne

            projection_equations[i] = replace(projection_equations[i], "x[$j]" => "state[$j]")

        end
    end

    for j = 1:nv
        for i = 1:ne

            projection_equations[i] = replace(projection_equations[i], "x[$(nx+j)]" => "x[$j]")

        end
    end

    for j = 1:ny
        for i = 1:ne

            projection_equations[i] = replace(projection_equations[i], "x[$(nx+nv+j)]" => "approx$j")

        end
    end

    for j = 1:ns
        for i = 1:ne

            projection_equations[i] = replace(projection_equations[i], "x[$(2*nv+j)]" => 0.0)

        end
    end

    jumps_to_be_approximated = Int64[]
    for i = 1:ny
        for j = 1:ne
            if occursin("approx$i",projection_equations[j]) == true
                push!(jumps_to_be_approximated,i)
            end
        end
    end

    jumps_to_be_approximated = unique(jumps_to_be_approximated)

    return projection_equations, jumps_to_be_approximated

end

function re_model_processed(re_model_primatives::REModelPrimatives{Q,T}) where {Q <: AbstractString,T <: AbstractFloat}

    # This function processes the primative information about the model to put
    # it in a form that the model-solvers can operate on.

    reordered_equations, states = reorder_equations(re_model_primatives)
    re_model_primatives.states  = states

    reorganized_equations, states, variables = reorganize_equations(reordered_equations,re_model_primatives)

    re_model_primatives.states         = states
    re_model_primatives.equations      = reorganized_equations
    re_model_primatives.variables      = variables
    re_model_primatives.lag_variables  = create_lag_variables(variables)
    re_model_primatives.lead_variables = create_lead_variables(variables)

    repackaged_equations   = repackage_equations(reorganized_equations,variables,re_model_primatives)
    steady_state_equations = create_steady_state_equations(reorganized_equations,re_model_primatives)

    static_function  = function faux_static_function() end
    dynamic_function = function faux_dynamic_function() end

    each_eqn_function = Array{Function}(undef, length(reorganized_equations))
    each_eqn_function .= function faux_eqn_function() end

    projection_equations, jumps_to_be_approximated = create_projection_equations(repackaged_equations,re_model_primatives)
    closure_function = function faux_closure_function() end
    closure_function_piecewise = function faux_closure_function() end

    number_states    = count_variables(states)
    number_jumps     = count_variables(re_model_primatives.jumps)
    number_shocks    = count_variables(re_model_primatives.shocks)
    number_variables = count_variables(variables)
    number_equations = count_variables(reorganized_equations)

    dsge = REModel(number_states, number_jumps, number_shocks, number_variables, variables, number_equations, steady_state_equations, repackaged_equations, faux_static_function, faux_static_function, faux_dynamic_function, each_eqn_function, projection_equations, closure_function, closure_function_piecewise, jumps_to_be_approximated)

    return dsge

end

function create_model_functions(model::REModel, path::AbstractString)

    # Takes the model's equations (which are read in as strings) and turns these
    # equations into functions.  There are functions for the static equations,
    # functions for the dynamic equations, and functions for each individual
    # dynamic equation.  These functions are inserted into the processed model
    # structure and they are saved as text files in the same folder as the model
    # file.

    static_equations = make_equations_equal_zero(model.steady_state_equations)

    static_string         = "function static_equations(x::Array{T,1}) where {T<:Number} \n \n"
    static_string         = string(static_string, "  f = Array{T}(undef,length(x)) \n \n")
    nlsolve_static_string = "function nlsolve_static_equations(f,x) \n \n"
    for i = 1:length(static_equations)
        static_string         = string(static_string, "  f[$i] = ", static_equations[i], "\n")
        nlsolve_static_string = string(nlsolve_static_string, "  f[$i] = ", static_equations[i], "\n")
    end

    static_string         = string(static_string, "\n  return f \n \n", "end")
    nlsolve_static_string = string(nlsolve_static_string, "\n", "end")

    static_model_path = replace(path, ".txt" => "_static.txt")

    model.static_function         = eval(Meta.parse(static_string))
    model.nlsolve_static_function = eval(Meta.parse(nlsolve_static_string))
    open(static_model_path, "w") do io write(io, string(nlsolve_static_string, " \n \n", static_string)) end

    ###

    dynamic_equations = make_equations_equal_zero(model.repackaged_equations)
    n_eqns = length(dynamic_equations)
    dynamic_string = "function dynamic_equations(x::Array{T,1}) where {T<:Number} \n \n"
    dynamic_string = string(dynamic_string, "  f = Array{T}(undef,$n_eqns) \n \n")
    for i = 1:length(dynamic_equations)
        dynamic_string = string(dynamic_string, "  f[$i] = ", dynamic_equations[i], "\n")
    end

    dynamic_string = string(dynamic_string, "\n  return f \n \n end \n")

    model.dynamic_function = eval(Meta.parse(dynamic_string))

    each_equation_string = Array{String}(undef, length(dynamic_equations))

    for i = 1:length(dynamic_equations)
        each_equation_string[i] = "function dynamic_eqn_$i(x::Array{T,1}) where {T<:Number} \n \n"
        each_equation_string[i] = string(each_equation_string[i], "  f = ", dynamic_equations[i], "\n", "\n  return f \n \n", "end \n")
        model.each_eqn_function[i] = eval(Meta.parse(each_equation_string[i]))
    end

    for i = 1:length(dynamic_equations)
        dynamic_string = string(dynamic_string, "\n", each_equation_string[i])
    end

    projection_equations = make_equations_equal_zero(model.nonlinear_equations)
    ny = model.number_jumps

    # For Chebyshev and Smolyak

    closure_string = "function closure_projection_equations(state,scaled_weights,order,domain,approximate) \n \n"
    closure_string = string(closure_string, "  function projection_equations(x::Array{T,1}) where {T<:Number} \n \n")
    weight_number = 1
    for i in model.jumps_approximated
        closure_string = string(closure_string, "    approx$i = approximate(scaled_weights[$weight_number],x[$ny+1:end],order,domain)", "\n")
        weight_number += 1
    end
    closure_string = string(closure_string, "\n", "    f = Array{T}(undef,$n_eqns) \n \n")
    for i = 1:length(projection_equations)
        closure_string = string(closure_string, "    f[$i] = ", projection_equations[i], "\n")
    end
    closure_string = string(closure_string, "\n    return f \n \n  end \n \n  return projection_equations \n \n", "end \n")

    model.closure_function = eval(Meta.parse(closure_string))

    # For piecewise linear

    ns = model.number_shocks # We need to separate the function generated for the stochastic and deterministic cases
    if ns == 0
        closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,approximate) \n \n"
    else
        closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,integrals,approximate) \n \n"
    end
    closure_pl_string = string(closure_pl_string, "  function projection_equations_pl(x::Array{T,1}) where {T<:Number} \n \n")
    for i in model.jumps_approximated
        if ns == 0
            closure_pl_string = string(closure_pl_string, "    approx$i = approximate(variables[$i],grid,x[$ny+1:end])", "\n")
        else
            closure_pl_string = string(closure_pl_string, "    approx$i = approximate(variables[$i],grid,x[$ny+1:end],integrals)", "\n")
        end
    end
    closure_pl_string = string(closure_pl_string, "\n", "    f = Array{T}(undef,$n_eqns) \n \n")
    for i = 1:length(projection_equations)
        closure_pl_string = string(closure_pl_string, "    f[$i] = ", projection_equations[i], "\n")
    end
    closure_pl_string = string(closure_pl_string, "\n    return f \n \n  end \n \n  return projection_equations_pl \n \n", "end \n")

    model.closure_function_piecewise = eval(Meta.parse(closure_pl_string))

    dynamic_string = string(dynamic_string, "\n", closure_string)
    dynamic_string = string(dynamic_string, "\n", closure_pl_string)

    dynamic_model_path = replace(path, ".txt" => "_dynamic.txt")
    open(dynamic_model_path, "w") do io write(io, dynamic_string) end

end

function get_re_model(model_array::Array{Q,1}, path::Q) where {Q <: AbstractString}

    # Creates the processed model structure for rational expectations models
    # (anticipating that other types of models may come later).

    re_model_primatives = get_re_model_primatives(model_array)
    re_model            = re_model_processed(re_model_primatives)

    create_model_functions(re_model, path)

    return re_model

end

function get_model(path::Q) where {Q <: AbstractString}

    # Main function used to open, read, and process a model file.  The structure
    # returned from this function presents the model in the form required for
    # the model-solvers.

    model_array = open_model_file(path)

    model = get_re_model(model_array, path)

    println("The model's variables are now in this order: ", model.variables)

    return model

end
