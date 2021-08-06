################# Parser functions ######################

function open_model_file(path::Q) where {Q <: AbstractString}

    #= Opens the model file and reads the contents, which are stored
       in a vector of strings. =#

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

function find_term(model_array::Array{Q,1},term::Q) where {Q <: AbstractString}

    #= Finds the position in the String-vector where model-terms (states, jumps,
       parameters, equations, etc) are located. =#

    locations = findall(y->y == term, model_array)
    if length(locations) == 1
        return locations[1]
    elseif length(locations) == 0
        error("The $term-designation does not appear in the model file.")
    else
        error("The $term-designation appears multiple times in the model file.")
    end

end

function find_end(model_array::Array{Q,1},startfrom::S) where {Q <: AbstractString,S <: Integer}

    #= Each of the key model terms must be followed by an 'end', whose location
       in the vector this function finds.  White space at the start of a line is
       stripped out. =#

    for end_location = startfrom:length(model_array)
        if startswith(strip(model_array[end_location]), "end") == true
            return end_location
        end
    end

    error("Model file is missing an 'end' term.")

end

function get_variables(model_array::Array{Q,1},term::Q) where {Q <: AbstractString}

    #= This function extracts the variable names and ensures that no names
       are repeated =#

    term_begin = find_term(model_array, term) + 1
    term_end   = find_end(model_array, term_begin) - 1

    if term_begin > term_end
        if term in ["shocks:","states:"]
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    end

    term_block = model_array[term_begin:term_end]

    # Extract the variable names

    # Remove any trailing variable separators: "," or ";".

    for i = 1:length(term_block)
        if endswith(term_block[i],union(",",";"))
            term_block[i] = term_block[i][1:end-1]
        end
    end

    # Extract the names and place them in a vector

    variables = String.(strip.(split(term_block[1], union(",",";"))))
    for i = 2:length(term_block)
        variables = [variables;String.(strip.(split(term_block[i], union(",",";"))))]
    end

    # Check whether names are repeated

    if length(variables) != length(unique(variables))
        error("Some $(term[1:end-1]) are repreated.")
    end

    # Check to ensure that variables contains non-empty string elements

    if length(variables) == 1 && variables[1] == ""
        if term in ["shocks:","states:"]
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    else
        return variables
    end

end

function combine_states_and_jumps(x::Array{Q,1},y::Array{Q,1}) where {Q <: AbstractString}

    #= This function combines the states, "x", with the jump variables, "y", to
       generate the model's variables. =#

    if length(x) == 0 # There are no states
        return y
    elseif length(intersect(x,y)) > 0
        error("Some states and jumps have the same name.")
    else
        return [x;y]
    end

end

function get_parameters_and_values(model_array::Array{Q,1},term::Q) where {Q <: AbstractString}

    #= This function extracts the names and associated values for each of the
       model's parameters.  The parameter names are sorted so that larger names
       come first. =#

    parametersbegin = find_term(model_array, term) + 1
    parametersend   = find_end(model_array, parametersbegin) - 1
    if parametersbegin > parametersend
        error("The model file contains no $(term[1:end-1])")
    end

    parameterblock = model_array[parametersbegin:parametersend]

    # Extract the parameter names and values

    # Remove any trailing separators: "," or ";".

    for i = 1:length(parameterblock)
        if endswith(parameterblock[i],union(",",";"))
            parameterblock[i] = parameterblock[i][1:end-1]
        end
    end

    revised_parameterblock = String.(strip.(split(parameterblock[1], union(",",";"))))
    for i = 2:length(parameterblock)
        revised_parameterblock = [revised_parameterblock;String.(strip.(split(parameterblock[i], union(",",";"))))]
    end

    # Extract the parameter names and values

    unassigned_parameter_index = 1
    unassigned_parameters = Array{Q}(undef,0)

    parameters = Array{Q}(undef,length(revised_parameterblock))
    values     = Array{Q}(undef,length(revised_parameterblock))
    for i = 1:length(revised_parameterblock)
        if occursin("=",revised_parameterblock[i]) == false
            parameters[i] = revised_parameterblock[i]
            values[i] = "p[$unassigned_parameter_index]" # p is a reserved name
            push!(unassigned_parameters,revised_parameterblock[i])
            unassigned_parameter_index += 1
        else
            pair = strip.(split(revised_parameterblock[i],"="))
            parameters[i] = pair[1]
            values[i] = pair[2]
        end
    end

    parameter_order   = sortperm(length.(parameters),rev = true)
    sorted_parameters = parameters[parameter_order]
    sorted_values     = values[parameter_order]

    return sorted_parameters, sorted_values, unassigned_parameters

end

function get_equations(model_array::Array{Q,1},term::Q) where {Q <: AbstractString}

    # Extract the model's equations.

    equationsbegin  = find_term(model_array, term) + 1
    equationsend    = find_end(model_array, equationsbegin) - 1
    if equationsbegin > equationsend
        error("The model file contains no $(term[1:end-1])")
    end

    equation_block  = model_array[equationsbegin:equationsend]

    # Extract the equations

    # Remove any trailing separators: "," or ";".

    for i = 1:length(equation_block)
        if endswith(equation_block[i],union(",",";"))
            equation_block[i] = equation_block[i][1:end-1]
        end
    end

    # Extract the equations and place them in a vector

    equations = String.(strip.(split(equation_block[1], union(",",";"))))
    for i = 2:length(equation_block)
        equations = [equations;String.(strip.(split(equation_block[i], union(",",";"))))]
    end

    for i = 1:length(equation_block)
        if occursin("=",equations[i]) == false
            error("Equation line $i does not contain an '=' sign.")
        end
    end

    return equations

end

function reorder_equations(equations::Array{Q,1},shocks::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1}) where {Q <: AbstractString}

    #= This function reorders the model's equations so that the equations
       containing the shock processes appear first. =#

    reordered_equations = copy(equations)
    reordered_states = copy(states)

    # Construct summary information about each equation

    shocks_number = zeros(Int64,length(equations))
    states_number = zeros(Int64,length(equations))
    jumps_number  = zeros(Int64,length(equations))
    for i = 1:length(equations)
        if length(shocks) != 0
            shocks_number[i] = sum(occursin.(shocks,equations[i]))
        end
        if length(states) != 0
            states_number[i] = sum(occursin.(states,equations[i]))
        end
        jumps_number[i] = sum(occursin.(jumps,equations[i]))
    end
    number_eqns_with_shocks = sum(shocks_number .!= 0)

    # Put the equations with no jumps in them at the top

    pos = 1
    for i = 1:length(equations)
        if jumps_number[i] == 0 
            if i != pos
                reordered_equations[pos], reordered_equations[i] = reordered_equations[i], reordered_equations[pos]
                shocks_number[pos], shocks_number[i] = shocks_number[i], shocks_number[pos]
                states_number[pos], states_number[i] = states_number[i], states_number[pos]
                jumps_number[pos], jumps_number[i]   = jumps_number[i], jumps_number[pos]
            end
            pos += 1
        end
    end

    # Put the shock equations with the fewest shocks at the top

    new_order = sortperm(shocks_number[1:number_eqns_with_shocks])
    reordered_equations[1:number_eqns_with_shocks] .= reordered_equations[new_order]
    shocks_number[1:number_eqns_with_shocks] .= shocks_number[new_order]
    states_number[1:number_eqns_with_shocks] .= states_number[new_order]
    jumps_number[1:number_eqns_with_shocks]  .= jumps_number[new_order]
    
    # Order the shock processes to try to make shock's variance matrix have non-zero diagonals.  Only applies to stochastic models.

    for j = 1:number_eqns_with_shocks
        if shocks_number[j] == 1
           for i = 1:length(shocks)
                if occursin(shocks[i],reordered_equations[j]) && j != i
                    reordered_equations[i], reordered_equations[j] = reordered_equations[j], reordered_equations[i]
                    shocks_number[j], shocks_number[i] = shocks_number[i], shocks_number[j]
                    states_number[j], states_number[i] = states_number[i], states_number[j]
                    jumps_number[j], jumps_number[i]   = jumps_number[i], jumps_number[j]
                    break
                end
            end
        end
    end

    #pos = 1
    #for i = 1:length(shocks)
    #    for j = pos:number_eqns_with_shocks
    #        if occursin(shocks[i],reordered_equations[j]) && j != pos
    #            reordered_equations[pos], reordered_equations[j] = reordered_equations[j], reordered_equations[pos]
    #            break
    #        end
    #    end
    #    pos += 1
    #end

    # Now sort out the order of the states in the system

    states_that_have_been_reordered = Int64[]
    pos = 1
    if number_eqns_with_shocks != 0 # Get the right ordering for states for stochastic models
        for i = 1:number_eqns_with_shocks
            for j = pos:length(states)
                if occursin(states[j],reordered_equations[i]) == true && (j in states_that_have_been_reordered) == false
                    reordered_states[pos], reordered_states[j] = reordered_states[j], reordered_states[pos]
                    push!(states_that_have_been_reordered,j)
                    pos += 1
                end
            end
        end
    else
        for i = 1:length(equations)
            if states_number[i] != 0 && jumps_number[i] == 0 # Get the right ordering for states for deterministic models
                for j = pos:length(states)
                    if occursin(states[j],reordered_equations[i]) == true && (j in states_that_have_been_reordered) == false
                        reordered_states[pos], reordered_states[j] = reordered_states[j], reordered_states[pos]
                        push!(states_that_have_been_reordered,j)
                        pos += 1
                    end
                end
            end
        end
    end

    return reordered_equations, reordered_states

end

function reorganize_equations(equations::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},variables::Array{Q,1},lag_variables::Array{Q,1}) where {Q <: AbstractString}

    #= This function replaces lagged variables with pseudo current variables,
       augmenting the state vector and the model's equations accordingly. =#

    reorganized_equations = copy(equations)
    reorganized_states    = copy(states)

    # First we determine if any equation contains a lagged variable.

    model_has_lags = false
    for i = 1:length(equations)
        for j in lag_variables
            if occursin(j,equations[i]) == true
                model_has_lags = true
                break
            end
        end
    end

    #= If a model contains lagged variables, then we introduce a pseudo variable
       in its place and augment the list of state variables and the set of model
       equations. =#

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

    reorganized_variables = [reorganized_states;jumps]

    return reorganized_equations, reorganized_states, reorganized_variables

end

function get_re_model_primatives(model_array::Array{Q,1}) where {Q <: AbstractString}

    #= This function takes the model-array read from a model file, extracts the
       critical model information, does some basic error checking, and returns
       it in a structure. =#

    states                        = get_variables(model_array,"states:")
    jumps                         = get_variables(model_array,"jumps:")
    shocks                        = get_variables(model_array,"shocks:")
    variables                     = combine_states_and_jumps(states,jumps)
    equations                     = get_equations(model_array,"equations:")
    (parameters, parametervalues, unassigned_parameters) = get_parameters_and_values(model_array,"parameters:")

    for i in [variables;parameters]
        if sum(occursin.(i,equations)) == false
            println("Warning: $i is not in any equation.")
        end
    end

    combined_names = [parameters;variables;shocks]
    if length(unique(combined_names)) != length(combined_names)
        error("Some parameters, variables, or shocks have the same name.")
    end

    reserved_names = ("exp","log","x","p")
    for name in reserved_names
        if name in combined_names
            error("$name cannot be a name for a variable, a shock, or a parameter.")
        end
    end

    lag_variables  = string.(variables,"(-1)")

    reordered_equations, states = reorder_equations(equations,shocks,states,jumps)
    reorganized_equations, states, variables = reorganize_equations(reordered_equations,states,jumps,variables,lag_variables)

    re_model_primatives = REModelPrimatives(states,jumps,shocks,variables,parameters,parametervalues,reorganized_equations,unassigned_parameters)

    return re_model_primatives

end

function repackage_equations(model::ModelPrimatives)

    #= This function is critical for repackaging the model's equations, replacing
       parameter names with values, and numbering variables. =#

    equations       = model.equations
    shocks          = model.shocks
    variables       = model.variables
    parameters      = model.parameters
    parametervalues = model.parametervalues

    repackaged_equations = copy(equations)

    if length(shocks) != 0
        combined_names = [variables;parameters;shocks]
    else
        combined_names = [variables;parameters]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]

    #= Now we go through every equation and replace future variables, variables, and
      shocks with a numbered element of a vector, "x".  We also replace parameter
      names with parameter values. =#

    for j in sorted_combined_names
        if j in variables
            variable_index = findfirst(isequal(j),variables)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],"$j(+1)" => "x[$(length(variables) + variable_index)]")
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(variable_index)]")
            end
        elseif j in parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => parametervalues[parameter_index])
            end
        elseif j in shocks # Okay even if there are no shocks
            shock_index = findfirst(isequal(j),shocks)
            for i = 1:length(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(2*length(variables) + shock_index)]")
            end
        end
    end

    return repackaged_equations

end

function create_steady_state_equations(model::ModelPrimatives)

    # Make the model static by replacing leads and lags with current variables
    # and setting shocks equal to zero.  Also, replace parameter names with
    # their associated value.

    equations       = model.equations
    variables       = model.variables
    shocks          = model.shocks
    parameters      = model.parameters
    parametervalues = model.parametervalues

    steady_state_equations = copy(equations)

    if length(shocks) != 0
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
        elseif j in shocks
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

function create_projection_equations(equations::Array{Q,1},model::ModelPrimatives) where {Q <: AbstractString}

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
    eqns_to_be_approximated  = Int64[]
    for i = 1:ny
        for j = 1:ne
            if occursin("approx$i",projection_equations[j]) == true
                push!(jumps_to_be_approximated,i)
                push!(eqns_to_be_approximated,j)
            end
        end
    end

    jumps_to_be_approximated = unique(jumps_to_be_approximated)

    return projection_equations, jumps_to_be_approximated, eqns_to_be_approximated

end

function create_processed_model_file(model::ModelPrimatives,path::Q) where {Q <: AbstractString}

    # Takes the model's primatives and turns these into a processed-model file.
    # This file is saved as a text file in the same folder as the model file.

    # First, get or construct all the information needed for the processed-model file

    repackaged_equations = repackage_equations(model)

    nonlinear_equations, jumps_to_be_approximated, eqns_to_be_approximated = create_projection_equations(repackaged_equations,model)
    projection_equations = make_equations_equal_zero(nonlinear_equations)

    steady_state_equations = create_steady_state_equations(model)
    static_equations       = make_equations_equal_zero(steady_state_equations)
    dynamic_equations      = make_equations_equal_zero(repackaged_equations)

    number_states    = length(model.states)
    number_jumps     = length(model.jumps)
    number_shocks    = length(model.shocks)
    number_variables = length(model.variables)
    number_equations = length(model.equations)

    variables = model.variables

    # Build up the string containing the processed model information that gets saved

    # First, add the model's summary information

    model_string = "nx = $number_states \n \n"
    model_string = string(model_string,"ny = $number_jumps \n \n")
    model_string = string(model_string,"ns = $number_shocks \n \n")
    model_string = string(model_string,"nv = $number_variables \n \n")
    model_string = string(model_string,"ne = $number_equations \n \n")

    model_string = string(model_string,"jumps_to_approximate = $jumps_to_be_approximated \n \n")
    model_string = string(model_string,"eqns_to_approximate = $eqns_to_be_approximated \n \n")
    model_string = string(model_string,"variables = $variables \n \n")

    # Second, add the model's static information

    if length(model.unassigned_parameters) != 0
        nlsolve_static_string = "function nlsolve_static_equations(f::Array{T,1},x::Array{T,1},p::Array{T1,1}) where {T<:Number, T1<:Real} \n \n"
        static_string         = "function static_equations(x::Array{T,1},p::Array{T1,1}) where {T<:Number, T1<:Real} \n \n"
    else
        nlsolve_static_string = "function nlsolve_static_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n"
        static_string         = "function static_equations(x::Array{T,1}) where {T<:Number} \n \n"
    end
    static_string         = string(static_string,"  f = Array{T,1}(undef,length(x)) \n \n")
    for i = 1:length(static_equations)
        static_string         = string(static_string,"  f[$i] = ",static_equations[i],"\n")
        nlsolve_static_string = string(nlsolve_static_string,"  f[$i] = ",static_equations[i],"\n")
    end

    static_string         = string(static_string,"\n  return f \n \n","end")
    nlsolve_static_string = string(nlsolve_static_string,"\n","end")

    model_string = string(model_string,nlsolve_static_string," \n \n",static_string," \n \n")

    # Third, add the model's dynamic information for perturbation solvers

    if length(model.unassigned_parameters) != 0
        dynamic_string = "function dynamic_equations(x::Array{T,1},p::Array{T1,1}) where {T<:Number, T1<:Real} \n \n"
    else
        dynamic_string = "function dynamic_equations(x::Array{T,1}) where {T<:Number} \n \n"
    end
    dynamic_string = string(dynamic_string,"  f = Array{T,1}(undef,$number_equations) \n \n")
    for i = 1:number_equations
        dynamic_string = string(dynamic_string,"  f[$i] = ",dynamic_equations[i],"\n")
    end

    dynamic_string = string(dynamic_string,"\n  return f \n \n","end \n")

    each_equation_string = Array{String}(undef,number_equations)

    for i = 1:number_equations
        if length(model.unassigned_parameters) != 0
            each_equation_string[i] = "function dynamic_eqn_$i(x::Array{T,1},p::Array{T1,1}) where {T<:Number, T1<:Real} \n \n"
        else
            each_equation_string[i] = "function dynamic_eqn_$i(x::Array{T,1}) where {T<:Number} \n \n"
        end
        each_equation_string[i] = string(each_equation_string[i],"  f = ", dynamic_equations[i],"\n","\n  return f \n \n","end \n")
    end

    for i = 1:number_equations
        dynamic_string = string(dynamic_string, "\n", each_equation_string[i])
    end

    individual_equations_string = "individual_equations = Array{Function}(undef,$number_equations) \n"
    for i = 1:number_equations
        individual_equations_string = string(individual_equations_string,"individual_equations[$i] = dynamic_eqn_$i","\n")
    end

    dynamic_string = string(dynamic_string,"\n",individual_equations_string)
    model_string = string(model_string,dynamic_string)

    # Fourth, add the model's dynamic information for projection solvers

    # For Chebyshev and Smolyak

    if length(model.unassigned_parameters) != 0
        closure_string = "function closure_projection_equations(state,scaled_weights,order,domain,approximate,p) \n \n"
    else
        closure_string = "function closure_projection_equations(state,scaled_weights,order,domain,approximate) \n \n"
    end
    closure_string = string(closure_string, "  function projection_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    weight_number = 1
    for i in jumps_to_be_approximated
        closure_string = string(closure_string, "    approx$i = approximate(scaled_weights[$weight_number],x[$number_jumps+1:end],order,domain)","\n")
        weight_number += 1
    end
    closure_string = string(closure_string,"\n","    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i = 1:length(projection_equations)
        closure_string = string(closure_string,"    f[$i] = ", projection_equations[i], "\n")
    end
    closure_string = string(closure_string,"\n    #return f \n \n  end \n \n  return projection_equations \n \n","end \n")

    # For piecewise linear

    if length(model.unassigned_parameters) != 0
        if number_shocks == 0  # We need to separate the function generated for the stochastic and deterministic cases
            closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,approximate,p) \n \n"
        else
            closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,integrals,approximate,p) \n \n"
        end
    else
        if number_shocks == 0  # We need to separate the function generated for the stochastic and deterministic cases
            closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,approximate) \n \n"
        else
            closure_pl_string = "function closure_projection_equations_pl(variables,grid,state,integrals,approximate) \n \n"
        end
    end

    closure_pl_string = string(closure_pl_string,"  function projection_equations_pl(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    for i in jumps_to_be_approximated
        if number_shocks == 0
            closure_pl_string = string(closure_pl_string,"    approx$i = approximate(variables[$i],grid,x[$number_jumps+1:end])","\n")
        else
            closure_pl_string = string(closure_pl_string,"    approx$i = approximate(variables[$i],grid,x[$number_jumps+1:end],integrals)","\n")
        end
    end

    closure_pl_string = string(closure_pl_string,"\n","    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i = 1:length(projection_equations)
        closure_pl_string = string(closure_pl_string,"    f[$i] = ",projection_equations[i], "\n")
    end

    closure_pl_string = string(closure_pl_string,"\n    #return f \n \n  end \n \n  return projection_equations_pl \n \n","end \n")

    model_string = string(model_string,"\n",closure_string)
    model_string = string(model_string,"\n",closure_pl_string)
    model_string = string(model_string,"\n","unassigned_parameters = $(model.unassigned_parameters)" )

    model_path = replace(path,".txt" => "_processed.txt")
    open(model_path, "w") do io write(io, model_string) end

end

function process_re_model(model_array::Array{Q,1},path::Q) where {Q <: AbstractString}

    # Creates the processed model structure for rational expectations models
    # (anticipating that other types of models may come later).

    re_model_primatives = get_re_model_primatives(model_array)

    create_processed_model_file(re_model_primatives,path)

    println("The model's variables are now in this order: ",re_model_primatives.variables)
    if length(re_model_primatives.unassigned_parameters) != 0
        println("The following parameters do not have values assigned: $(re_model_primatives.unassigned_parameters)")
    end

end

function process_model(path::Q) where {Q <: AbstractString}

    # Main function used to open, read, and process a model file.  The processed model
    # in written to a file that contains all the information needed for the model
    # solvers.

    model_array = open_model_file(path)

    process_re_model(model_array, path)

end

function retrieve_processed_model(path::Q) where {Q <: AbstractString}

    if !occursin("_processed",path)
        path = replace(path,".txt" => "_processed.txt")
    end

    include(path) # The information included is placed in the global scope

    if length(unassigned_parameters) != 0
        dsge_model = REModelPartial(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_projection_equations,closure_projection_equations_pl,unassigned_parameters)
    else
        dsge_model = REModel(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_projection_equations,closure_projection_equations_pl)
    end

    return dsge_model

end

function assign_parameters(model,param::Array{T,1}) where {T <: Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    jumps_approx = model.jumps_approximated
    eqns_approx = model.eqns_approximated
    vars = model.variables

    nlsse(f,x) = model.nlsolve_static_function(f,x,param)
    sf(x)      = model.static_function(x,param)
    df(x)      = model.dynamic_function(x,param)

    ief        = Array{Function}(undef,ne)
    for i = 1:ne
        ffie(x) = model.each_eqn_function[i](x,param)
        ief[i] = ffie
    end

    cf(state,scaled_weights,order,domain,approximate) = model.closure_function(state,scaled_weights,order,domain,approximate,param)
    cfpl_stoch(variables,grid,state,integrals,approximate)  = model.closure_function_piecewise(variables,grid,state,integrals,approximate,param)
    cfpl_det(variables,grid,state,approximate)  = model.closure_function_piecewise(variables,grid,state,approximate,param)

    if ns != 0
        newmod = REModel(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,vars,nlsse,sf,df,ief,cf,cfpl_stoch)
        return newmod
    else
        newmod = REModel(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,vars,nlsse,sf,df,ief,cf,cfpl_det)
        return newmod
    end

end
