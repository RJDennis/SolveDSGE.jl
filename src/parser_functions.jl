################# Parser functions ######################

"""
Opens the model file and reads the contents, which are stored in a vector of strings.

This function also deals with comments, empty lines, and blank lines.

Internal function; not exposed to users.
"""
function open_model_file(path::Q) where {Q<:AbstractString}

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

    # Identify lines that contain only spaces

    for i in eachindex(model_array)
        if unique(model_array[i]) == [' '] # line contains only spaces
            model_array[i] = "" # Make it a blank line
        end
    end

    # Remove blank lines

    model_array = model_array[model_array.!=""]

    return model_array

end

"""
Finds the position in the model's string-vector where model-terms (states, jumps,
parameters, equations, etc) are located.

Internal function; not exposed to users.
"""
function find_term(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    locations = findall(y -> contains(y,term),model_array)
    if length(locations) == 1
        return locations[1]
    elseif length(locations) == 0
        error("The $term-designation does not appear in the model file.")
    else
        error("The $term-designation appears multiple times in the model file.")
    end

end

"""
Finds the position in the model's string-vector where the model-term 'end' appears.  The
function searches the model's string-vector beginning from line 'startfrom' and returns 
the line where the next 'end' occurs, signifying the end of a specification block.

White space at the start of a line is ignored.

Internal function; not exposed to users.
"""
function find_end(model_array::Array{Q,1},startfrom::S) where {Q<:AbstractString,S<:Integer}

    for end_location = startfrom:length(model_array)
        if startswith(strip(model_array[end_location]),"end") == true
            return end_location
        end
    end

    error("Model file is missing an 'end' term.")

end

"""
Extracts the names from within a specification block of a model's string-vector 
and checks that no names are repeated.  Used to capture the names of jumps:, 
states:, and shocks:. 

Internal function; not exposed to users.
"""
function get_variables(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    term_begin = find_term(model_array,term) + 1 # Starts the line after the term is located
    term_end = find_end(model_array,term_begin) - 1 # Ends the line before 'end' is located

    if term_begin > term_end
        if term in ["shocks:","states:"]
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    end

    term_block = model_array[term_begin:term_end]

    # Remove any trailing variable separators: "," or ";".

    for i in eachindex(term_block)
        if endswith(term_block[i],union(",",";"))
            term_block[i] = term_block[i][1:end-1]
        end
    end

    # Extract the names and place them in a vector

    variables = String.(setdiff(strip.(split(term_block[1],union(",",";"))),[""]))
    for i = 2:length(term_block)
        variables = [variables; String.(setdiff(strip.(split(term_block[i],union(",",";"))),[""]))]
    end

    # Check whether names are repeated

    if length(variables) != length(unique(variables))
        error("Some $(term[1:end-1]) are repreated.")
    end

    # Check to ensure that variables contains non-empty string elements

    if length(variables) == 1 && variables[1] == ""
        if term in ["shocks:", "states:"]
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    else
        return variables
    end

end

"""
Combines the states, "x", with the jump variables, "y", into a single vector to
generate a vector of the model's variables.

Internal function; not exposed to users.
"""
function combine_states_and_jumps(x::Array{Q,1},y::Array{Q,1}) where {Q<:AbstractString}


    if length(x) == 0 # There are no states
        return y
    elseif length(intersect(x,y)) > 0
        error("Some states and jumps have the same name.")
    else
        return [x; y]
    end

end

"""
Extracts the names and values for each of the model's parameters.  The parameter names are 
sorted so that larger names come first.

Internal function; not exposed to users.
"""
function get_parameters_and_values(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    parametersbegin = find_term(model_array,term) + 1 # Starts the line after the term is located
    parametersend = find_end(model_array,parametersbegin) - 1 # Ends the line before 'end' is located
    if parametersbegin > parametersend
        error("The model file contains no $(term[1:end-1])")
    end

    parameterblock = model_array[parametersbegin:parametersend]

    # Remove any trailing separators: "," or ";".

    for i in eachindex(parameterblock)
        if endswith(parameterblock[i],union(",",";"))
            parameterblock[i] = parameterblock[i][1:end-1]
        end
    end

    revised_parameterblock = String.(strip.(split(parameterblock[1],union(",",";"))))
    for i = 2:length(parameterblock)
        revised_parameterblock = [revised_parameterblock; String.(strip.(split(parameterblock[i],union(",",";"))))]
    end

    # Extract the parameter names and values, making note of any parameters with unassigned values

    unassigned_parameter_index = 1
    unassigned_parameters = Array{Q}(undef,0)

    parameters = Array{Q}(undef,length(revised_parameterblock))
    values = Array{Q}(undef,length(revised_parameterblock))
    for i in eachindex(revised_parameterblock)
        if occursin("=",revised_parameterblock[i]) == false
            parameters[i] = revised_parameterblock[i]
            values[i] = "p[$unassigned_parameter_index]" # p is a reserved name
            push!(unassigned_parameters,revised_parameterblock[i])
            unassigned_parameter_index += 1
        else
            pair = strip.(split(revised_parameterblock[i],"="))
            parameters[i] = pair[1]
            values[i]     = pair[2]
        end
    end

    parameter_order = sortperm(length.(parameters),rev = true)
    sorted_parameters = parameters[parameter_order]
    sorted_values = values[parameter_order]

    return sorted_parameters, sorted_values, unassigned_parameters

end

"""
Extracts the model's equations.

Internal function; not exposed to users.
"""
function get_equations(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    equationsbegin = find_term(model_array,term) + 1 # Starts the line after the term is located
    equationsend = find_end(model_array,equationsbegin) - 1 # Ends the line before 'end' is located
    if equationsbegin > equationsend
        error("The model file contains no $(term[1:end-1])")
    end

    equation_block = model_array[equationsbegin:equationsend]

    # Remove any trailing separators: "," or ";".

    for i in eachindex(equation_block)
        if endswith(equation_block[i],union(",",";"))
            equation_block[i] = equation_block[i][1:end-1]
        end
    end

    # Extract the equations and place them in a vector

    equations = String.(strip.(split(equation_block[1],union(",",";"))))
    for i = 2:length(equation_block)
        equations = [equations; String.(strip.(split(equation_block[i],union(",",";"))))]
    end

    # For every model equation...

    for i in eachindex(equation_block)
        if occursin("[",equations[i]) == true # Replace open square bracket with open round parenthesis
            equations[i] = replace(equations[i],"[" => "(")
        elseif occursin("]",equations[i]) == true # Replace close square bracket with close round parenthesis
            equations[i] = replace(equations[i],"]" => ")")
        elseif occursin("{",equations[i]) == true # Replace open curly brace with open round parenthesis
            equations[i] = replace(equations[i],"{" => "(")
        elseif occursin("}",equations[i]) == true # Replace close curly brace with close round parenthesis
            equations[i] = replace(equations[i],"}" => ")")
        end
        if occursin("=",equations[i]) == false # Check that each equation contains an equals sign
            error("Equation line $i does not contain an '=' sign.")
        elseif length(findall("(",equations[i])) != length(findall(")",equations[i]))
            error("Equation line $i has unbalanced parentheses")
        end
    end

    return equations

end

"""
Extracts the solvers: designation, which determines which solvers can be applied to this model.

Internal function; not exposed to users.
"""
function get_solvers(model_array::Array{Q,1}, term::Q) where {Q<:AbstractString}

    locations = findall(y -> contains(y, term), model_array)
    if length(locations) == 0
        return "Any" # If solvers: is not specified it defaults to "Any"
    elseif length(locations) > 1
        error("The $term-designation appears multiple times in the model file.")
    else
        solvers_line = locations[1]
    end

    if occursin("Linear", model_array[solvers_line]) || occursin("linear", model_array[solvers_line])
        return "Linear"
    elseif occursin("Perturbation", model_array[solvers_line]) || occursin("perturbation", model_array[solvers_line])
        return "Perturbation"
    elseif occursin("Projection", model_array[solvers_line]) || occursin("projection", model_array[solvers_line])
        return "Projection"
    elseif occursin("Any", model_array[solvers_line]) || occursin("any", model_array[solvers_line])
        return "Any"
    else
        error("The specified solver is unknown")
    end
end

"""
Reorders the model's equation so that the shock processes are at the top.

Internal function; not exposed to users.
"""
function reorder_equations(equations::Array{Q,1},shocks::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},parameters::Array{Q,1}) where {Q<:AbstractString}

    if isempty(shocks) == false # Case for stochastic models
        reordered_equations, reordered_states, reordered_shocks = _reorder_equations(equations,shocks,states,jumps,parameters)
        return reordered_equations, reordered_states, reordered_shocks
    else # Case for deterministic models
        reordered_equations, reordered_states = _reorder_equations(equations,states,jumps,parameters)
        return reordered_equations, reordered_states, shocks
    end

end

"""
Reorders the equations of a deteministic model so that the shock processes are at the top.

Internal function; not exposed to users.
"""
function _reorder_equations(equations::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},parameters::Array{Q,1}) where {Q<:AbstractString}

     reordered_equations = copy(equations)
     reordered_states = copy(states)
   
    # Construct summary information about each equation
   
    junk_equations = copy(equations)
    states_left_in_eqns = copy(equations)
    combined_names = [states; jumps; parameters; ["exp", "log"]]
    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]
   
    states_number = zeros(Int64,length(equations))
    jumps_number  = zeros(Int64,length(equations))
   
    for i in eachindex(junk_equations)
        for j in sorted_combined_names
            if occursin(j,junk_equations[i]) == true
                if j in states
                    states_number[i] += 1
                elseif j in jumps
                    jumps_number[i] += 1
                    states_left_in_eqns[i] = replace(states_left_in_eqns[i],j => ":")
                else
                    states_left_in_eqns[i] = replace(states_left_in_eqns[i],j => ":")
                end
                junk_equations[i] = replace(junk_equations[i],j => ":")
            end
        end
    end
    number_eqns_with_no_jumps = sum(jumps_number .== 0)

    # Put the equations with no jumps in them at the top

    ind = sortperm(jumps_number)
    reordered_equations .= reordered_equations[ind]
    states_left_in_eqns .= states_left_in_eqns[ind]
    states_number .= states_number[ind]
    jumps_number  .= jumps_number[ind]

    # Now sort out the order of the states in the system

    states_that_have_been_ordered = Int64[]
    for k in eachindex(states)
        for j = 1:number_eqns_with_no_jumps
            if states_number[j] == k
                for i in eachindex(states)
                    if occursin(reordered_states[i],states_left_in_eqns[j]) == true && j != i && (i in states_that_have_been_ordered) == false
                        reordered_states[i], reordered_states[j] = reordered_states[j], reordered_states[i]
                        push!(states_that_have_been_ordered,j)
                        break
                    end
                end
            end
        end
    end

    return reordered_equations, reordered_states

end

"""
Reorders the equations of a stochastic model so that the shock processes are at the top.

Internal function; not exposed to users.
"""
function _reorder_equations(equations::Array{Q,1},shocks::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},parameters::Array{Q,1}) where {Q<:AbstractString}

    reordered_equations = copy(equations)
    reordered_states = copy(states)
    reordered_shocks = copy(shocks)

    # Construct summary information about each equation

    junk_equations = copy(equations)
    states_left_in_eqns = copy(equations)
    combined_names = [states; jumps; parameters; shocks; ["exp", "log"]]
    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]

    shocks_number = zeros(Int64,length(equations))
    states_number = zeros(Int64,length(equations))
    jumps_number  = zeros(Int64,length(equations))

    for i in eachindex(junk_equations)
        for j in sorted_combined_names
            if occursin(j,junk_equations[i]) == true
                if j in states
                    states_number[i] += 1
                elseif j in jumps
                    jumps_number[i] += 1
                    states_left_in_eqns[i] = replace(states_left_in_eqns[i],j => ":")
                elseif j in shocks
                    shocks_number[i] += 1
                    states_left_in_eqns[i] = replace(states_left_in_eqns[i],j => ":")
                else
                    states_left_in_eqns[i] = replace(states_left_in_eqns[i],j => ":")
                end
                junk_equations[i] = replace(junk_equations[i],j => ":")
            end
        end
    end
    number_eqns_with_shocks = sum(shocks_number .!= 0)

    # Put the equations with no jumps in them at the top

    ind = sortperm(jumps_number)
    reordered_equations .= reordered_equations[ind]
    states_left_in_eqns .= states_left_in_eqns[ind]
    states_number .= states_number[ind]
    jumps_number  .= jumps_number[ind]
    shocks_number .= shocks_number[ind]

    # Put the shock equations with the fewest shocks at the top

    ind = sortperm(shocks_number[1:number_eqns_with_shocks])
    reordered_equations[1:number_eqns_with_shocks] .= reordered_equations[ind]
    states_left_in_eqns[1:number_eqns_with_shocks] .= states_left_in_eqns[ind]
    shocks_number[1:number_eqns_with_shocks] .= shocks_number[ind]
    states_number[1:number_eqns_with_shocks] .= states_number[ind]
    jumps_number[1:number_eqns_with_shocks] .= jumps_number[ind]

    # Order the shock processes to try to make shock's variance-covariance matrix have non-zero diagonals.

    shocks_that_have_been_ordered = Int64[]
    for k in eachindex(shocks)
        for j = 1:number_eqns_with_shocks
            if shocks_number[j] == k
                for i in eachindex(shocks)
                    if occursin(reordered_shocks[i],reordered_equations[j]) == true && j == i && (i in shocks_that_have_been_ordered) == false
                        push!(shocks_that_have_been_ordered,i)
                        break
                    elseif occursin(reordered_shocks[i],reordered_equations[j]) == true && j != i && (i in shocks_that_have_been_ordered) == false
                        reordered_shocks[i], reordered_shocks[j]       = reordered_shocks[j], reordered_shocks[i]
                        states_left_in_eqns[i], states_left_in_eqns[j] = states_left_in_eqns[j], states_left_in_eqns[i]
                        push!(shocks_that_have_been_ordered,i)
                        break
                    end
                end
            end
        end
    end

    # Now sort out the order of the states in the system

    states_that_have_been_ordered = Int64[]
    for k in eachindex(states)
        for j = 1:number_eqns_with_shocks
            if states_number[j] == k
                for i in eachindex(states)
                    if occursin(reordered_states[i],states_left_in_eqns[j]) == true && j != i && (i in states_that_have_been_ordered) == false
                        reordered_states[i], reordered_states[j] = reordered_states[j], reordered_states[i]
                        push!(states_that_have_been_ordered,j)
                        break
                    end
                end
            end
        end
    end

    return reordered_equations, reordered_states, reordered_shocks

end

"""
Replaces lagged variables with pseudo current variables, augmenting the state vector and 
the model's equations accordingly.

Internal function; not exposed to users.
"""
function deal_with_lags(equations::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},variables::Array{Q,1}) where {Q<:AbstractString}

    lag_variables = string.(variables,"(-1)")

    reorganized_equations = copy(equations)
    reorganized_states = copy(states)

    # First we determine if any equation contains a lagged variable.

    model_has_lags = false
    for i in eachindex(equations)
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
        for j in eachindex(lag_variables)
            if sum(occursin.(lag_variables[j],equations)) != 0
                for i in eachindex(equations)
                    reorganized_equations[i] = replace(reorganized_equations[i],lag_variables[j] => string(variables[j],"lag"))
                end
                reorganized_states = [reorganized_states; string(variables[j],"lag")]
                reorganized_equations = [reorganized_equations; string(variables[j],"lag(+1) = ",variables[j])]
            end
        end
    end

    reorganized_variables = [reorganized_states; jumps]

    return reorganized_equations, reorganized_states, reorganized_variables

end

"""
Takes the model-array read from the model file for a rational expectations model, 
extracts the critical model information, does some basic error checking, and 
returns the model's key information in a structure.

Internal function; not exposed to users.
"""
function get_re_model_primatives(model_array::Array{Q,1}) where {Q<:AbstractString}

    states = get_variables(model_array,"states:")
    jumps = get_variables(model_array,"jumps:")
    shocks = get_variables(model_array,"shocks:")
    variables = combine_states_and_jumps(states,jumps)
    equations = get_equations(model_array,"equations:")
    (parameters, parametervalues, unassigned_parameters) = get_parameters_and_values(model_array,"parameters:")
    solvers = get_solvers(model_array,"solvers:")

    for i in [variables; parameters]
        if i in variables
            if sum(occursin.(i,equations)) == false
                error("Variable $i is not in any equation.")
            end
        else
            if sum(occursin.(i,equations)) + sum(occursin.(i,parametervalues)) == 0
                println("Warning: parameter $i is not in any equation.")
            end
        end
    end

    combined_names = [parameters; variables; shocks]
    if length(unique(combined_names)) != length(combined_names)
        error("Some parameters, variables, or shocks have the same name.")
    end

    reserved_names = ("deriv", "exp", "log", "x", "p", ":", ";")
    for name in reserved_names
        if name in combined_names
            error("$name cannot be the name for a variable, a shock, or a parameter.")
        end
    end

    reordered_equations, states, shocks = reorder_equations(equations,shocks,states,jumps,parameters)
    reorganized_equations, states, variables = deal_with_lags(reordered_equations,states,jumps,variables)

    re_model_primatives = REModelPrimatives(states,jumps,shocks,variables,parameters,parametervalues,reorganized_equations,unassigned_parameters,solvers)

    return re_model_primatives

end

"""
Repackages the model's equations, replaces parameter names with parameter values, and
determines the order of the variables in the system.

Internal function; not exposed to users.
"""
function repackage_equations(model::DSGEModelPrimatives)

    equations = model.equations
    shocks = model.shocks
    variables = model.variables
    parameters = model.parameters
    parametervalues = model.parametervalues

    repackaged_equations       = copy(equations)
    repackaged_parametervalues = copy(parametervalues)

    if length(shocks) != 0
        combined_names = [variables; parameters; shocks]
    else
        combined_names = [variables; parameters]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]
    sorted_parameters     = parameters[sortperm(length.(parameters),rev = true)]

    #= First we go through every equation and replace exp with : and log with ;.  
       This is to guard them during variables and parameter substitution. =#

    for i in eachindex(repackaged_equations)
        if occursin("exp",equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],"exp" => ":")
        elseif occursin("log",equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],"log" => ";")
        end
    end

    #= Next we go through every parameter expression and replace exp with : and log with ;.  
    This is to guard them during variables and parameter substitution. =#

    for i in eachindex(parametervalues)
        if occursin("exp",parametervalues[i]) == true
            repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],"exp" => ":")
        elseif occursin("log",parametervalues[i]) == true
            repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],"log" => ";")
        end
    end

    #= Here we go through all parameters and deal with parameters depending on other
       parameters =#

    loops = 0 # Counts the number of loops over the parameters
    while true
        count = 0 # counts whether parameter values are still being assigned
        for j in sorted_parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i in eachindex(repackaged_parametervalues)
                if occursin(j,repackaged_parametervalues[i]) == true
                    repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],j => string("(",repackaged_parametervalues[parameter_index],")"))
                    count += 1
                end
            end
        end
        loops += 1
        if count == 0
            break
        end
        if loops > length(parameters)-1
            error("There is a circularity in the parameter definitions")
        end
    end

    #= Now we go through every equation and replace future variables, variables, and
       shocks with a numbered element of a vector, "x".  We also replace parameter names 
       with parameter values. =#

    for j in sorted_combined_names
        if j in variables
            variable_index = findfirst(isequal(j),variables)
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],"$j(+1)" => "x[$(length(variables) + variable_index)]")
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(variable_index)]")
            end
        elseif j in parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => string("(",j,")"))
                repackaged_equations[i] = replace(repackaged_equations[i],j => repackaged_parametervalues[parameter_index])
            end
        elseif j in shocks # Okay even if there are no shocks
            shock_index = findfirst(isequal(j),shocks)
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => "x[$(2*length(variables) + shock_index)]")
            end
        end
    end

    #= Finally, go back through every equation and restore exp and log where necessary =#

    for i in eachindex(repackaged_equations)
        if occursin(":",repackaged_equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],":" => "exp")
        end
        if occursin(";",repackaged_equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],";" => "log")
        end
    end

    return repackaged_equations

end

"""
Make the model static by replacing leads and lags with current variables
and setting shocks equal to zero.  Also, replaces parameter names with
their associated value.

Internal function; not exposed to users.
"""
function create_steady_state_equations(model::DSGEModelPrimatives)

    equations = model.equations
    variables = model.variables
    shocks = model.shocks
    parameters = model.parameters
    parametervalues = model.parametervalues

    steady_state_equations       = copy(equations)
    steady_state_parametervalues = copy(parametervalues)

    if length(shocks) != 0
        combined_names = [variables; parameters; shocks]
    else
        combined_names = [variables; parameters]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]
    sorted_parameters     = parameters[sortperm(length.(parameters),rev = true)]

    #= First we go through every equation and replace exp with : and log with ;.  
    This is to guard them during variables and parameter substitution. =#

    for i in eachindex(steady_state_equations)
        if occursin("exp",equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],"exp" => ":")
        elseif occursin("log",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],"log" => ";")
        end
    end

    #= Next we go through every parameter expression and replace exp with : and log with ;.  
    This is to guard them during variables and parameter substitution. =#

    for i in eachindex(parametervalues)
        if occursin("exp",parametervalues[i]) == true
            steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],"exp" => ":")
        elseif occursin("log",parametervalues[i]) == true
            steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],"log" => ";")
        end
    end

    #= Now we take care of the fact that some model parameters may be functions of deeper
    behavioral parameters =#

    loops = 0 # Counts the number of loops over the parameters
    while true
        count = 0 # counts whether parameter values are still being assigned
        for j in sorted_parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i in eachindex(steady_state_parametervalues)
                if occursin(j,steady_state_parametervalues[i]) == true
                    steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],j => string("(",steady_state_parametervalues[parameter_index],")"))
                    count += 1
                end
            end
        end
        loops += 1
        if count == 0
            break
        end
        if loops > length(parameters)-1
            error("There is a circularity in the parameter definitions")
        end
    end

    # Now we go through every equation and replace future variables, variables, and
    # shocks with a numbered element of a vector, "x".  We also replace parameter
    # names with parameter values

    for j in sorted_combined_names
        if j in variables
            variable_index = findfirst(isequal(j),variables)
            for i in eachindex(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],"$j(-1)" => "x[$(variable_index)]")
                steady_state_equations[i] = replace(steady_state_equations[i],"$j(+1)" => "x[$(variable_index)]")
                steady_state_equations[i] = replace(steady_state_equations[i],j => "x[$(variable_index)]")
            end
        elseif j in parameters
            parameter_index = findfirst(isequal(j),parameters)
            for i in eachindex(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => string("(",j,")"))
                steady_state_equations[i] = replace(steady_state_equations[i],j => steady_state_parametervalues[parameter_index])
            end
        elseif j in shocks
            shock_index = findfirst(isequal(j),shocks)
            for i in eachindex(equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => 0.0)
            end
        end
    end

    #= Finally, go back through every equation and restore exp and log where necessary =#

    for i in eachindex(steady_state_equations)
        if occursin(":",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],":" => "exp")
        end
        if occursin(";",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],";" => "log")
        end
    end

    return steady_state_equations

end

"""
Reexpress the model's equations such that they each equal zero.

Internal function; not exposed to users.
"""
function make_equations_equal_zero(equations::Array{Q,1}) where {Q<:AbstractString}

    zeroed_equations = similar(equations)

    for i in eachindex(equations)
        pair = strip.(split(equations[i],"="))
        zeroed_equations[i] = string(pair[1]," - (",pair[2],")")
    end

    return zeroed_equations

end

"""
Specific to projection methods.  This function expresses the equations such that they can be solved 
using the projection solvers.

Creates an inventory of which variables have functions that need approximating and in which 
equations they reside.

Internal function; not exposed to users.
"""
function create_projection_equations(equations::Array{Q,1},model::DSGEModelPrimatives) where {Q<:AbstractString}

    projection_equations = copy(equations)

    nx = length(model.states)
    ny = length(model.jumps)
    ns = length(model.shocks)
    ne = length(equations)
    nv = nx + ny

    for j = 1:nx
        for i = 1:ne
            projection_equations[i] = replace(projection_equations[i],"x[$j]" => "state[$j]")
        end
    end

    for j = 1:nv
        for i = 1:ne
            projection_equations[i] = replace(projection_equations[i],"x[$(nx+j)]" => "x[$j]")
        end
    end

    for j = 1:ny
        for i = 1:ne
            projection_equations[i] = replace(projection_equations[i],"x[$(nx+nv+j)]" => "approx$j")
        end
    end

    for j = 1:ns
        for i = 1:ne
            projection_equations[i] = replace(projection_equations[i],"x[$(2*nv+j)]" => 0.0)
        end
    end

    jumps_to_be_approximated = Int64[]
    eqns_to_be_approximated = Int64[]
    for i = 1:ny
        for j = 1:ne
            if occursin("approx$i",projection_equations[j]) == true
                push!(jumps_to_be_approximated,i)
                push!(eqns_to_be_approximated,j)
            end
        end
    end

    jumps_to_be_approximated = unique(jumps_to_be_approximated)
    eqns_to_be_approximated = sort(unique(eqns_to_be_approximated))

    return projection_equations, jumps_to_be_approximated, eqns_to_be_approximated

end

"""
Creates and saves the file for the processed model.

Internal funtion; not exposed to users.
"""
function create_processed_model_file(model::DSGEModelPrimatives, path::Q) where {Q<:AbstractString}

    # Takes the model's primatives and turns these into a processed-model file.
    # This file is saved as a text file in the same folder as the model file.

    # First, get or construct all the information needed for the processed-model file

    repackaged_equations = repackage_equations(model)

    nonlinear_equations, jumps_to_be_approximated, eqns_to_be_approximated = create_projection_equations(repackaged_equations, model)
    projection_equations = make_equations_equal_zero(nonlinear_equations)

    steady_state_equations = create_steady_state_equations(model)
    static_equations = make_equations_equal_zero(steady_state_equations)
    dynamic_equations = make_equations_equal_zero(repackaged_equations)

    number_states = length(model.states)
    number_jumps = length(model.jumps)
    number_shocks = length(model.shocks)
    number_variables = length(model.variables)
    number_equations = length(model.equations)

    #variables = model.variables
    #variables = sort(Dict(model.variables[i] => i for i = 1:number_variables), byvalue = :true)
    variables = OrderedDict(model.variables[i] => i for i = 1:number_variables)

    # Build up the string containing the processed model information that gets saved

    # First, add the model's summary information

    model_string = "nx = $number_states \n \n"
    model_string = string(model_string, "ny = $number_jumps \n \n")
    model_string = string(model_string, "ns = $number_shocks \n \n")
    model_string = string(model_string, "nv = $number_variables \n \n")
    model_string = string(model_string, "ne = $number_equations \n \n")

    model_string = string(model_string, "jumps_to_approximate = $jumps_to_be_approximated \n \n")
    model_string = string(model_string, "eqns_to_approximate = $eqns_to_be_approximated \n \n")
    model_string = string(model_string, "variables = $variables \n \n")

    # Second, add the model's static information

    if length(model.unassigned_parameters) != 0
        nlsolve_static_string = "function nlsolve_static_equations(f::Array{T,1},x::Array{T,1},p::Array{T1,1}) where {T<:Number,T1<:Real} \n \n"
        static_string = "function static_equations(x::Array{T,1},p::Array{T1,1}) where {T<:Number,T1<:Real} \n \n"
    else
        nlsolve_static_string = "function nlsolve_static_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n"
        static_string = "function static_equations(x::Array{T,1}) where {T<:Number} \n \n"
    end
    static_string = string(static_string, "  f = Array{T,1}(undef,length(x)) \n \n")
    for i in eachindex(static_equations)
        nlsolve_static_string = string(nlsolve_static_string, "  f[$i] = ", static_equations[i], "\n")
        static_string = string(static_string, "  f[$i] = ", static_equations[i], "\n")
    end

    nlsolve_static_string = string(nlsolve_static_string, "\n", "end")
    static_string = string(static_string, "\n  return f \n \n", "end")

    model_string = string(model_string, nlsolve_static_string, " \n \n", static_string, " \n \n")

    # Third, add the model's dynamic information for perturbation solvers

    if length(model.unassigned_parameters) != 0
        dynamic_string = "function dynamic_equations(x::Array{T,1},p::Array{T1,1}) where {T<:Number,T1<:Real} \n \n"
    else
        dynamic_string = "function dynamic_equations(x::Array{T,1}) where {T<:Number} \n \n"
    end
    dynamic_string = string(dynamic_string, "  f = Array{T,1}(undef,$number_equations) \n \n")
    for i = 1:number_equations
        dynamic_string = string(dynamic_string, "  f[$i] = ", dynamic_equations[i], "\n")
    end

    dynamic_string = string(dynamic_string, "\n  return f \n \n", "end \n")

    each_equation_string = Array{String}(undef, number_equations)

    for i = 1:number_equations
        if length(model.unassigned_parameters) != 0
            each_equation_string[i] = "function dynamic_eqn_$i(x::Array{T,1},p::Array{T1,1}) where {T<:Number,T1<:Real} \n \n"
        else
            each_equation_string[i] = "function dynamic_eqn_$i(x::Array{T,1}) where {T<:Number} \n \n"
        end
        each_equation_string[i] = string(each_equation_string[i], "  f = ", dynamic_equations[i], "\n", "\n  return f \n \n", "end \n")
    end

    for i = 1:number_equations
        dynamic_string = string(dynamic_string, "\n", each_equation_string[i])
    end

    individual_equations_string = "individual_equations = Array{Function}(undef,$number_equations) \n"
    for i = 1:number_equations
        individual_equations_string = string(individual_equations_string, "individual_equations[$i] = dynamic_eqn_$i", "\n")
    end

    dynamic_string = string(dynamic_string, "\n", individual_equations_string)
    model_string = string(model_string, dynamic_string)

    # Fourth, add the model's dynamic information for projection solvers

    # For Chebyshev

    if length(model.unassigned_parameters) != 0
        closure_cheb_string = "function closure_chebyshev_equations(state,scaled_weights,order,domain,p) \n \n"
    else
        closure_cheb_string = "function closure_chebyshev_equations(state,scaled_weights,order,domain) \n \n"
    end
    closure_cheb_string = string(closure_cheb_string, "  function chebyshev_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    weight_number = 1
    for i in jumps_to_be_approximated
        closure_cheb_string = string(closure_cheb_string, "    approx$i = chebyshev_evaluate(scaled_weights[$weight_number],x[$number_jumps+1:end],order,domain)", "\n")
        weight_number += 1
    end
    closure_cheb_string = string(closure_cheb_string, "\n", "    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i in eachindex(projection_equations)
        closure_cheb_string = string(closure_cheb_string, "    f[$i] = ", projection_equations[i], "\n")
    end
    closure_cheb_string = string(closure_cheb_string, "\n    #return f \n \n  end \n \n  return chebyshev_equations \n \n", "end \n")

    # For Smolyak

    if length(model.unassigned_parameters) != 0
        closure_smol_string = "function closure_smolyak_equations(state,scaled_weights,order,domain,p) \n \n"
    else
        closure_smol_string = "function closure_smolyak_equations(state,scaled_weights,order,domain) \n \n"
    end
    closure_smol_string = string(closure_smol_string, "  function smolyak_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    closure_smol_string = string(closure_smol_string, "    poly = smolyak_polynomial(x[$number_jumps+1:end],order,domain) \n")

    weight_number = 1
    for i in jumps_to_be_approximated
        closure_smol_string = string(closure_smol_string, "    approx$i = smolyak_evaluate(scaled_weights[$weight_number],poly)", "\n")
        weight_number += 1
    end
    closure_smol_string = string(closure_smol_string, "\n", "    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i in eachindex(projection_equations)
        closure_smol_string = string(closure_smol_string, "    f[$i] = ", projection_equations[i], "\n")
    end
    closure_smol_string = string(closure_smol_string, "\n    #return f \n \n  end \n \n  return smolyak_equations \n \n", "end \n")

    # For Hyperbolic-cross

    if length(model.unassigned_parameters) != 0
        closure_hcross_string = "function closure_hcross_equations(state,scaled_weights,order,domain,p) \n \n"
    else
        closure_hcross_string = "function closure_hcross_equations(state,scaled_weights,order,domain) \n \n"
    end
    closure_hcross_string = string(closure_hcross_string, "  function hcross_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    closure_hcross_string = string(closure_hcross_string, "    poly = hyperbolic_cross_polynomial(x[$number_jumps+1:end],order,domain) \n")

    weight_number = 1
    for i in jumps_to_be_approximated
        closure_hcross_string = string(closure_hcross_string, "    approx$i = hyperbolic_cross_evaluate(scaled_weights[$weight_number],poly)", "\n")
        weight_number += 1
    end
    closure_hcross_string = string(closure_hcross_string, "\n", "    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i in eachindex(projection_equations)
        closure_hcross_string = string(closure_hcross_string, "    f[$i] = ", projection_equations[i], "\n")
    end
    closure_hcross_string = string(closure_hcross_string, "\n    #return f \n \n  end \n \n  return hcross_equations \n \n", "end \n")

    # For piecewise linear

    if length(model.unassigned_parameters) != 0
        if number_shocks == 0  # We need to separate the function generated for the stochastic and deterministic cases
            closure_pl_string = "function closure_piecewise_equations(variables,grid,state,p) \n \n"
        else
            closure_pl_string = "function closure_piecewise_equations(variables,grid,state,integrals,p) \n \n"
        end
    else
        if number_shocks == 0  # We need to separate the function generated for the stochastic and deterministic cases
            closure_pl_string = "function closure_piecewise_equations(variables,grid,state) \n \n"
        else
            closure_pl_string = "function closure_piecewise_equations(variables,grid,state,integrals) \n \n"
        end
    end

    closure_pl_string = string(closure_pl_string, "  function piecewise_equations(f::Array{T,1},x::Array{T,1}) where {T<:Number} \n \n")
    for i in jumps_to_be_approximated
        if number_shocks == 0
            closure_pl_string = string(closure_pl_string, "    approx$i = piecewise_linear_evaluate(variables[$i],grid,x[$number_jumps+1:end])", "\n")
        else
            closure_pl_string = string(closure_pl_string, "    approx$i = piecewise_linear_evaluate(variables[$i],grid,x[$number_jumps+1:end],integrals)", "\n")
        end
    end

    closure_pl_string = string(closure_pl_string, "\n", "    #f = Array{T,1}(undef,$number_equations) \n \n")
    for i in eachindex(projection_equations)
        closure_pl_string = string(closure_pl_string, "    f[$i] = ", projection_equations[i], "\n")
    end

    closure_pl_string = string(closure_pl_string, "\n    #return f \n \n  end \n \n  return piecewise_equations \n \n", "end \n")

    model_string = string(model_string, "\n", closure_cheb_string)
    model_string = string(model_string, "\n", closure_smol_string)
    model_string = string(model_string, "\n", closure_hcross_string)
    model_string = string(model_string, "\n", closure_pl_string)
    model_string = string(model_string, "\n", "unassigned_parameters = $(model.unassigned_parameters) \n")
    model_string = string(model_string, "\n", """solvers = "$(model.solvers)" """)

    model_path = replace(path, ".txt" => "_processed.txt")
    open(model_path, "w") do io
        write(io, model_string)
    end

end

"""
Opens, and processes the contents of a model file for a rational expectations model.

This functions exists anticipating that models that depart from rational expections 
will be introduced at some point 

Internal function; not exposed to users
"""
function process_re_model(model_array::Array{Q,1},path::Q) where {Q<:AbstractString}

    # Creates the processed model structure for rational expectations models
    # (anticipating that other types of models may come later).

    re_model_primatives = get_re_model_primatives(model_array)

    create_processed_model_file(re_model_primatives,path)

    println("The model's variables are now in this order: ",re_model_primatives.variables)
    println("The model's shocks are now in this order:    ",re_model_primatives.shocks)
    if length(re_model_primatives.unassigned_parameters) != 0
        println("The following parameters do not have values assigned: $(re_model_primatives.unassigned_parameters)")
    end

end

"""
Opens, and processes the contents of a model file.

Exported function.
"""
function process_model(path::Q) where {Q<:AbstractString}

    # Main function used to open, read, and process a model file.  The processed model
    # in written to a file that contains all the information needed for the model
    # solvers.

    model_array = open_model_file(path)

    process_re_model(model_array,path)

end

"""
Retrives and stores in a structure the information from a processed model file.

Exported function.
"""
function retrieve_processed_model(path::Q) where {Q<:AbstractString}

    if !occursin("_processed",path)
        path = replace(path,".txt" => "_processed.txt")
    end

    include(path) # The information included is placed in the global scope, but then put in a struct

    if length(unassigned_parameters) != 0
        dsge_model = REModelPartial(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations,unassigned_parameters,solvers)
    else
        dsge_model = REModel(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations,solvers)
    end

    return dsge_model

end

"""
Assigns values to parameters previously without values.

Exported function.
"""
function assign_parameters(model::REModelPartial,param::Array{T,1}) where {T<:Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    jumps_approx = model.jumps_approximated
    eqns_approx = model.eqns_approximated
    vars = model.variables
    solvers = model.solvers

    nlsse(f,x) = model.nlsolve_static_function(f,x,param)
    sf(x) = model.static_function(x,param)
    df(x) = model.dynamic_function(x,param)

    ief = Array{Function}(undef,ne)
    for i = 1:ne
        ffie(x) = model.each_eqn_function[i](x,param)
        ief[i] = ffie
    end

    cf_cheb(state,scaled_weights,order,domain) = model.closure_function_chebyshev(state,scaled_weights,order,domain,param)
    cf_smol(state,scaled_weights,order,domain) = model.closure_function_smolyak(state,scaled_weights,order,domain,param)
    cf_hcross(state,scaled_weights,order,domain) = model.closure_function_hcross(state,scaled_weights,order,domain,param)
    cfpl_stoch(variables,grid,state,integrals) = model.closure_function_piecewise(variables,grid,state,integrals,param)
    cfpl_det(variables,grid,state) = model.closure_function_piecewise(variables,grid,state,param)

    if ns != 0
        newmod = REModel(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,vars,nlsse,sf,df,ief,cf_cheb,cf_smol,cf_hcross,cfpl_stoch,solvers)
        return newmod
    else
        newmod = REModel(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,vars,nlsse,sf,df,ief,cf_cheb,cf_smol,cf_hcross,cfpl_det,solvers)
        return newmod
    end

end

function assign_parameters(model::REModelPartial,paramdict::Dict{Q,T}) where {T<:Number,Q<:AbstractString}

    param = zeros(length(paramdict))

    for i in eachindex(model.unassigned_parameters)
        param[i] = paramdict[model.unassigned_parameters[i]]
    end

    newmod = assign_parameters(model,param)

    return newmod

end