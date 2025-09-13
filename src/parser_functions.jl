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

    # Remove sections that have been commented out via #= =#.  Does not allow multiple comment blocks on the same line.

    comment_block_start  = Int[]
    comment_block_finish = Int[]
    for i in eachindex(model_array)

        if occursin("#=",model_array[i]) == true
            push!(comment_block_start,i)
        end

        if occursin("=#",model_array[i]) == true
            push!(comment_block_finish,i)
        end

    end

    if length(comment_block_start) != length(comment_block_finish)
        error("Your comment blocks are formulated incorrectly")
    end

    if any(<(0), comment_block_finish-comment_block_start)
        error("Your comment blocks are formulated incorrectly")
    end

    for i in eachindex(comment_block_start)
        if comment_block_start[i] == comment_block_finish[i] # The comment block is contained within one line
            components_start  = split(model_array[comment_block_start[i]],"#=")
            components_finish = split(model_array[comment_block_finish[i]],"=#")
            if components_start[1] == ""
                model_array[comment_block_start[i]] = components_finish[2]
            elseif components_finish[2] == ""
                model_array[comment_block_start[i]] = components_start[1]
            else
                model_array[comment_block_start[i]] = prod(string.(components_start[1],components_finish[2]))
            end
        else # The comment block is spread over multiple lines.
            for j in comment_block_start[i]:comment_block_finish[i]
                if j == comment_block_start[i]
                    components = split(model_array[comment_block_start[i]],"#=")
                    if components[1] == ""
                        model_array[j] = ""
                    else
                        model_array[j] = string(components[1])
                    end
                elseif j == comment_block_finish[i]
                    components = split(model_array[comment_block_start[i]],"=#")
                    if length(components) ==1 || components[2] == "" 
                        model_array[j] = ""
                    else
                        model_array[j] = string(components[2])
                    end
                else
                    model_array[j] = ""
                end
            end
        end
    end

    # Remove lines that have been commented out via #. "The fat# hairy# cat." becomes "The fat."

    for i in eachindex(model_array)
        if occursin("#",model_array[i]) == true
            components = split(model_array[i],"#")
            model_array[i] = string(components[1])
        end
    end

    # Identify lines that contain only spaces

    for i in eachindex(model_array)
        if unique(model_array[i]) == [' '] # line contains only space characters
            model_array[i] = "" # Make it a blank line
        end
    end

    # Remove blank lines

    model_array = model_array[model_array .!= ""]

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
        error("The $term-designation appears multiple times in the model file, but should appear only once.")
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

    for end_location in startfrom:length(model_array)
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

    model_block_begins = find_term(model_array,term) + 1              # Starts the line after the term is located
    model_block_ends   = find_end(model_array,model_block_begins) - 1 # Ends the line before 'end' is located

    if model_block_begins > model_block_ends
        if term in ("shocks:","states:")
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    end

    model_block = model_array[model_block_begins:model_block_ends]

    # Remove any trailing variable separators: "," or ";".

    for i in eachindex(model_block)
        if endswith(model_block[i],union(",",";"))
            model_block[i] = model_block[i][1:end-1]
        end
    end

    # Extract the names and place them in a vector

    variables = String.(strip.(split(model_block[1],union(",",";"))))
    variables = [i for i in variables if i != ""]
    for i = 2:length(model_block)
        new_variables = String.(strip.(split(model_block[i],union(",",";"))))
        new_variables = [i for i in new_variables if i != ""]
        variables     = [variables; new_variables]
    end

    # Check whether names are repeated

    if length(variables) != length(unique(variables))
        error("Some $(term[1:end-1]) names are repreated.")
    end

    # Check to ensure that variables contains non-empty string elements

    if length(variables) == 1 && variables[1] == ""
        if term in ("shocks:", "states:")
            return String[]
        else
            error("The model file contains no $(term[1:end-1]).")
        end
    else
        return variables
    end

end

"""
Extracts the names and values for each of the model's parameters.  The parameter names are 
sorted so that larger names come first.

Internal function; not exposed to users.
"""
function get_parameters_and_values(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    parameter_block_begin = find_term(model_array,term) + 1                 # Starts the line after the term is located
    parameter_block_end   = find_end(model_array,parameter_block_begin) - 1 # Ends the line before 'end' is located
    if parameter_block_begin > parameter_block_end
        error("The model file contains no parameters.")
    end

    parameter_block = model_array[parameter_block_begin:parameter_block_end]

    # Remove any trailing separators: "," or ";".

    for i in eachindex(parameter_block)
        if endswith(parameter_block[i],union(",",";"))
            parameter_block[i] = parameter_block[i][1:end-1]
        end
    end

    revised_parameter_block = String.(strip.(split(parameter_block[1],union(",",";"))))
    revised_parameter_block = [i for i in revised_parameter_block if i != ""]
    for i = 2:length(parameter_block)
        new_revised_parameter_block = String.(strip.(split(parameter_block[i],union(",",";"))))
        new_revised_parameter_block = [i for i in new_revised_parameter_block if i != ""]
        revised_parameter_block     = [revised_parameter_block; new_revised_parameter_block]
    end

    # Extract the parameter names and values, making note of any parameters with unassigned values

    unassigned_parameter_index = 1
    unassigned_parameters      = Array{Q}(undef,0)

    parameters = Array{Q}(undef,length(revised_parameter_block))
    values = Array{Q}(undef,length(revised_parameter_block))
    for i in eachindex(revised_parameter_block)
        if occursin("=",revised_parameter_block[i]) == false # No value has been assigned to the parameter
            parameters[i] = revised_parameter_block[i]
            values[i]     = "p[$unassigned_parameter_index]" # p is a reserved name
            push!(unassigned_parameters,parameters[i])
            unassigned_parameter_index += 1
        else
            pair = strip.(split(revised_parameter_block[i],"="))
            if pair[2] == "" # No value has been assigned to the parameter
                parameters[i] = pair[1]
                values[i]     = "p[$unassigned_parameter_index]" # p is a reserved name
                push!(unassigned_parameters,parameters[i])
                unassigned_parameter_index += 1
            else
              parameters[i] = pair[1]
              values[i]     = pair[2]
            end
        end
    end

    # Check whether names are repeated

    if length(parameters) != length(unique(parameters))
        error("Some parameter names are repreated.")
    end

    return parameters, values, unassigned_parameters

end

"""
Extracts the model's equations.

Internal function; not exposed to users.
"""
function get_equations(model_array::Array{Q,1},term::Q) where {Q<:AbstractString}

    equation_block_begin = find_term(model_array,term) + 1                # Starts the line after the term is located
    equation_block_end   = find_end(model_array,equation_block_begin) - 1 # Ends the line before 'end' is located
    if equation_block_begin > equation_block_end
        error("The model file contains no equations.")
    end

    equation_block = model_array[equation_block_begin:equation_block_end]

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

    # For every model equation: 1) unify the bracketing; 2) check whether the open- close-parentheses are balanced.

    for i in eachindex(equation_block)

        if occursin("[",equations[i]) == true # Replace open square bracket with open round parenthesis
            equations[i] = replace(equations[i],"[" => "(")
        elseif occursin("]",equations[i]) == true # Replace close square bracket with close round parenthesis
            equations[i] = replace(equations[i],"]" => ")")
        end

        if occursin("=",equations[i]) == false # Check that each equation contains an equals sign
            error("Equation line $i does not contain an '=' sign.")
        end

        n_open_paren   = length(findall("(",equations[i]))
        n_closed_paren = length(findall(")",equations[i]))
        if n_open_paren != n_closed_paren # Check whether parentheses are balanced
            error("Equation line $i has $n_open_paren open parentheses and $n_closed_paren closed parentheses.")
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
        error("The specified solver is unknown.")
    end
end

# At this point all of the model file's contents has been read in, now we need to do some consistency checking on its contents.

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
        reordered_shocks = copy(shocks)
        return reordered_equations, reordered_states, reordered_shocks
    end

end

function catalogue_equations(equations,states,jumps,parameters) # Deterministic models

    junk_equations = copy(equations)

    if length(states) > 0
        combined_names = [states;jumps;parameters;["log", "exp", "deriv"]]
    else
        combined_names = [jumps;parameters;["log", "exp", "deriv"]]
    end

    combined_names .= combined_names[sortperm(length.(combined_names),rev=true)] # sort the names according to their length, longest to shortest

    n_eqns = length(equations) # number of equations

    # Initialise containers to store equation information

    states_info = [Array{String}(undef,0) for _ in 1:n_eqns]

    n_states_info = zeros(Int,n_eqns)
    n_jumps_info  = zeros(Int,n_eqns)

    n_future_states_info = zeros(Int,n_eqns)
    n_future_jumps_info  = zeros(Int,n_eqns)

    n_lag_states_info = zeros(Int,n_eqns)
    n_lag_jumps_info  = zeros(Int,n_eqns)

    # Begin equation analysis

    for name in combined_names
        for i in 1:n_eqns
            if occursin(name,junk_equations[i]) == true
                if name in parameters
                    junk_equations[i] = replace(junk_equations[i],name => ":")
                elseif name in states
                    if occursin("$name(-1)",junk_equations[i]) == true
                        n_lag_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(-1)" => ":")
                    end
                    if occursin("$name(+1)",junk_equations[i]) == true
                        n_future_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(+1)" => ":")
                    end
                    if occursin(name,junk_equations[i]) == true
                        n_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i], name => ":")
                        push!(states_info[i],name)
                    end
                elseif name in jumps
                    if occursin("$name(-1)",junk_equations[i]) == true
                        n_lag_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(-1)" => ":")
                    end
                    if occursin("$name(+1)",junk_equations[i]) == true
                        n_future_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(+1)" => ":")
                    end
                    if occursin(name,junk_equations[i]) == true
                        n_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i], name => ":")
                    end
                else
                    junk_equations[i] = replace(junk_equations[i], name => ":")
                end
            end
        end
    end

    return n_lag_states_info, n_states_info, n_future_states_info, n_lag_jumps_info, n_jumps_info, n_future_jumps_info, states_info

end

function catalogue_equations(equations,shocks,states,jumps,parameters) # Stochastic models

    junk_equations = copy(equations)

    if length(states) > 0
        combined_names = [shocks;states;jumps;parameters;["log", "exp", "deriv"]]
    else
        combined_names = [shocks;jumps;parameters;["log", "exp", "deriv"]]
    end

    combined_names .= combined_names[sortperm(length.(combined_names),rev=true)] # sort the names according to their length, longest to shortest

    n_eqns = length(equations) # number of equations

    # Initialise containers to store equation information

    shocks_info = [Array{String}(undef,0) for _ in 1:n_eqns]
    states_info = [Array{String}(undef,0) for _ in 1:n_eqns]

    n_shocks_info = zeros(Int,n_eqns)
    n_states_info = zeros(Int,n_eqns)
    n_jumps_info  = zeros(Int,n_eqns)

    n_future_states_info = zeros(Int,n_eqns)
    n_future_jumps_info  = zeros(Int,n_eqns)

    n_lag_states_info = zeros(Int,n_eqns)
    n_lag_jumps_info  = zeros(Int,n_eqns)

    # Begin equation analysis

    for name in combined_names
        for i in 1:n_eqns
            if occursin(name,junk_equations[i]) == true
                if name in parameters
                    junk_equations[i] = replace(junk_equations[i],name => ":")
                elseif name in states
                    if occursin("$name(-1)",junk_equations[i]) == true
                        n_lag_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(-1)" => ":")
                    end
                    if occursin("$name(+1)",junk_equations[i]) == true
                        n_future_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(+1)" => ":")
                    end
                    if occursin(name,junk_equations[i]) == true
                        n_states_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i], name => ":")
                        push!(states_info[i],name)
                    end
                elseif name in jumps
                    if occursin("$name(-1)",junk_equations[i]) == true
                        n_lag_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(-1)" => ":")
                    end
                    if occursin("$name(+1)",junk_equations[i]) == true
                        n_future_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i],"$name(+1)" => ":")
                    end
                    if occursin(name,junk_equations[i]) == true
                        n_jumps_info[i] += 1
                        junk_equations[i] = replace(junk_equations[i], name => ":")
                    end
                elseif name in shocks
                    n_shocks_info[i] += 1
                    junk_equations[i] = replace(junk_equations[i], name => ":")
                    push!(shocks_info[i],name)
                else
                    junk_equations[i] = replace(junk_equations[i], name => ":")
                end
            end
        end
    end

    return n_lag_states_info, n_states_info, n_future_states_info, n_lag_jumps_info, n_jumps_info, n_future_jumps_info, states_info, n_shocks_info, shocks_info

end

"""
Reorders the equations of a deteministic model so that the shock processes are at the top.

Internal function; not exposed to users.
"""
function _reorder_equations(equations::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},parameters::Array{Q,1}) where {Q<:AbstractString}

    n_eqns = length(equations)

    reordered_equations = copy(equations)
   
    # Retrieve summary information about each equation
   
    n_lag_states_info, n_states_info, n_future_states_info, n_lag_jumps_info, n_jumps_info, n_future_jumps_info, states_info = catalogue_equations(equations,states,jumps,parameters)

    # Put the equations for predetermined variables at the top

    jumps_locations = n_jumps_info + n_future_jumps_info
    ind = sortperm(jumps_locations)

    reordered_equations  .= reordered_equations[ind]
    n_lag_states_info    .= n_lag_states_info[ind]
    n_states_info        .= n_states_info[ind]
    n_future_states_info .= n_future_states_info[ind]
    n_lag_jumps_info     .= n_lag_jumps_info[ind]
    n_jumps_info         .= n_jumps_info[ind]
    n_future_jumps_info  .= n_future_jumps_info[ind]
    states_info          .= states_info[ind]

    # Now sort out the order of the states in the system.

    reordered_states = states_info[1]
    for i in 2:n_eqns
        if length(states_info[i]) > 0
            new_states = setdiff(states_info[i],reordered_states)
            reordered_states = [reordered_states;new_states]
        end
    end

    return reordered_equations, reordered_states

end

"""
Reorders the equations of a stochastic model so that the shock processes are at the top.

Internal function; not exposed to users.
"""
function _reorder_equations(equations::Array{Q,1},shocks::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1},parameters::Array{Q,1}) where {Q<:AbstractString}

    n_eqns = length(equations)

    reordered_equations = copy(equations)
   
    # Retrieve summary information about each equation
   
    n_lag_states_info, n_states_info, n_future_states_info, n_lag_jumps_info, n_jumps_info, n_future_jumps_info, states_info, n_shocks_info, shocks_info = catalogue_equations(equations,shocks,states,jumps,parameters)

    # Put equations containing shocks at the top

    shocks_locations = n_shocks_info .> 0
    ind = sortperm(shocks_locations,rev=true)

    reordered_equations  .= reordered_equations[ind]
    n_lag_states_info    .= n_lag_states_info[ind]
    n_states_info        .= n_states_info[ind]
    n_future_states_info .= n_future_states_info[ind]
    n_lag_jumps_info     .= n_lag_jumps_info[ind]
    n_jumps_info         .= n_jumps_info[ind]
    n_future_jumps_info  .= n_future_jumps_info[ind]
    states_info          .= states_info[ind]

    n_shocks_info        .= n_shocks_info[ind]
    shocks_info          .= shocks_info[ind]

    # Wityhin the block that contains shocks, reordered from smallest to largest number of shocks

    n_eqns_with_shocks = sum(shocks_locations)

    n_sub_shocks_info = n_shocks_info[1:n_eqns_with_shocks]
    sub_shocks_info   = shocks_info[1:n_eqns_with_shocks]
 
    ind = sortperm(n_sub_shocks_info)

    reordered_equations[1:n_eqns_with_shocks]  .= reordered_equations[1:n_eqns_with_shocks][ind]
    n_lag_states_info[1:n_eqns_with_shocks]    .= n_lag_states_info[1:n_eqns_with_shocks][ind]
    n_states_info[1:n_eqns_with_shocks]        .= n_states_info[1:n_eqns_with_shocks][ind]
    n_future_states_info[1:n_eqns_with_shocks] .= n_future_states_info[1:n_eqns_with_shocks][ind]
    n_lag_jumps_info[1:n_eqns_with_shocks]     .= n_lag_jumps_info[1:n_eqns_with_shocks][ind]
    n_jumps_info[1:n_eqns_with_shocks]         .= n_jumps_info[1:n_eqns_with_shocks][ind]
    n_future_jumps_info[1:n_eqns_with_shocks]  .= n_future_jumps_info[1:n_eqns_with_shocks][ind]
    states_info[1:n_eqns_with_shocks]          .= states_info[1:n_eqns_with_shocks][ind]

    n_shocks_info[1:n_eqns_with_shocks]        .= n_shocks_info[1:n_eqns_with_shocks][ind]
    shocks_info[1:n_eqns_with_shocks]          .= shocks_info[1:n_eqns_with_shocks][ind]

    # Sort out the order of the shocks in the system

    reordered_shocks = shocks_info[1]
    for i in 2:n_eqns
        if length(shocks_info[i]) > 0
            new_shocks = setdiff(shocks_info[i],reordered_shocks)
            reordered_shocks = [reordered_shocks;new_shocks]
        end
    end

    # Now sort out the order of the states in the system.

    reordered_states = states_info[1]
    for i in 2:n_eqns
        if length(states_info[i]) > 0
            new_states = setdiff(states_info[i],reordered_states)
            reordered_states = [reordered_states;new_states]
        end
    end

    return reordered_equations, reordered_states, reordered_shocks

end

"""
Replaces lagged variables with pseudo current variables, augmenting the state vector and 
the model's equations accordingly.

Internal function; not exposed to users.
"""
function deal_with_lags(equations::Array{Q,1},states::Array{Q,1},jumps::Array{Q,1}) where {Q<:AbstractString}

    if length(states) > 0
        variables = [states;jumps]
    else
        variables = [jumps;]
    end
    sorted_variables = variables[sortperm(length.(variables),rev=true)]

    n_variables = length(variables)

    lag_variables = string.(sorted_variables,"(-1)")
    
    reorganized_equations = copy(equations)

    # First we determine if any equation contains a lagged variable.

    model_has_lags = false
    for i in eachindex(reorganized_equations)
        for j in lag_variables
            if occursin(j,reorganized_equations[i]) == true
                model_has_lags = true
            end
        end
    end

    #= If a model contains lagged variables, then we introduce a pseudo variable
       in its place and augment the list of state variables and the set of model
       equations. =#

    if model_has_lags == false
        return equations, states, jumps, variables
    else
        new_states = String[]
        new_eqns   = String[]
        for j in 1:n_variables
            flag = false
            for i in eachindex(equations)
                if occursin(lag_variables[j],reorganized_equations[i]) == true
                    reorganized_equations[i] = replace(reorganized_equations[i],lag_variables[j] => string(sorted_variables[j],"lag"))
                    flag = true
                end
            end
            if flag == true
                push!(new_states,string(sorted_variables[j],"lag"))
                push!(new_eqns,string(sorted_variables[j],"lag(+1) = ",sorted_variables[j]))
            end
        end
    
        reorganized_states    = [states; new_states]
        reorganized_equations = [reorganized_equations; new_eqns]
        reorganized_variables = [reorganized_states; jumps]

        return reorganized_equations, reorganized_states, jumps, reorganized_variables

    end

end

"""
Takes the model-array read from the model file for a rational expectations model, 
extracts the critical model information, does some basic error checking, and 
returns the model's key information in a structure.

Internal function; not exposed to users.
"""
function get_re_model_primatives(model_array::Array{Q,1}) where {Q<:AbstractString}

    solvers   = get_solvers(model_array,"solvers:")

    states    = get_variables(model_array,"states:") # Repeated state names will error at this function call
    jumps     = get_variables(model_array,"jumps:")  # Repeated jump names will error at this function call
    shocks    = get_variables(model_array,"shocks:") # Repeated shock names will error at this function call; will be empty array for deterministic models

    equations = get_equations(model_array,"equations:")

    parameters, parametervalues, unassigned_parameters = get_parameters_and_values(model_array,"parameters:") # Repeated parameter names will error at this function call
    
    # Check that every state variable enters at least one equation

    for i in states
        if sum(occursin.(i,equations)) == false
            error("The state variable $i is not in any equation.")
        end
    end

    # Check that every jump variable enters at least one equation

    for i in jumps
        if sum(occursin.(i,equations)) == false
            error("The jump variable $i is not in any equation.")
        end
    end

    # Check that every parameter enters at least one equation or is used to define another parameter

    for i in parameters
        if sum(occursin.(i,equations)) + sum(occursin.(i,parametervalues)) == 0 # Recall that at this stage parametervalues can be convolutions of parameters
            println("Warning: The parameter $i is not in any equation.")
        end
    end

    # Check whether there are name overlaps across shocks, states, jumps, parameters

    if length(shocks) > 0
        combined_names = [parameters; states; jumps; shocks]
    else
        combined_names = [parameters; states; jumps]
    end
    if length(unique(combined_names)) != length(combined_names)
        error("Some parameters, states, jumps, or shocks have the same name.")
    end

    # Check whether names conflict with reserved names

    reserved_names = ("deriv", "exp", "log", "x", "p", ":", ";", "&")
    for name in reserved_names
        if name in combined_names
            error("$name is reserved and cannot be the name of a variable, a shock, or a parameter.")
        end
    end

    reordered_equations, reordered_states, reordered_shocks           = reorder_equations(equations,shocks,states,jumps,parameters)
    reorganized_equations, expanded_states, jumps, expanded_variables = deal_with_lags(reordered_equations,reordered_states,jumps)

    re_model_primatives = REModelPrimatives(expanded_states,jumps,reordered_shocks,expanded_variables,parameters,parametervalues,reorganized_equations,unassigned_parameters,solvers)

    return re_model_primatives

end

"""
Repackages the model's equations, replaces parameter names with parameter values.

Internal function; not exposed to users.
"""
function repackage_equations(model::DSGEModelPrimatives)

    equations       = model.equations
    shocks          = model.shocks
    variables       = model.variables
    parameters      = model.parameters
    parametervalues = model.parametervalues

    repackaged_equations       = copy(equations)
    repackaged_parametervalues = copy(parametervalues)

    if length(shocks) != 0
        combined_names = [variables; parameters; shocks; ["deriv", "exp", "log"]]
    else
        combined_names = [variables; parameters; ["deriv", "exp", "log"]]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]
    sorted_parameters     = parameters[sortperm(length.(parameters),rev = true)] # not needed --- parameters have already been sorted by length?

    #= First we go through every parameter expression and replace exp with : and log with ;.  
       This is to guard them during variable and parameter substitution. =#

    for i in eachindex(parametervalues)
        if occursin("exp",parametervalues[i]) == true
            repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],"exp" => ":")
        elseif occursin("log",parametervalues[i]) == true
            repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],"log" => ";")
        elseif occursin("deriv",parametervalues[i]) == true
            repackaged_parametervalues[i] = replace(repackaged_parametervalues[i],"deriv" => "&")
        end
    end

    #= Go through all parameters and deal with parameters depending on other parameters =#

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
            error("There is a circularity in the parameter definitions.")
        end
    end

    #= Now we go through every equation and replace future variables, variables, and
       shocks with a numbered element of a vector, "x".  We guard against overwriting 
       "deriv", "exp", and "log".  We also replace parameter names 
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
        elseif j == "deriv"
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => "&")
            end
        elseif j == "exp"
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => ":")
            end
        elseif j == "log"
            for i in eachindex(repackaged_equations)
                repackaged_equations[i] = replace(repackaged_equations[i],j => ";")
            end
        end
    end

    #= Finally, go back through every equation and restore exp, log, and deriv where necessary =#

    for i in eachindex(repackaged_equations)
        if occursin(":",repackaged_equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],":" => "exp")
        end
        if occursin(";",repackaged_equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],";" => "log")
        end
        if occursin("&",repackaged_equations[i]) == true
            repackaged_equations[i] = replace(repackaged_equations[i],"&" => "deriv")
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

    equations       = model.equations
    variables       = model.variables
    shocks          = model.shocks
    parameters      = model.parameters
    parametervalues = model.parametervalues

    steady_state_equations       = copy(equations)
    steady_state_parametervalues = copy(parametervalues)

    if length(shocks) != 0
        combined_names = [variables; parameters; shocks; ["deriv", "exp", "log"]]
    else
        combined_names = [variables; parameters; ["deriv", "exp", "log"]]
    end

    sorted_combined_names = combined_names[sortperm(length.(combined_names),rev = true)]
    sorted_parameters     = parameters[sortperm(length.(parameters),rev = true)]

    #= First we go through every parameter expression and replace exp with :, log with ;, and deriv with &.  
       This is to guard them during variables and parameter substitution. =#

    for i in eachindex(parametervalues)
        if occursin("exp",parametervalues[i]) == true
            steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],"exp" => ":")
        elseif occursin("log",parametervalues[i]) == true
            steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],"log" => ";")
        elseif occursin("deriv",parametervalues[i]) == true
            steady_state_parametervalues[i] = replace(steady_state_parametervalues[i],"deriv" => "&")
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
            error("There is a circularity in the parameter definitions.")
        end
    end

    # Now we go through every equation and replace future variables, variables, and
    # shocks with a numbered element of a vector, "x".  We guard against overwriting 
    # "deriv", "exp" and "log".  We also replace parameter names with parameter values.

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
        elseif j == "deriv"
            for i in eachindex(steady_state_equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => "&")
            end
        elseif j == "exp"
            for i in eachindex(steady_state_equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => ":")
            end
        elseif j == "log"
            for i in eachindex(steady_state_equations)
                steady_state_equations[i] = replace(steady_state_equations[i],j => ";")
            end
        end
    end

    #= Finally, go back through every equation and restore exp, log, and deriv where necessary =#

    for i in eachindex(steady_state_equations)
        if occursin(":",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],":" => "exp")
        end
        if occursin(";",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],";" => "log")
        end
        if occursin("&",steady_state_equations[i]) == true
            steady_state_equations[i] = replace(steady_state_equations[i],"&" => "deriv")
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
    eqns_with_jumps          = Int64[]
    for i = 1:ny
        for j = 1:ne
            if occursin("approx$i",projection_equations[j]) == true
                push!(jumps_to_be_approximated,i)
                push!(eqns_with_jumps,j)
            end
        end
    end

    jumps_to_be_approximated = unique(jumps_to_be_approximated)

    eqns_with_derivs              = Int64[]
    derivs_to_be_approximated_num = Int64[]
    derivs_to_be_approximated_den = Int64[]
    for j = 1:ne
        if occursin("deriv{",projection_equations[j]) == true
            push!(eqns_with_derivs,j)
            s = findall("{",projection_equations[j])
            f = findall("}",projection_equations[j])
            n_derivs = length(s)
            for jj = 1:n_derivs
                snippet = projection_equations[j][s[jj][1]:f[jj][1]]
                snippet1,snippet2 = split(snippet,"|")
                for i in 1:nv
                    if occursin("x[$(nx+i)]",snippet1) == true
                        push!(derivs_to_be_approximated_num,nx+i)
                        projection_equations[j] = replace(projection_equations[j],"deriv$snippet" => "deriv$(nx+i)")
                    end
                end
                for i in 1:nx
                    if occursin("state[$i]",snippet2) == true
                        push!(derivs_to_be_approximated_den,i)
                    end
                end
            end
        end
    end

    eqns_with_jumps = sort(unique(eqns_with_jumps))

    return projection_equations, jumps_to_be_approximated, eqns_with_jumps, derivs_to_be_approximated_num, derivs_to_be_approximated_den, eqns_with_derivs

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

    nonlinear_equations, jumps_to_be_approximated, eqns_with_jumps, derivs_to_be_approximated_num, derivs_to_be_approximated_den, eqns_with_derivs = create_projection_equations(repackaged_equations, model)
    
    projection_equations = make_equations_equal_zero(nonlinear_equations)

    steady_state_equations = create_steady_state_equations(model)
    static_equations = make_equations_equal_zero(steady_state_equations)
    dynamic_equations = make_equations_equal_zero(repackaged_equations)

    number_states = length(model.states)
    number_jumps = length(model.jumps)
    number_shocks = length(model.shocks)
    number_variables = length(model.variables)
    number_equations = length(model.equations)

    variables = OrderedDict(model.variables[i] => i for i = 1:number_variables)

    #unassigned_parameters = copy(model.unassigned_parameters)

    # Build up the string containing the processed model information that gets saved

    # First, add the model's summary information

    model_string = "nx = $number_states \n \n"
    model_string = string(model_string, "ny = $number_jumps \n \n")
    model_string = string(model_string, "ns = $number_shocks \n \n")
    model_string = string(model_string, "nv = $number_variables \n \n")
    model_string = string(model_string, "ne = $number_equations \n \n")

    model_string = string(model_string, "jumps_to_approximate = $jumps_to_be_approximated \n \n")
    model_string = string(model_string, "eqns_to_approximate = $eqns_with_jumps \n \n")
    model_string = string(model_string, "derivs_to_approximate_num = $derivs_to_be_approximated_num \n \n")
    model_string = string(model_string, "derivs_to_approximate_den = $derivs_to_be_approximated_den \n \n")
    model_string = string(model_string, "eqns_with_derivs = $eqns_with_derivs \n \n")
    model_string = string(model_string, "variables = $variables \n \n")
    model_string = string(model_string, "unassigned_parameters = $(model.unassigned_parameters) \n \n")
    model_string = string(model_string, """solvers = "$(model.solvers)" \n \n""")

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

    for i in 1:length(derivs_to_be_approximated_num)
       closure_cheb_string = string(closure_cheb_string, "    deriv$(derivs_to_be_approximated_num[i])  = chebyshev_derivative(scaled_weights[$weight_number],x[$number_jumps+1:end],$(derivs_to_be_approximated_den[i]),order,domain)", "\n")
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

    for i in 1:length(derivs_to_be_approximated_num)
        closure_smol_string = string(closure_smol_string, "    deriv$(derivs_to_be_approximated_num[i])  = smolyak_derivative(scaled_weights[$weight_number],x[$number_jumps+1:end],order,domain,$(derivs_to_be_approximated_den[i]))", "\n")
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

    for i in 1:length(derivs_to_be_approximated_num)
        closure_hcross_string = string(closure_hcross_string, "    deriv$(derivs_to_be_approximated_num[i])  = hyperbolic_cross_derivative(scaled_weights[$weight_number],x[$number_jumps+1:end],order,domain,$(derivs_to_be_approximated_den[i]))", "\n")
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
            closure_pl_string = "function closure_piecewise_equations(variables,grid,statearams) \n \n"
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

    for i in 1:length(derivs_to_be_approximated_num)
        if number_shocks == 0 || derivs_to_be_approximated_num[i] <= number_jumps
            closure_pl_string = string(closure_pl_string, "    deriv$(derivs_to_be_approximated_num[i])  = piecewise_linear_derivative(variables[$(derivs_to_be_approximated_num[i])],grid,x[$number_jumps+1:end],$(derivs_to_be_approximated_den[i]))", "\n")
        else
            closure_pl_string = string(closure_pl_string, "    deriv$(derivs_to_be_approximated_num[i])  = piecewise_linear_derivative(variables[$(derivs_to_be_approximated_num[i])],grid,x[$number_jumps+1:end],integrals,$(derivs_to_be_approximated_den[i]))", "\n")
        end
        weight_number += 1
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
    
    model_path = replace(path, ".txt" => "_processed.txt")
    open(model_path, "w") do io
        write(io, model_string)
    end

end

"""
Opens, and processes the contents of a model file, creating the processed model file.

Signature
=========
```
process_model(path)
```
"""
function process_model(path::Q) where {Q<:AbstractString}

    # Main function used to open, read, and process a model file.  The processed model
    # in written to a file that contains all the information needed for the model
    # solvers.

    model_array = open_model_file(path)

    # Creates the processed model structure for rational expectations models

    re_model_primatives = get_re_model_primatives(model_array)

    create_processed_model_file(re_model_primatives,path)

    println("The model's variables are now in this order: ",re_model_primatives.variables)
    println("The model's shocks are in this order:        ",re_model_primatives.shocks)
    if length(re_model_primatives.unassigned_parameters) != 0
        println("The following parameters do not have values assigned: $(re_model_primatives.unassigned_parameters)")
    end

end

"""
Retrives and stores in a model structure the information extracted from a processed model file.

Signature
=========
```
model = retrieve_processed_model()
```
"""
function retrieve_processed_model()

    global nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,derivs_to_approximate_num,derivs_to_approximate_den,eqns_with_derivs,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations,unassigned_parameters

    if length(unassigned_parameters) != 0
      if solvers == "Any"
        dsge_model = REModelPartialAny(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,derivs_to_approximate_num,derivs_to_approximate_den,eqns_with_derivs,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations,unassigned_parameters)
      elseif solvers == "Projection"
        dsge_model = REModelPartialProj(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,derivs_to_approximate_num,derivs_to_approximate_den,eqns_with_derivs,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations,unassigned_parameters)
      elseif solvers == "Perturbation"
        dsge_model = REModelPartialPert(nx,ny,ns,nv,ne,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,unassigned_parameters)
      elseif solvers == "Linear"
        dsge_model = REModelPartialLinear(nx,ny,ns,nv,ne,variables,nlsolve_static_equations,static_equations,dynamic_equations,unassigned_parameters)
      end
    else
      if solvers == "Any"
        dsge_model = REModelAny(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,derivs_to_approximate_num,derivs_to_approximate_den,eqns_with_derivs,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations)
      elseif solvers == "Projection"
        dsge_model = REModelProj(nx,ny,ns,nv,ne,jumps_to_approximate,eqns_to_approximate,derivs_to_approximate_num,derivs_to_approximate_den,eqns_with_derivs,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations,closure_chebyshev_equations,closure_smolyak_equations,closure_hcross_equations,closure_piecewise_equations)
      elseif solvers == "Perturbation"
        dsge_model = REModelPert(nx,ny,ns,nv,ne,variables,nlsolve_static_equations,static_equations,dynamic_equations,individual_equations)
      elseif solvers == "Linear"
        dsge_model = REModelLinear(nx,ny,ns,nv,ne,variables,nlsolve_static_equations,static_equations,dynamic_equations)
      end
    end

    return dsge_model

end

create_model_structure = retrieve_processed_model

"""
Assigns values to parameters in a partially specified model.  `param` can be 
either a vector or a dictionary.

Signature
=========
```
new_mod = assign_parameters(model,param)
```
"""
function assign_parameters(model::REModelPartialLinear,param::Array{T,1}) where {T<:Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    vars = model.variables

    nlsse(f,x) = model.nlsolve_static_function(f,x,param)
    sf(x) = model.static_function(x,param)
    df(x) = model.dynamic_function(x,param)

    newmod = REModelLinear(nx,ny,ns,nv,ne,vars,nlsse,sf,df)

    return newmod

end

function assign_parameters(model::REModelPartialPert,param::Array{T,1}) where {T<:Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    vars = model.variables

    nlsse(f,x) = model.nlsolve_static_function(f,x,param)
    sf(x) = model.static_function(x,param)
    df(x) = model.dynamic_function(x,param)

    ief = Array{Function}(undef,ne)
    for i = 1:ne
        ffie(x) = model.each_eqn_function[i](x,param)
        ief[i] = ffie
    end

    newmod = REModelPert(nx,ny,ns,nv,ne,vars,nlsse,sf,df,ief)

    return newmod

end

function assign_parameters(model::REModelPartialProj,param::Array{T,1}) where {T<:Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    jumps_approx = model.jumps_approximated
    eqns_approx = model.eqns_approximated
    vars = model.variables
    derivs_approx_num = model.derivs_approximated_num
    derivs_approx_den = model.derivs_approximated_den
    eqns_derivs = model.eqns_with_derivs

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
        newmod = REModelProj(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,derivs_approx_num,derivs_approx_den,eqns_derivs,vars,nlsse,sf,df,ief,cf_cheb,cf_smol,cf_hcross,cfpl_stoch)
        return newmod
    else
        newmod = REModelProj(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,derivs_approx_den,eqns_derivs,vars,nlsse,sf,df,cf_cheb,cf_smol,cf_hcross,cfpl_det)
        return newmod
    end

end

function assign_parameters(model::REModelPartialAny,param::Array{T,1}) where {T<:Number}

    nx = model.number_states
    ny = model.number_jumps
    ns = model.number_shocks
    nv = model.number_variables
    ne = model.number_equations
    jumps_approx = model.jumps_approximated
    eqns_approx = model.eqns_approximated
    vars = model.variables
    derivs_approx_num = model.derivs_approximated_num
    derivs_approx_den = model.derivs_approximated_den
    eqns_derivs = model.eqns_with_derivs

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
        newmod = REModelAny(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,derivs_approx_num,derivs_approx_den,eqns_derivs,vars,nlsse,sf,df,ief,cf_cheb,cf_smol,cf_hcross,cfpl_stoch)
        return newmod
    else
        newmod = REModelAny(nx,ny,ns,nv,ne,jumps_approx,eqns_approx,derivs_approx_num,derivs_approx_den,eqns_derivs,vars,nlsse,sf,df,ief,cf_cheb,cf_smol,cf_hcross,cfpl_det)
        return newmod
    end

end

function assign_parameters(model::M,paramdict::Dict{Q,T}) where {M<:REModelPartial,T<:Number,Q<:AbstractString}

    param = zeros(length(paramdict))

    for i in eachindex(model.unassigned_parameters)
        param[i] = paramdict[model.unassigned_parameters[i]]
    end

    newmod = assign_parameters(model,param)

    return newmod

end