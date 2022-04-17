##################### Analysis functions #########################

function decision_rule(soln::R) where {R<:Union{FirstOrderSolutionDet,FirstOrderSolutionStoch}}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end

        y = soln.gbar + soln.gx*(state - soln.hbar)

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:SecondOrderSolutionDet}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.gxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:SecondOrderSolutionStoch}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*soln.gss[i] + (1/2)*(soln.gxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:ThirdOrderSolutionDet}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.gxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:ThirdOrderSolutionStoch}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*soln.gss[i] + (1/2)*(soln.gxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*soln.gsss[i] + (3/6)*(soln.gssx[i:i,:]*(state-soln.hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:FourthOrderSolutionDet}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.gxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)))[1] + (1/24)*(soln.gxxxx[i:i,:]*kron(kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)),(state-soln.hbar)))[1]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:FourthOrderSolutionStoch}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)
        ny = length(soln.gbar)

        y = zeros(ny)

        for i = 1:ny
            y[i] = soln.gbar[i] + (soln.gx[i:i,:]*(state-soln.hbar))[1] + (1/2)*soln.gss[i] + (1/2)*(soln.gxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (3/6)*(soln.gssx[i:i,:]*(state-soln.hbar))[1] + (1/6)*(soln.gxxx[i:i,:]*kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)))[1] + (6/24)*(soln.gssxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (1/24)*(soln.gxxxx[i:i,:]*kron(kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)),(state-soln.hbar)))[1] + (1/24)*soln.gssss[i]
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:Union{ChebyshevSolutionDet,ChebyshevSolutionStoch}}

    nx = length(soln.nodes)
    nv = length(soln.variables)
    ny = nv - nx

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,ny)
    if soln.node_generator == chebyshev_nodes
        for i = 1:ny
            w[i] = chebyshev_weights(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:ny
            w[i] = chebyshev_weights_extrema(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:ny
            w[i] = chebyshev_weights_extended(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        y = zeros(ny)

        for i = 1:ny
            y[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:Union{SmolyakSolutionDet,SmolyakSolutionStoch}}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,ny)
    smol_iim = smolyak_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:ny
        weights[i] = smolyak_weights(soln.variables[nx+i],smol_iim)
    end

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        y    = zeros(ny)
        poly = smolyak_polynomial(state,soln.multi_index,soln.domain)
        for i = 1:ny
            y[i] = smolyak_evaluate(weights[i],poly)
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:Union{HyperbolicCrossSolutionDet,HyperbolicCrossSolutionStoch}}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx

    T = eltype(soln.variables[1])

    weights    = Array{Array{T,1},1}(undef,ny)
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)
    for i = 1:ny
        weights[i] = hyperbolic_cross_weights(soln.variables[nx+i],hcross_iim)
    end

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        y    = zeros(ny)
        poly = hyperbolic_cross_polynomial(state,soln.multi_index,soln.domain)

        for i = 1:ny
            y[i] = hyperbolic_cross_evaluate(weights[i],poly)
        end

        return y

    end

    return create_decision_rule

end

function decision_rule(soln::R) where {R<:Union{PiecewiseLinearSolutionDet,PiecewiseLinearSolutionStoch}}

    function create_decision_rule(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.nodes)
            error("state vector has incorrect size")
        end
        nx = length(soln.nodes)
        nv = length(soln.variables)
        ny = nv - nx

        y = zeros(ny)

        for i = 1:ny
            y[i] = piecewise_linear_evaluate(soln.variables[nx+i],soln.nodes,state)
        end

        return y

    end

    return create_decision_rule

end

function state_transition(soln::R) where {R<:FirstOrderSolutionDet}

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end

        x_update = soln.hbar + soln.hx*(state - soln.hbar)

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:FirstOrderSolutionStoch}

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end

        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = soln.hbar + soln.hx*(state - soln.hbar) + soln.k*shocks

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:SecondOrderSolutionDet}

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.hxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:SecondOrderSolutionStoch}

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (soln.k[i:i,:]*shocks)[1] + (1/2)*soln.hss[i] + (1/2)*(soln.hxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:ThirdOrderSolutionDet}

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.hxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:ThirdOrderSolutionStoch}

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (soln.k[i:i,:]*shocks)[1] + (1/2)*soln.hss[i] + (1/2)*(soln.hxx[i:i,:]*kron((state - soln.hbar),(state - soln.hbar)))[1] + (1/6)*soln.hsss[i] + (3/6)*(soln.hssx[i:i,:]*(state-soln.hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state - soln.hbar),(state - soln.hbar)),(state - soln.hbar)))[1]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:FourthOrderSolutionDet}

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (1/2)*(soln.hxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)))[1] + (1/24)*(soln.hxxxx[i:i,:]*kron(kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)),(state-soln.hbar)))[1]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:FourthOrderSolutionStoch}

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != length(soln.hbar)
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end
        nx = length(soln.hbar)

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = soln.hbar[i] + (soln.hx[i:i,:]*(state-soln.hbar))[1] + (soln.k[i:i,:]*shocks)[1] + (1/2)*soln.hss[i] + (1/2)*(soln.hxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (3/6)*(soln.hssx[i:i,:]*(state-soln.hbar))[1] + (1/6)*(soln.hxxx[i:i,:]*kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)))[1] + (6/24)*(soln.hssxx[i:i,:]*kron((state-soln.hbar),(state-soln.hbar)))[1] + (1/24)*(soln.hxxxx[i:i,:]*kron(kron(kron((state-soln.hbar),(state-soln.hbar)),(state-soln.hbar)),(state-soln.hbar)))[1] + (1/24)*soln.hssss[i]
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:ChebyshevSolutionDet}

    nx = length(soln.nodes)

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,nx)
    if soln.node_generator == chebyshev_nodes
        for i = 1:nx
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nx
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nx
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:ChebyshevSolutionStoch}

    nx = length(soln.nodes)
    ns = size(soln.k,2)

    T = eltype(soln.variables[1])

    w = Array{Array{T,nx},1}(undef,nx)
    if soln.node_generator == chebyshev_nodes
        for i = 1:nx
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nx
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nx
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = chebyshev_evaluate(w[i],state,soln.order,soln.domain)
        end
        x_update[1:ns] .+= soln.k*shocks

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:SmolyakSolutionDet}

    nx = size(soln.grid,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    smol_iim = smolyak_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nx
        weights[i] = smolyak_weights(soln.variables[i],smol_iim)
    end

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)
        poly     = smolyak_polynomial(state,soln.multi_index,soln.domain)

        for i = 1:nx
            x_update[i] = smolyak_evaluate(weights[i],poly)
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:SmolyakSolutionStoch}

    nx = size(soln.grid,2)
    ns = size(soln.k,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    smol_iim = smolyak_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nx
        weights[i] = smolyak_weights(soln.variables[i],smol_iim)
    end

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)
        poly     = smolyak_polynomial(state,soln.multi_index,soln.domain)

        for i = 1:nx
            x_update[i] = smolyak_evaluate(weights[i],poly)
        end
        x_update[1:ns] .+= soln.k*shocks

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:HyperbolicCrossSolutionDet}

    nx = size(soln.grid,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nx
        weights[i] = hyperbolic_cross_weights(soln.variables[i],hcross_iim)
    end

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)
        poly = hyperbolic_cross_polynomial(state,soln.multi_index,soln.domain)

        for i = 1:nx
            x_update[i] = hyperbolic_cross_evaluate(weights[i],poly)
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:HyperbolicCrossSolutionStoch}

    nx = size(soln.grid,2)
    ns = size(soln.k,2)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,nx)
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nx
        weights[i] = hyperbolic_cross_weights(soln.variables[i],hcross_iim)
    end

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)
        poly = hyperbolic_cross_polynomial(state,soln.multi_index,soln.domain)

        for i = 1:nx
            x_update[i] = hyperbolic_cross_evaluate(weights[i],poly)
        end
        x_update[1:ns] .+= soln.k*shocks

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:PiecewiseLinearSolutionDet}

    nx = length(soln.nodes)

    function create_state_transition(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = piecewise_linear_evaluate(soln.variables[i],soln.nodes,state)
        end

        return x_update

    end

    return create_state_transition

end

function state_transition(soln::R) where {R<:PiecewiseLinearSolutionStoch}

    nx = length(soln.nodes)
    ns = size(soln.k,2)

    function create_state_transition(state::AbstractArray{T,1},shocks::AbstractArray{T,1}) where {T<:AbstractFloat}

        if length(state) != nx
            error("state vector has incorrect size")
        end
        if length(shocks) != size(soln.k,2)
            error("shocks vector has incorrect size")
        end

        x_update = zeros(nx)

        for i = 1:nx
            x_update[i] = piecewise_linear_evaluate(soln.variables[i],soln.nodes,state)
        end
        x_update[1:ns] .+= soln.k*shocks

        return x_update

    end

    return create_state_transition

end

function expected_jumps(soln::R) where {R<:Union{PerturbationSolutionDet,ProjectionSolutionDet}}

    dec_rules = decision_rule(soln)
    state_trans = state_transition(soln)

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = dec_rules(state_trans(state))

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R<:PerturbationSolutionStoch}

    dec_rules = decision_rule(soln)
    state_trans = state_transition(soln)

    ns = size(soln.k,2)

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = dec_rules(state_trans(state,zeros(ns)))

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R<:ChebyshevSolutionStoch}

    nx = length(soln.nodes)
    nv = length(soln.variables)
    ny = nv - nx
    ns = size(soln.k,2)

    state_trans = state_transition(soln)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,nx},1}(undef,ny)
    scaled_weights = Array{Array{T,nx},1}(undef,ny)
    if soln.node_generator == chebyshev_nodes
        for i = 1:ny
            weights[i] = chebyshev_weights(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:ny
            weights[i] = chebyshev_weights_extrema(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:ny
            weights[i] = chebyshev_weights_extended(soln.variables[nx+i],soln.nodes,soln.order,soln.domain)
        end
    end

    if typeof(soln.integrals) == Array{T,ns}
        for i = 1:ny
            for j = 1:ns
                scaled_weights[i] = soln.integrals .* weights[i]
            end
        end
    else
        for i = 1:ny
            for j = 1:ns
                index = [1:ndims(weights[i]);]
                index[1],index[j] = index[j],index[1]
                scaled_weights[i] = permutedims(soln.integrals[j] .* permutedims(weights[i],index),index)
            end
        end
    end

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = zeros(ny)

        for i = 1:ny

            y[i] = chebyshev_evaluate(scaled_weights[i],state_trans(state,zeros(ns)),soln.order,soln.domain)

        end

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R<:SmolyakSolutionStoch}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx
    ns = size(soln.k,2)

    state_trans = state_transition(soln)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,ny)
    for i = 1:ny
        weights[i] = smolyak_weights(soln.variables[nx+i],soln.grid,soln.multi_index,soln.domain) .* soln.scale_factor
    end

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = zeros(ny)

        for i = 1:ny

            y[i] = smolyak_evaluate(weights[i],state_trans(state,zeros(ns)),soln.multi_index,soln.domain)

        end

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R<:HyperbolicCrossSolutionStoch}

    nx = size(soln.grid,2)
    nv = length(soln.variables)
    ny = nv - nx
    ns = size(soln.k,2)

    state_trans = state_transition(soln)

    T = eltype(soln.variables[1])

    weights = Array{Array{T,1},1}(undef,ny)
    for i = 1:ny
        weights[i] = hyperbolic_cross_weights(soln.variables[nx+i],soln.grid,soln.multi_index,soln.domain) .* soln.scale_factor
    end

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = zeros(ny)

        for i = 1:ny

            y[i] = hyperbolic_cross_evaluate(weights[i],state_trans(state,zeros(ns)),soln.multi_index,soln.domain)

        end

        return y

    end

    return create_expected_jumps

end

function expected_jumps(soln::R) where {R<:PiecewiseLinearSolutionStoch}

    dec_rules = decision_rule(soln)
    state_trans = state_transition(soln)

    ns = size(soln.k,2)

    function create_expected_jumps(state::AbstractArray{T,1}) where {T<:AbstractFloat}

        y = dec_rules(state_trans(state,zeros(ns)))

        return y

    end

    return create_expected_jumps

end

function state_space_eqm(soln::R) where {R<:ModelSolution}

    dec_rules = decision_rule(soln)
    transition_eqn = state_transition(soln)
    forecast_eqn = expected_jumps(soln)

    eqm_dynamics = StateSpaceEqm(dec_rules,transition_eqn,forecast_eqn)

    return eqm_dynamics

end

function compute_mean(soln::R) where {R<:PerturbationSolution}

    if typeof(soln) <: FirstOrderSolutionStoch

        mean_states = zeros(length(soln.hbar))
        mean_jumps  = zeros(length(soln.gbar))

    elseif typeof(soln) <: SecondOrderSolutionStoch

        nx = length(soln.hbar)
        ny = length(soln.gbar)

        mean_states = (I-soln.hx)\((1/2)*soln.hss + (1/2)*hxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k')))
        mean_jumps  = soln.gx*mean_states + (1/2)*soln.gss + (1/2)*gxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k'))

    elseif typeof(soln) <: ThirdOrderSolutionStoch

        mean_states = (I-soln.hx)\((1/2)*soln.hss + (1/2)*soln.hxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k')))
        mean_jumps  = soln.gx*mean_states + (1/2)*soln.gss + (1/2)*soln.gxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k'))

        term1 = (I-kron(kron(soln.hx,soln.hx),soln.hx))\(kron(kron(soln.k,soln.k),soln.k)*vec(soln.skewness))
        term2 = (I-kron(soln.hx,soln.hx))\(kron(soln.hx,(1/2)*soln.hxx)*term1)

        mean_states += (I-soln.hx)\(soln.hxx*term2 + (1/6)*soln.hxxx*term1 + (1/2)*soln.hssx + (1/6)*soln.hsss)
        mean_jumps  += soln.gx*((I-soln.hx)\(soln.hxx*term2 + (1/6)*soln.hxxx*term1 + (1/6)*soln.hsss)) + soln.gxx*term2 + (1/6)*soln.gxxx*term1 + (1/2)*soln.gssx + (1/6)*soln.gsss

    elseif typeof(soln) <: FourthOrderSolutionStoch

        mean_states = (I-soln.hx)\((1/2)*soln.hss + (1/2)*soln.hxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k')))
        mean_jumps  = soln.gx*mean_states + (1/2)*soln.gss + (1/2)*soln.gxx*((I-kron(soln.hx,soln.hx))\vec(soln.k*soln.sigma*soln.k'))

        mean_states += (I-soln.hx)\((1/2)*soln.hssx + (1/2)*soln.hssxx + (1/24)*soln.hssss)
        mean_jumps  += soln.gx*((I-soln.hx)\((1/2)*soln.hssx + (1/2)*soln.hssxx + (1/24)*soln.hssss)) + (1/2)*soln.gssx + (1/2)*soln.gssxx + (1/24)*soln.gssss

    else # All deterministic cases are handled here

        mean_states = zeros(length(soln.hbar))
        mean_jumps = zeros(length(soln.gbar))

    end

    return mean_states .+ soln.hbar,mean_jumps .+ soln.gbar

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R<:FirstOrderSolutionDet,T<:Real,S<:Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length] .+ soln.hbar; simulated_jumps_f[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S; rndseed = 123456) where {R<:FirstOrderSolutionStoch,T<:Real,S<:Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.sigma,2)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(ns)
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
    end

    return [simulated_states_f[:,1:sim_length] .+ soln.hbar; simulated_jumps_f[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R<:SecondOrderSolutionDet,T<:Real,S<:Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_s = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i] = soln.hx*simulated_states_s[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
    end

    simulated_states = simulated_states_f + simulated_states_s
    simulated_jumps = simulated_jumps_f + simulated_jumps_s

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S; rndseed = 123456) where {R<:SecondOrderSolutionStoch,T<:Real,S<:Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.sigma,2)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss
    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state - soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_s = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(ns)
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i] = soln.hx*simulated_states_s[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + hss
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + gss
    end

    simulated_states = simulated_states_f + simulated_states_s
    simulated_jumps = simulated_jumps_f + simulated_jumps_s

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R<:ThirdOrderSolutionDet,T<:Real,S<:Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hxxx = (1/6)*soln.hxxx
    gxxx = (1/6)*soln.gxxx

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_s = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    simulated_states_t = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_t = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i] = soln.hx*simulated_states_s[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_states_t[:,i] = soln.hx*simulated_states_t[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
        simulated_jumps_t[:,i-1] = soln.gx*simulated_states_t[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t
    simulated_jumps = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S; rndseed = 123456) where {R<:ThirdOrderSolutionStoch,T<:Real,S<:Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.sigma,2)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss
    hxxx = (1/6)*soln.hxxx
    gxxx = (1/6)*soln.gxxx
    hssx = (3/6)*soln.hssx
    gssx = (3/6)*soln.gssx
    hsss = (1/6)*soln.hsss
    gsss = (1/6)*soln.gsss

    simulated_states_f = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_f = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_s = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    simulated_states_t = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps_t = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i] = soln.hx*simulated_states_f[:,i-1] + soln.k*randn(ns)
        simulated_jumps_f[:,i-1] = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i] = soln.hx*simulated_states_s[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + hss
        simulated_jumps_s[:,i-1] = soln.gx*simulated_states_s[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + gss
        simulated_states_t[:,i] = soln.hx*simulated_states_t[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + hssx*simulated_states_f[:,i-1] + hsss
        simulated_jumps_t[:,i-1] = soln.gx*simulated_states_t[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + gssx*simulated_states_f[:,i-1] + gsss
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t
    simulated_jumps = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R<:FourthOrderSolutionDet,T<:Real,S<:Integer}

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hxx   = (1/2)*soln.hxx
    gxx   = (1/2)*soln.gxx
    hxxx  = (1/6)*soln.hxxx
    gxxx  = (1/6)*soln.gxxx
    hxxxx = (1/24)*soln.hxxxx
    gxxxx = (1/24)*soln.gxxxx

    simulated_states_f = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f  = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s  = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    simulated_states_t = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_t  = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] .= zeros(nx)

    simulated_states_fo = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_fo  = Array{T,2}(undef,ny,sim_length)
    simulated_states_fo[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i]   = soln.hx*simulated_states_f[:,i-1]
        simulated_jumps_f[:,i-1]  = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]   = soln.hx*simulated_states_s[:,i-1]  + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_jumps_s[:,i-1]  = soln.gx*simulated_states_s[:,i-1]  + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1])
        simulated_states_t[:,i]   = soln.hx*simulated_states_t[:,i-1]  + hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
        simulated_jumps_t[:,i-1]  = soln.gx*simulated_states_t[:,i-1]  + gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
        simulated_states_fo[:,i]  = soln.hx*simulated_states_fo[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_t[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_s[:,i-1]) + hxxxx*kron(kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
        simulated_jumps_fo[:,i-1] = soln.gx*simulated_states_fo[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_t[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_s[:,i-1]) + gxxxx*kron(kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]),simulated_states_f[:,i-1])
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t + simulated_states_fo
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t + simulated_jumps_fo

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S; rndseed = 123456) where {R<:FourthOrderSolutionStoch,T<:Real,S<:Integer}

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.sigma,2)

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    hxx   = (1/2)*soln.hxx
    gxx   = (1/2)*soln.gxx
    hss   = (1/2)*soln.hss
    gss   = (1/2)*soln.gss
    hxxx  = (1/6)*soln.hxxx
    gxxx  = (1/6)*soln.gxxx
    hssx  = (3/6)*soln.hssx
    gssx  = (3/6)*soln.gssx
    hxxxx = (1/24)*soln.hxxxx
    gxxxx = (1/24)*soln.gxxxx
    hssxx = (6/24)*soln.hssxx
    gssxx = (6/24)*soln.gssxx
    hssss = (1/24)*soln.hssss
    gssss = (1/24)*soln.gssss

    simulated_states_f = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_f  = Array{T,2}(undef,ny,sim_length)
    simulated_states_f[:,1] .= initial_state .- soln.hbar

    simulated_states_s = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_s  = Array{T,2}(undef,ny,sim_length)
    simulated_states_s[:,1] .= zeros(nx)

    simulated_states_t = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_t  = Array{T,2}(undef,ny,sim_length)
    simulated_states_t[:,1] .= zeros(nx)

    simulated_states_fo = Array{T,2}(undef,nx,sim_length+1)
    simulated_jumps_fo  = Array{T,2}(undef,ny,sim_length)
    simulated_states_fo[:,1] .= zeros(nx)

    @views for i = 2:sim_length+1
        simulated_states_f[:,i]   = soln.hx*simulated_states_f[:,i-1]  + soln.k*randn(ns)
        simulated_jumps_f[:,i-1]  = soln.gx*simulated_states_f[:,i-1]
        simulated_states_s[:,i]   = soln.hx*simulated_states_s[:,i-1]  + hxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + hss
        simulated_jumps_s[:,i-1]  = soln.gx*simulated_states_s[:,i-1]  + gxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + gss
        simulated_states_t[:,i]   = soln.hx*simulated_states_t[:,i-1]  + hxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + hssx*simulated_states_f[:,i-1] 
        simulated_jumps_t[:,i-1]  = soln.gx*simulated_states_t[:,i-1]  + gxx*kron(simulated_states_f[:,i-1],simulated_states_s[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + gssx*simulated_states_f[:,i-1] 
        simulated_states_fo[:,i]  = soln.hx*simulated_states_fo[:,i-1] + hxx*kron(simulated_states_f[:,i-1],simulated_states_t[:,i-1]) + hxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_s[:,i-1]) + hssx*simulated_states_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]),simulated_states_f[:,i-1])  + hssxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + hssss
        simulated_jumps_fo[:,i-1] = soln.gx*simulated_states_fo[:,i-1] + gxx*kron(simulated_states_f[:,i-1],simulated_states_t[:,i-1]) + gxxx*kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_s[:,i-1]) + gssx*simulated_states_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]),simulated_states_f[:,i-1]),simulated_states_f[:,i-1]) + gssxx*kron(simulated_states_f[:,i-1],simulated_states_f[:,i-1]) + gssss
    end

    simulated_states = simulated_states_f + simulated_states_s + simulated_states_t + simulated_states_fo
    simulated_jumps  = simulated_jumps_f + simulated_jumps_s + simulated_jumps_t + simulated_jumps_fo

    return [simulated_states[:,1:sim_length] .+ soln.hbar; simulated_jumps[:,1:end] .+ soln.gbar]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S) where {R<:ProjectionSolutionDet,T<:Real,S<:Integer}

    eqm = state_space_eqm(soln)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] .= initial_state

    for i = 2:sim_length+1
        simulated_states[:,i] .= eqm.h(simulated_states[:,i-1])
        simulated_jumps[:,i-1] .= eqm.g(simulated_states[:,i-1])
    end

    return [simulated_states[:,1:sim_length]; simulated_jumps[:,1:end]]

end

function simulate(soln::R,initial_state::Array{T,1},sim_length::S; rndseed = 123456) where {R<:ProjectionSolutionStoch,T<:Real,S<:Integer}

    Random.seed!(rndseed)

    eqm = state_space_eqm(soln)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(initial_state) != nx
        error("The number of inital values for the states must equal the number of states")
    end

    simulated_states = Array{T,2}(undef,nx,sim_length + 1)
    simulated_jumps = Array{T,2}(undef,ny,sim_length)
    simulated_states[:,1] .= initial_state

    for i = 2:sim_length+1
        simulated_states[:,i] .= eqm.h(simulated_states[:,i-1],randn(ns))
        simulated_jumps[:,i-1] .= eqm.g(simulated_states[:,i-1])
    end

    return [simulated_states[:,1:sim_length]; simulated_jumps[:,1:end]]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},rndseed = 123456) where {R<:FirstOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f = zeros(nx,n + 1)
    simulated_jumps_pos_f = zeros(ny,n)
    simulated_states_pos_f[:,1] .= soln.k*innovation_vector

    simulated_states_neg_f = zeros(nx,n + 1)
    simulated_jumps_neg_f = zeros(ny,n)
    simulated_states_neg_f[:,1] .= -soln.k*innovation_vector

    @views for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n]; simulated_jumps_pos_f],[simulated_states_neg_f[:,1:n]; simulated_jumps_neg_f]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},rndseed = 123456) where {R<:FirstOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f = zeros(nx,n + 1)
    simulated_jumps_pos_f = zeros(ny,n)
    simulated_states_pos_f[:,1] .= initial_state - soln.hbar + soln.k*innovation_vector

    simulated_states_neg_f = zeros(nx,n + 1)
    simulated_jumps_neg_f = zeros(ny,n)
    simulated_states_neg_f[:,1] .= initial_state - soln.hbar - soln.k*innovation_vector

    @views for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n]+soln.hbar; simulated_jumps_pos_f+soln.gbar],[simulated_states_neg_f[:,1:n]+soln.hbar; simulated_jumps_neg_f+soln.gbar]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:FirstOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f = zeros(nx,n + 1)
    simulated_jumps_pos_f = zeros(ny,n)
    simulated_states_pos_f[:,1] .= soln.k*innovation_vector

    simulated_states_neg_f = zeros(nx,n + 1)
    simulated_jumps_neg_f = zeros(ny,n)
    simulated_states_neg_f[:,1] .= -soln.k*innovation_vector

    @views for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n]; simulated_jumps_pos_f],[simulated_states_neg_f[:,1:n]; simulated_jumps_neg_f]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:FirstOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    simulated_states_pos_f = zeros(nx,n + 1)
    simulated_jumps_pos_f = zeros(ny,n)
    simulated_states_pos_f[:,1] .= initial_state - soln.hbar + soln.k*innovation_vector

    simulated_states_neg_f = zeros(nx,n + 1)
    simulated_jumps_neg_f = zeros(ny,n)
    simulated_states_neg_f[:,1] .= initial_state - soln.hbar - soln.k*innovation_vector

    @views for i = 2:n+1
        simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1]
        simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
        simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1]
        simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
    end

    return [simulated_states_pos_f[:,1:n]; simulated_jumps_pos_f],[simulated_states_neg_f[:,1:n]; simulated_jumps_neg_f]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:SecondOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss

    sample = simulate(soln,soln.hss,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n + 1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n + 1)
        simulated_jumps_pos_s  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n + 1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n + 1)
        simulated_jumps_neg_s  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n + 1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n + 1)
        simulated_jumps_base_s  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)] - soln.hbar
        simulated_states_pos_f[:,1]  .= initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f + simulated_jumps_pos_s - simulated_jumps_base_f - simulated_jumps_base_s)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f + simulated_jumps_neg_s - simulated_jumps_base_f - simulated_jumps_base_s)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:SecondOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n + 1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n + 1)
        simulated_jumps_pos_s  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n + 1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n + 1)
        simulated_jumps_neg_s  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n + 1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n + 1)
        simulated_jumps_base_s  = zeros(ny,n)

        simulated_states_pos_f[:,1]  .= initial_state - soln.hbar + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.hbar - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state - soln.hbar

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f + simulated_jumps_pos_s - simulated_jumps_base_f - simulated_jumps_base_s)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s - simulated_states_base_f - simulated_states_base_s)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f + simulated_jumps_neg_s - simulated_jumps_base_f - simulated_jumps_base_s)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:ThirdOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hxxx = (1/6)*soln.hxxx
    gxxx = (1/6)*soln.gxxx
    hssx = (3/6)*soln.hssx
    gssx = (3/6)*soln.gssx
    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss
    hsss = (1/6)*soln.hsss
    gsss = (1/6)*soln.gsss

    sample = simulate(soln,soln.hss,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n + 1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n + 1)
        simulated_jumps_pos_s  = zeros(ny,n)
        simulated_states_pos_t = zeros(nx,n + 1)
        simulated_jumps_pos_t  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n + 1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n + 1)
        simulated_jumps_neg_s  = zeros(ny,n)
        simulated_states_neg_t = zeros(nx,n + 1)
        simulated_jumps_neg_t  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n + 1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n + 1)
        simulated_jumps_base_s  = zeros(ny,n)
        simulated_states_base_t = zeros(nx,n + 1)
        simulated_jumps_base_t  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)] - soln.hbar
        simulated_states_pos_f[:,1]  .= initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss
            simulated_states_pos_t[:,i]  = soln.hx*simulated_states_pos_t[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] + hsss
            simulated_jumps_pos_t[:,i-1] = soln.gx*simulated_states_pos_t[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] + gsss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss
            simulated_states_neg_t[:,i]  = soln.hx*simulated_states_neg_t[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] + hsss
            simulated_jumps_neg_t[:,i-1] = soln.gx*simulated_states_neg_t[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] + gsss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
            simulated_states_base_t[:,i]  = soln.hx*simulated_states_base_t[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssx*simulated_states_base_f[:,i-1] + hsss
            simulated_jumps_base_t[:,i-1] = soln.gx*simulated_states_base_t[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssx*simulated_states_base_f[:,i-1] + gsss
        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s + simulated_states_pos_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f + simulated_jumps_pos_s + simulated_jumps_pos_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s + simulated_states_neg_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f + simulated_jumps_neg_s + simulated_jumps_neg_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:ThirdOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx = (1/2)*soln.hxx
    gxx = (1/2)*soln.gxx
    hxxx = (1/6)*soln.hxxx
    gxxx = (1/6)*soln.gxxx
    hssx = (3/6)*soln.hssx
    gssx = (3/6)*soln.gssx
    hss = (1/2)*soln.hss
    gss = (1/2)*soln.gss
    hsss = (1/6)*soln.hsss
    gsss = (1/6)*soln.gsss

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f = zeros(nx,n + 1)
        simulated_jumps_pos_f  = zeros(ny,n)
        simulated_states_pos_s = zeros(nx,n + 1)
        simulated_jumps_pos_s  = zeros(ny,n)
        simulated_states_pos_t = zeros(nx,n + 1)
        simulated_jumps_pos_t  = zeros(ny,n)

        simulated_states_neg_f = zeros(nx,n + 1)
        simulated_jumps_neg_f  = zeros(ny,n)
        simulated_states_neg_s = zeros(nx,n + 1)
        simulated_jumps_neg_s  = zeros(ny,n)
        simulated_states_neg_t = zeros(nx,n + 1)
        simulated_jumps_neg_t  = zeros(ny,n)

        simulated_states_base_f = zeros(nx,n + 1)
        simulated_jumps_base_f  = zeros(ny,n)
        simulated_states_base_s = zeros(nx,n + 1)
        simulated_jumps_base_s  = zeros(ny,n)
        simulated_states_base_t = zeros(nx,n + 1)
        simulated_jumps_base_t  = zeros(ny,n)

        simulated_states_pos_f[:,1]  .= initial_state - soln.hbar + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.hbar - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state - soln.hbar

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            simulated_states_pos_f[:,i]  = soln.hx*simulated_states_pos_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1] = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]  = soln.hx*simulated_states_pos_s[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1] = soln.gx*simulated_states_pos_s[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss
            simulated_states_pos_t[:,i]  = soln.hx*simulated_states_pos_t[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] + hsss
            simulated_jumps_pos_t[:,i-1] = soln.gx*simulated_states_pos_t[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] + gsss

            simulated_states_neg_f[:,i]  = soln.hx*simulated_states_neg_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1] = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]  = soln.hx*simulated_states_neg_s[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1] = soln.gx*simulated_states_neg_s[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss
            simulated_states_neg_t[:,i]  = soln.hx*simulated_states_neg_t[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] + hsss
            simulated_jumps_neg_t[:,i-1] = soln.gx*simulated_states_neg_t[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] + gsss

            simulated_states_base_f[:,i]  = soln.hx*simulated_states_base_f[:,i-1] + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1] = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]  = soln.hx*simulated_states_base_s[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1] = soln.gx*simulated_states_base_s[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
            simulated_states_base_t[:,i]  = soln.hx*simulated_states_base_t[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssx*simulated_states_base_f[:,i-1] + hsss
            simulated_jumps_base_t[:,i-1] = soln.gx*simulated_states_base_t[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssx*simulated_states_base_f[:,i-1] + gsss
        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s + simulated_states_pos_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f + simulated_jumps_pos_s + simulated_jumps_pos_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s + simulated_states_neg_t - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f + simulated_jumps_neg_s + simulated_jumps_neg_t - simulated_jumps_base_f - simulated_jumps_base_s - simulated_jumps_base_t)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:FourthOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx   = (1/2)*soln.hxx
    gxx   = (1/2)*soln.gxx
    hss   = (1/2)*soln.hss
    gss   = (1/2)*soln.gss
    hxxx  = (1/6)*soln.hxxx
    gxxx  = (1/6)*soln.gxxx
    hssx  = (3/6)*soln.hssx
    gssx  = (3/6)*soln.gssx
    hxxxx = (1/24)*soln.hxxxx
    gxxxx = (1/24)*soln.gxxxx
    hssxx = (6/24)*soln.hssxx
    gssxx = (6/24)*soln.gssxx
    hssss = (1/24)*soln.hssss
    gssss = (1/24)*soln.gssss

    sample = simulate(soln,soln.hss,5*reps+100)

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f  = zeros(nx,n+1)
        simulated_jumps_pos_f   = zeros(ny,n)
        simulated_states_pos_s  = zeros(nx,n+1)
        simulated_jumps_pos_s   = zeros(ny,n)
        simulated_states_pos_t  = zeros(nx,n+1)
        simulated_jumps_pos_t   = zeros(ny,n)
        simulated_states_pos_fo = zeros(nx,n+1)
        simulated_jumps_pos_fo  = zeros(ny,n)

        simulated_states_neg_f  = zeros(nx,n+1)
        simulated_jumps_neg_f   = zeros(ny,n)
        simulated_states_neg_s  = zeros(nx,n+1)
        simulated_jumps_neg_s   = zeros(ny,n)
        simulated_states_neg_t  = zeros(nx,n+1)
        simulated_jumps_neg_t   = zeros(ny,n)
        simulated_states_neg_fo = zeros(nx,n+1)
        simulated_jumps_neg_fo  = zeros(ny,n)

        simulated_states_base_f  = zeros(nx,n+1)
        simulated_jumps_base_f   = zeros(ny,n)
        simulated_states_base_s  = zeros(nx,n+1)
        simulated_jumps_base_s   = zeros(ny,n)
        simulated_states_base_t  = zeros(nx,n+1)
        simulated_jumps_base_t   = zeros(ny,n)
        simulated_states_base_fo = zeros(nx,n+1)
        simulated_jumps_base_fo  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)] - soln.hbar
        simulated_states_pos_f[:,1]  .= initial_state + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state

        innovations = randn(size(soln.k,2),n+1)

        @views for i = 2:n+1

            simulated_states_pos_f[:,i]   = soln.hx*simulated_states_pos_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1]  = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]   = soln.hx*simulated_states_pos_s[:,i-1]  + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1]  = soln.gx*simulated_states_pos_s[:,i-1]  + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss
            simulated_states_pos_t[:,i]   = soln.hx*simulated_states_pos_t[:,i-1]  + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] 
            simulated_jumps_pos_t[:,i-1]  = soln.gx*simulated_states_pos_t[:,i-1]  + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] 
            simulated_states_pos_fo[:,i]  = soln.hx*simulated_states_pos_fo[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_t[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_s[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hssss
            simulated_jumps_pos_fo[:,i-1] = soln.gx*simulated_states_pos_fo[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_t[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_s[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gssss
    
            simulated_states_neg_f[:,i]   = soln.hx*simulated_states_neg_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1]  = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]   = soln.hx*simulated_states_neg_s[:,i-1]  + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1]  = soln.gx*simulated_states_neg_s[:,i-1]  + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss
            simulated_states_neg_t[:,i]   = soln.hx*simulated_states_neg_t[:,i-1]  + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] 
            simulated_jumps_neg_t[:,i-1]  = soln.gx*simulated_states_neg_t[:,i-1]  + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] 
            simulated_states_neg_fo[:,i]  = soln.hx*simulated_states_neg_fo[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_t[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_s[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hssss
            simulated_jumps_neg_fo[:,i-1] = soln.gx*simulated_states_neg_fo[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_t[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_s[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gssss

            simulated_states_base_f[:,i]   = soln.hx*simulated_states_base_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1]  = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]   = soln.hx*simulated_states_base_s[:,i-1]  + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1]  = soln.gx*simulated_states_base_s[:,i-1]  + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
            simulated_states_base_t[:,i]   = soln.hx*simulated_states_base_t[:,i-1]  + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_s[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssx*simulated_states_base_f[:,i-1] 
            simulated_jumps_base_t[:,i-1]  = soln.gx*simulated_states_base_t[:,i-1]  + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_s[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssx*simulated_states_base_f[:,i-1] 
            simulated_states_base_fo[:,i]  = soln.hx*simulated_states_base_fo[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_t[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_s[:,i-1]) + hssx*simulated_states_base_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hssss
            simulated_jumps_base_fo[:,i-1] = soln.gx*simulated_states_base_fo[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_t[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_s[:,i-1]) + gssx*simulated_states_base_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gssss

        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s + simulated_states_pos_t + simulated_states_pos_fo - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t - simulated_states_base_fo)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f  + simulated_jumps_pos_s  + simulated_jumps_pos_t  + simulated_jumps_pos_fo  - simulated_jumps_base_f  - simulated_jumps_base_s  - simulated_jumps_base_t  - simulated_jumps_base_fo)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s + simulated_states_neg_t + simulated_states_neg_fo - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t - simulated_states_base_fo)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f  + simulated_jumps_neg_s  + simulated_jumps_neg_t  + simulated_jumps_neg_fo  - simulated_jumps_base_f  - simulated_jumps_base_s  - simulated_jumps_base_t  - simulated_jumps_base_fo)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S;rndseed = 123456) where {R<:FourthOrderSolutionStoch,S<:Integer,T<:Real}

    if length(innovation_vector) > size(soln.k,2)
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < size(soln.k,2)
        error("Each shock needs an innovation (even if it's zero).")
    end

    Random.seed!(rndseed)

    nx = length(soln.hbar)
    ny = length(soln.gbar)

    hxx   = (1/2)*soln.hxx
    gxx   = (1/2)*soln.gxx
    hss   = (1/2)*soln.hss
    gss   = (1/2)*soln.gss
    hxxx  = (1/6)*soln.hxxx
    gxxx  = (1/6)*soln.gxxx
    hssx  = (3/6)*soln.hssx
    gssx  = (3/6)*soln.gssx
    hxxxx = (1/24)*soln.hxxxx
    gxxxx = (1/24)*soln.gxxxx
    hssxx = (6/24)*soln.hssxx
    gssxx = (6/24)*soln.gssxx
    hssss = (1/24)*soln.hssss
    gssss = (1/24)*soln.gssss

    impulses_states_pos = zeros(nx,n+1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n+1)
    impulses_jumps_neg  = zeros(ny,n)

    for j = 1:reps
        simulated_states_pos_f  = zeros(nx,n+1)
        simulated_jumps_pos_f   = zeros(ny,n)
        simulated_states_pos_s  = zeros(nx,n+1)
        simulated_jumps_pos_s   = zeros(ny,n)
        simulated_states_pos_t  = zeros(nx,n+1)
        simulated_jumps_pos_t   = zeros(ny,n)
        simulated_states_pos_fo = zeros(nx,n+1)
        simulated_jumps_pos_fo  = zeros(ny,n)

        simulated_states_neg_f  = zeros(nx,n+1)
        simulated_jumps_neg_f   = zeros(ny,n)
        simulated_states_neg_s  = zeros(nx,n+1)
        simulated_jumps_neg_s   = zeros(ny,n)
        simulated_states_neg_t  = zeros(nx,n+1)
        simulated_jumps_neg_t   = zeros(ny,n)
        simulated_states_neg_fo = zeros(nx,n+1)
        simulated_jumps_neg_fo  = zeros(ny,n)

        simulated_states_base_f  = zeros(nx,n+1)
        simulated_jumps_base_f   = zeros(ny,n)
        simulated_states_base_s  = zeros(nx,n+1)
        simulated_jumps_base_s   = zeros(ny,n)
        simulated_states_base_t  = zeros(nx,n+1)
        simulated_jumps_base_t   = zeros(ny,n)
        simulated_states_base_fo = zeros(nx,n+1)
        simulated_jumps_base_fo  = zeros(ny,n)

        simulated_states_pos_f[:,1]  .= initial_state - soln.hbar + soln.k*innovation_vector
        simulated_states_neg_f[:,1]  .= initial_state - soln.hbar - soln.k*innovation_vector
        simulated_states_base_f[:,1] .= initial_state - soln.hbar

        innovations = randn(size(soln.k,2),n+1)

        @views for i = 2:n+1

            simulated_states_pos_f[:,i]   = soln.hx*simulated_states_pos_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_pos_f[:,i-1]  = soln.gx*simulated_states_pos_f[:,i-1]
            simulated_states_pos_s[:,i]   = soln.hx*simulated_states_pos_s[:,i-1]  + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hss
            simulated_jumps_pos_s[:,i-1]  = soln.gx*simulated_states_pos_s[:,i-1]  + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gss
            simulated_states_pos_t[:,i]   = soln.hx*simulated_states_pos_t[:,i-1]  + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] 
            simulated_jumps_pos_t[:,i-1]  = soln.gx*simulated_states_pos_t[:,i-1]  + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_s[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] 
            simulated_states_pos_fo[:,i]  = soln.hx*simulated_states_pos_fo[:,i-1] + hxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_t[:,i-1]) + hxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_s[:,i-1]) + hssx*simulated_states_pos_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + hssxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + hssss
            simulated_jumps_pos_fo[:,i-1] = soln.gx*simulated_states_pos_fo[:,i-1] + gxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_t[:,i-1]) + gxxx*kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_s[:,i-1]) + gssx*simulated_states_pos_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]),simulated_states_pos_f[:,i-1]) + gssxx*kron(simulated_states_pos_f[:,i-1],simulated_states_pos_f[:,i-1]) + gssss
    
            simulated_states_neg_f[:,i]   = soln.hx*simulated_states_neg_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_neg_f[:,i-1]  = soln.gx*simulated_states_neg_f[:,i-1]
            simulated_states_neg_s[:,i]   = soln.hx*simulated_states_neg_s[:,i-1]  + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hss
            simulated_jumps_neg_s[:,i-1]  = soln.gx*simulated_states_neg_s[:,i-1]  + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gss
            simulated_states_neg_t[:,i]   = soln.hx*simulated_states_neg_t[:,i-1]  + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] 
            simulated_jumps_neg_t[:,i-1]  = soln.gx*simulated_states_neg_t[:,i-1]  + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_s[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] 
            simulated_states_neg_fo[:,i]  = soln.hx*simulated_states_neg_fo[:,i-1] + hxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_t[:,i-1]) + hxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_s[:,i-1]) + hssx*simulated_states_neg_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + hssxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + hssss
            simulated_jumps_neg_fo[:,i-1] = soln.gx*simulated_states_neg_fo[:,i-1] + gxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_t[:,i-1]) + gxxx*kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_s[:,i-1]) + gssx*simulated_states_neg_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]),simulated_states_neg_f[:,i-1]) + gssxx*kron(simulated_states_neg_f[:,i-1],simulated_states_neg_f[:,i-1]) + gssss

            simulated_states_base_f[:,i]   = soln.hx*simulated_states_base_f[:,i-1]  + soln.k*innovations[:,i]
            simulated_jumps_base_f[:,i-1]  = soln.gx*simulated_states_base_f[:,i-1]
            simulated_states_base_s[:,i]   = soln.hx*simulated_states_base_s[:,i-1]  + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hss
            simulated_jumps_base_s[:,i-1]  = soln.gx*simulated_states_base_s[:,i-1]  + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gss
            simulated_states_base_t[:,i]   = soln.hx*simulated_states_base_t[:,i-1]  + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_s[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssx*simulated_states_base_f[:,i-1] 
            simulated_jumps_base_t[:,i-1]  = soln.gx*simulated_states_base_t[:,i-1]  + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_s[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssx*simulated_states_base_f[:,i-1] 
            simulated_states_base_fo[:,i]  = soln.hx*simulated_states_base_fo[:,i-1] + hxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_t[:,i-1]) + hxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_s[:,i-1]) + hssx*simulated_states_base_f[:,i-1] + hxxxx*kron(kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + hssxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + hssss
            simulated_jumps_base_fo[:,i-1] = soln.gx*simulated_states_base_fo[:,i-1] + gxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_t[:,i-1]) + gxxx*kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_s[:,i-1]) + gssx*simulated_states_base_f[:,i-1] + gxxxx*kron(kron(kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]),simulated_states_base_f[:,i-1]) + gssxx*kron(simulated_states_base_f[:,i-1],simulated_states_base_f[:,i-1]) + gssss

        end

        impulses_states_pos .+= (simulated_states_pos_f + simulated_states_pos_s + simulated_states_pos_t + simulated_states_pos_fo - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t - simulated_states_base_fo)
        impulses_jumps_pos  .+= (simulated_jumps_pos_f  + simulated_jumps_pos_s  + simulated_jumps_pos_t  + simulated_jumps_pos_fo  - simulated_jumps_base_f  - simulated_jumps_base_s  - simulated_jumps_base_t  - simulated_jumps_base_fo)
        impulses_states_neg .+= (simulated_states_neg_f + simulated_states_neg_s + simulated_states_neg_t + simulated_states_neg_fo - simulated_states_base_f - simulated_states_base_s - simulated_states_base_t - simulated_states_base_fo)
        impulses_jumps_neg  .+= (simulated_jumps_neg_f  + simulated_jumps_neg_s  + simulated_jumps_neg_t  + simulated_jumps_neg_fo  - simulated_jumps_base_f  - simulated_jumps_base_s  - simulated_jumps_base_t  - simulated_jumps_base_fo)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:ChebyshevSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{eltype(soln.domain),N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1] .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1] .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1] .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = chebyshev_evaluate(w[j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_states_neg[j,i]  = chebyshev_evaluate(w[j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_states_base[j,i] = chebyshev_evaluate(w[j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_jumps_neg[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_jumps_base[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:ChebyshevSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    N = ndims(soln.weights[1])

    w = Array{Array{eltype(soln.domain),N},1}(undef,length(soln.variables))
    if soln.node_generator == chebyshev_nodes
        for i = 1:nv
            w[i] = chebyshev_weights(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extrema
        for i = 1:nv
            w[i] = chebyshev_weights_extrema(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    elseif soln.node_generator == chebyshev_extended
        for i = 1:nv
            w[i] = chebyshev_weights_extended(soln.variables[i],soln.nodes,soln.order,soln.domain)
        end
    end

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base = zeros(ny,n)

        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = chebyshev_evaluate(w[j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_states_neg[j,i]  = chebyshev_evaluate(w[j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_states_base[j,i] = chebyshev_evaluate(w[j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_pos[:,i-1],soln.order,soln.domain)
                simulated_jumps_neg[j,i-1]  = chebyshev_evaluate(w[nx+j],simulated_states_neg[:,i-1],soln.order,soln.domain)
                simulated_jumps_base[j,i-1] = chebyshev_evaluate(w[nx+j],simulated_states_base[:,i-1],soln.order,soln.domain)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:SmolyakSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    smol_iim = smolyak_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],smol_iim)
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            poly_pos  = smolyak_polynomial(simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
            poly_neg  = smolyak_polynomial(simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
            poly_base = smolyak_polynomial(simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            for j = 1:nx
                simulated_states_pos[j,i]  = smolyak_evaluate(w[j],poly_pos)
                simulated_states_neg[j,i]  = smolyak_evaluate(w[j],poly_neg)
                simulated_states_base[j,i] = smolyak_evaluate(w[j],poly_base)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = smolyak_evaluate(w[nx+j],poly_pos)
                simulated_jumps_neg[j,i-1]  = smolyak_evaluate(w[nx+j],poly_neg)
                simulated_jumps_base[j,i-1] = smolyak_evaluate(w[nx+j],poly_base)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:SmolyakSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    smol_iim = smolyak_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)
    
    for i = 1:nv
        w[i] = smolyak_weights(soln.variables[i],soln.grid,smol_iim)
    end

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base  = zeros(ny,n)

        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            poly_pos  = smolyak_polynomial(simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
            poly_neg  = smolyak_polynomial(simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
            poly_base = smolyak_polynomial(simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            for j = 1:nx
                simulated_states_pos[j,i]  = smolyak_evaluate(w[j],poly_pos)
                simulated_states_neg[j,i]  = smolyak_evaluate(w[j],poly_neg)
                simulated_states_base[j,i] = smolyak_evaluate(w[j],poly_base)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = smolyak_evaluate(w[nx+j],poly_pos)
                simulated_jumps_neg[j,i-1]  = smolyak_evaluate(w[nx+j],poly_neg)
                simulated_jumps_base[j,i-1] = smolyak_evaluate(w[nx+j],poly_base)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:HyperbolicCrossSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nv
        w[i] = hyperbolic_cross_weights(soln.variables[i],hcross_iim)
    end

    estimated_steady_state = zeros(nx)
    for i = 1:nx
        estimated_steady_state[i] = soln.variables[i][i]
    end

    sample = simulate(soln,estimated_steady_state,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base  = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            poly_pos  = hyperbolic_cross_polynomial(simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
            poly_neg  = hyperbolic_cross_polynomial(simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
            poly_base = hyperbolic_cross_polynomial(simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            for j = 1:nx
                simulated_states_pos[j,i]  = hyperbolic_cross_evaluate(w[j],poly_pos)
                simulated_states_neg[j,i]  = hyperbolic_cross_evaluate(w[j],poly_neg)
                simulated_states_base[j,i] = hyperbolic_cross_evaluate(w[j],poly_base)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = hyperbolic_cross_evaluate(w[nx+j],poly_pos)
                simulated_jumps_neg[j,i-1]  = hyperbolic_cross_evaluate(w[nx+j],poly_neg)
                simulated_jumps_base[j,i-1] = hyperbolic_cross_evaluate(w[nx+j],poly_base)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:HyperbolicCrossSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    w = Array{Array{eltype(soln.domain),1},1}(undef,length(soln.variables))
    hcross_iim = hyperbolic_cross_inverse_interpolation_matrix(soln.grid,soln.multi_index,soln.domain)

    for i = 1:nv
        w[i] = hyperbolic_cross_weights(soln.variables[i],hcross_iim)
    end

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos  = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg  = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos  = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg  = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base  = zeros(ny,n)

        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            poly_pos  = hyperbolic_cross_polynomial(simulated_states_pos[:,i-1],soln.multi_index,soln.domain)
            poly_neg  = hyperbolic_cross_polynomial(simulated_states_neg[:,i-1],soln.multi_index,soln.domain)
            poly_base = hyperbolic_cross_polynomial(simulated_states_base[:,i-1],soln.multi_index,soln.domain)
            for j = 1:nx
                simulated_states_pos[j,i]  = hyperbolic_cross_evaluate(w[j],poly_pos)
                simulated_states_neg[j,i]  = hyperbolic_cross_evaluate(w[j],poly_neg)
                simulated_states_base[j,i] = hyperbolic_cross_evaluate(w[j],poly_base)
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = hyperbolic_cross_evaluate(w[nx+j],poly_pos)
                simulated_jumps_neg[j,i-1]  = hyperbolic_cross_evaluate(w[nx+j],poly_neg)
                simulated_jumps_base[j,i-1] = hyperbolic_cross_evaluate(w[nx+j],poly_base)
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:PiecewiseLinearSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    estimated_steady_state = vec((soln.domain[1,:] + soln.domain[2,:]))/2

    sample = simulate(soln,estimated_steady_state,5*reps + 100)

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base = zeros(ny,n)

        initial_state = sample[1:nx,rand(101:5*reps+100)]
        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_states_neg[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_states_base[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_base[:,i-1])
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_jumps_neg[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_jumps_base[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_base[:,i-1])
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function impulses(soln::R,n::S,initial_state::Array{T,1},innovation_vector::Array{T,1},reps::S; rndseed = 123456) where {R<:PiecewiseLinearSolutionStoch,S<:Integer,T<:Real}

    Random.seed!(rndseed)

    nv = length(soln.variables)
    nx = size(soln.domain,2)
    ns = size(soln.k,2)
    ny = nv - nx

    if length(innovation_vector) > ns
        error("There are more innovations than shocks.")
    elseif length(innovation_vector) < ns
        error("Each shock needs an innovation (even if it's zero).")
    end

    impulses_states_pos = zeros(nx,n + 1)
    impulses_jumps_pos = zeros(ny,n)
    impulses_states_neg = zeros(nx,n + 1)
    impulses_jumps_neg = zeros(ny,n)

    for l = 1:reps
        simulated_states_pos = zeros(nx,n + 1)
        simulated_jumps_pos = zeros(ny,n)

        simulated_states_neg = zeros(nx,n + 1)
        simulated_jumps_neg = zeros(ny,n)

        simulated_states_base = zeros(nx,n + 1)
        simulated_jumps_base = zeros(ny,n)

        simulated_states_pos[:,1]    .= initial_state
        simulated_states_pos[1:ns,1] .+= soln.k*innovation_vector
        simulated_states_neg[:,1]    .= initial_state
        simulated_states_neg[1:ns,1] .-= soln.k*innovation_vector
        simulated_states_base[:,1]   .= initial_state

        innovations = randn(size(soln.k,2),n + 1)

        @views for i = 2:n+1
            for j = 1:nx
                simulated_states_pos[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_states_neg[j,i]  = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_states_base[j,i] = piecewise_linear_evaluate(soln.variables[j],soln.nodes,simulated_states_base[:,i-1])
            end
            simulated_states_pos[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_neg[1:ns,i]  .+= soln.k*innovations[:,i]
            simulated_states_base[1:ns,i] .+= soln.k*innovations[:,i]
            for j = 1:ny
                simulated_jumps_pos[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_pos[:,i-1])
                simulated_jumps_neg[j,i-1]  = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_neg[:,i-1])
                simulated_jumps_base[j,i-1] = piecewise_linear_evaluate(soln.variables[nx+j],soln.nodes,simulated_states_base[:,i-1])
            end
        end

        impulses_states_pos .+= (simulated_states_pos - simulated_states_base)
        impulses_jumps_pos  .+= (simulated_jumps_pos - simulated_jumps_base)
        impulses_states_neg .+= (simulated_states_neg - simulated_states_base)
        impulses_jumps_neg  .+= (simulated_jumps_neg - simulated_jumps_base)

    end

    impulses_states_pos = impulses_states_pos/reps
    impulses_jumps_pos  = impulses_jumps_pos/reps
    impulses_states_neg = impulses_states_neg/reps
    impulses_jumps_neg  = impulses_jumps_neg/reps

    return [impulses_states_pos[:,1:n]; impulses_jumps_pos],[impulses_states_neg[:,1:n]; impulses_jumps_neg]

end

function approximate_density(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:Real,S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    if point < a || point > b
        error("'point' must be between 'a' and 'b'")
    end

    n = 0

    c = zeros(order + 1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b - a))*cos((i - 1)*pi*(sample[j] - a)/(b - a))
            end
        end
    end
    c = c/n

    f = c[1]/2.0
    for i = 2:order+1
        f += c[i]*cos((i - 1)*pi*(point - a)/(b - a))
    end

    return f

end

function approximate_density(sample::Array{T,1},order::S,a::T,b::T) where {T<:Real,S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    n = 0

    c = zeros(order + 1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b - a))*cos((i - 1)*pi*(sample[j] - a)/(b - a))
            end
        end
    end
    c = c/n

    points = range(a,b,length = minimum([round(Int,n/100),100]))
    ff = zeros(length(points))
    for j in eachindex(ff)
        f = c[1]/2.0
        for i = 2:order+1
            f += c[i]*cos((i - 1)*pi*(points[j] - a)/(b - a))
        end
        ff[j] = f
    end

    return collect(points),ff

end

function approximate_distribution(sample::Array{T,1},point::T,order::S,a::T,b::T) where {T<:Real,S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    if point < a || point > b
        error("'point' must be between 'a' and 'b'")
    end

    n = 0

    c = zeros(order + 1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b - a))*cos((i - 1)*pi*(sample[j] - a)/(b - a))
            end
        end
    end
    c = c/n

    F = (point - a)/(b - a)
    for i = 2:order+1
        F += (c[i]*(b - a)/(pi*(i - 1)))*sin((i - 1)*pi*(point - a)/(b - a))
    end

    return F

end

function approximate_distribution(sample::Array{T,1},order::S,a::T,b::T) where {T<:Real,S<:Integer}

    if a >= b
        error("'a' must be less than 'b'")
    end

    n = 0

    c = zeros(order + 1)
    for j in eachindex(sample)
        if sample[j] >= a && sample[j] <= b
            n += 1
            for i = 1:order+1
                c[i] += (2.0/(b - a))*cos((i - 1)*pi*(sample[j] - a)/(b - a))
            end
        end
    end
    c = c/n

    points = range(a,b,length = minimum([round(Int,n/100),100]))
    FF = zeros(length(points))
    for j in eachindex(FF)
        F = (points[j] - a)/(b - a)
        for i = 2:order+1
            F += (c[i]*(b - a)/(pi*(i - 1)))*sin((i - 1)*pi*(points[j] - a)/(b - a))
        end
        FF[j] = F
    end

    return collect(points),FF

end

function compare_solutions(solna::R1,solnb::R2,domain::Array{T,2},n::S,seed::S = 123456) where {T<:Real,S<:Integer,R1<:ModelSolution,R2<:ModelSolution}

    if typeof(solna) <: PerturbationSolution && typeof(solnb) <: PerturbationSolution
        if length(solna.gbar) != length(solnb.gbar) || length(solna.hbar) != length(solnb.hbar)
            error("The solutions are not of the same model.")
        end
        nx = length(solna.hbar)
        ny = length(solna.gbar)
        sup_errors = zeros(ny)
    elseif typeof(solna) <: PerturbationSolution && typeof(solnb) <: ProjectionSolution
        if length(solna.hbar) != size(solnb.domain,2) || length(solna.gbar) != length(solnb.variables) - size(solnb.domain,2)
            error("The solutions are not of the same model.")
        end
        nx = length(solna.hbar)
        ny = length(solna.gbar)
        sup_errors = zeros(ny)
    elseif typeof(solna) <: ProjectionSolution && typeof(solnb) <: PerturbationSolution
        if length(solnb.hbar) != size(solna.domain,2) || length(solnb.gbar) != length(solna.variables) - size(solna.domain,2)
            error("The solutions are not of the same model.")
        end
        nx = length(solnb.hbar)
        ny = length(solnb.gbar)
        sup_errors = zeros(ny)
    else
        if size(solna.domain) != size(solnb.domain) || length(solna.variables) != length(solnb.variables)
            error("The solutions are not of the same model.")
        end
        nx = size(solna.domain,2)
        ny = length(solna.variables) - nx
        sup_errors = zeros(ny)
    end

    eqma = state_space_eqm(solna)
    eqmb = state_space_eqm(solnb)

    Random.seed!(seed)

    state = domain[2,:] .+ rand(nx,n) .* (domain[1,:] - domain[2,:])
    varsa = zeros(ny,n)
    varsb = zeros(ny,n)
    for i = 1:n
        varsa[:,i] .= eqma.g(state[:,i])
        varsb[:,i] .= eqmb.g(state[:,i])
    end
    for j = 1:ny
        sup_errors[j] = maximum(abs,varsa[j,:] - varsb[j,:])
    end

    return sup_errors

end

function euler_errors(model::REModel,soln::R,domain::Union{Array{T,2},Array{T,1}},npoints::S,seed::S = 123456) where {S<:Integer,T<:AbstractFloat,R<:PerturbationSolutionDet}

    if model.solvers == "Perturbation"
        error("The model hase been restricted to perturbation solvers only")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    f = zeros(nx+ny)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = domain[2,:] .+ rand(nx,npoints) .* (domain[1,:] - domain[2,:])

    @views for i = 1:npoints
        point = [states[:,i]; dynamics.g(states[:,i]); dynamics.h(states[:,i]); dynamics.gh(states[:,i])]
        f .= model.dynamic_function(point)
        euler_errors[:,i] .= f[model.eqns_approximated]
    end

    return euler_errors, states

end

function euler_errors(model::REModel,soln::R,domain::Union{Array{T,2},Array{T,1}},npoints::S,seed::S = 123456) where {S<:Integer,T<:AbstractFloat,R<:PerturbationSolutionStoch}

    if model.solvers == "Perturbation"
        error("The model hase been restricted to perturbation solvers only")
    end

    nx = length(soln.hbar)
    ny = length(soln.gbar)
    ns = size(soln.k,2)

    shocks = zeros(ns)
    f = zeros(nx+ny)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = domain[2,:] .+ rand(nx,npoints) .* (domain[1,:] - domain[2,:])

    @views for i = 1:npoints
        point = [states[:,i]; dynamics.g(states[:,i]); dynamics.h(states[:,i],shocks); dynamics.gh(states[:,i]); shocks]
        f .= model.dynamic_function(point)
        euler_errors[:,i] .= f[model.eqns_approximated]
    end

    return euler_errors, states

end

function euler_errors(model::REModel,soln::R,npoints::S,seed::S = 123456) where {S<:Integer,R<:ProjectionSolutionDet}

    if model.solvers == "Perturbation"
        error("The model hase been restricted to perturbation solvers only")
    end

    nx = size(soln.domain,2)
    nv = length(soln.variables)
    ny = nv - nx

    f = zeros(nv)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = soln.domain[2,:] .+ rand(nx,npoints) .* (soln.domain[1,:] - soln.domain[2,:])

    @views for i = 1:npoints
        point = [states[:,i]; dynamics.g(states[:,i]); dynamics.h(states[:,i]); dynamics.gh(states[:,i])]
        f .= model.dynamic_function(point)
        euler_errors[:,i] .= f[model.eqns_approximated]
    end

    return euler_errors, states

end

function euler_errors(model::REModel,soln::R,npoints::S,seed::S = 123456) where {S<:Integer,R<:ProjectionSolutionStoch}

    if model.solvers == "Perturbation"
        error("The model hase been restricted to perturbation solvers only")
    end

    nx = size(soln.domain,2)
    nv = length(soln.variables)
    ny = nv - nx

    shocks = zeros(ns)
    f = zeros(nv)

    dynamics = state_space_eqm(soln)

    euler_errors = zeros(length(model.eqns_approximated),npoints)

    Random.seed!(seed)
    states = soln.domain[2,:] .+ rand(nx,npoints).*(soln.domain[1,:] - soln.domain[2,:])

    @views for i = 1:npoints
        point = [states[:,i]; dynamics.g(states[:,i]); dynamics.h(states[:,i],shocks); dynamics.gh(states[:,i]); shocks]
        f .= model.dynamic_function(point)
        euler_errors[:,i] .= f[model.eqns_approximated]
    end

    return euler_errors, states

end

function den_haan_marcet(model::REModel,soln::R,steady_state::Array{T,1},seed::S = 123456) where {T<:AbstractFloat,S<:Integer,R<:Union{PerturbationSolutionStoch,ProjectionSolutionStoch}}

    nv = model.number_variables
    ns = model.number_shocks

    num_approx_eqns = length(model.eqns_approximated)

    shocks         = zeros(ns)
    f              = Array{T,1}(undef,nv)

    sim_length       = 3_500
    retained_length  = 3_000
    discarded_length = 500

    simulated_data = Array{T,2}(undef,nv,retained_length+1)

    num_tests = 1_000
    dhm_statistics = zeros(num_tests)

    Random.seed!(seed)
    for i = 1:num_tests
        random_seed = rand(1:10_000_000)
        errors  = zeros(num_approx_eqns)
        weights = zeros(num_approx_eqns,num_approx_eqns)
        simulated_data .= simulate(soln,steady_state[1:nx],sim_length+1,rndseed = random_seed)[:,(discarded_length+1):end]
        @views for j = 1:retained_length
            point = [simulated_data[:,j]; simulated_data[:,j+1]; shocks]
            f .= model.dynamic_function(point)
            errors  .+= f[model.eqns_approximated]
            weights .+= f[model.eqns_approximated]*f[model.eqns_approximated]'
        end
        dhm_statistics[i] = (errors'/weights)*errors
    end

    sort!(dhm_statistics,rev=true)

    dhm_stat = DenHaanMarcetStatistic(dhm_statistics[11],dhm_statistics[51],dhm_statistics[101],num_approx_eqns)

    return dhm_stat

end

function check_taylor_convergence(model::REModel,ss::Array{T,1},z::Array{T,1},degree::S) where {T<:Real,S<:Integer}

    if degree > 7
        error("The degree is too large (too many derivatives)")
    end
    if degree < 0
        error("The degree cannot be negative")
    end

    ne = model.number_equations
    nv = model.number_variables
    ns = model.number_shocks

    model_equations = model.each_eqn_function
    
    ss    = [ss;ss]
    z     = [z;z]
    point = [ss;zeros(ns)]

    taylor_terms = zeros(ne-ns,degree+1)

    taylor_terms[:,1] = model.dynamic_function(point)[ns+1:end]
    if degree == 0
        return taylor_terms
    else
        println("Skipping the shock processes.")
        @sync @qthreads for i = 1:ne-ns
            println("Working on equation ", ns+i)
            kprod = 1.0

            d(x) = ForwardDiff.gradient(model_equations[ns+i],x,ForwardDiff.GradientConfig(model_equations[ns+i],x,ForwardDiff.Chunk{2}()))[1:2*nv]
            dd(x) = ForwardDiff.jacobian(d,x,ForwardDiff.JacobianConfig(d,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddd(x) = ForwardDiff.jacobian(dd,x,ForwardDiff.JacobianConfig(dd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            dddd(x) = ForwardDiff.jacobian(ddd,x,ForwardDiff.JacobianConfig(ddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddddd(x) = ForwardDiff.jacobian(dddd,x,ForwardDiff.JacobianConfig(dddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            dddddd(x) = ForwardDiff.jacobian(ddddd,x,ForwardDiff.JacobianConfig(ddddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddddddd(x) = ForwardDiff.jacobian(dddddd,x,ForwardDiff.JacobianConfig(dddddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]

            for j = 1:degree
                kprod = kron(kprod,(z - ss))
                if j == 1
                    taylor_terms[i,j+1] = abs.(d(point))'kprod
                elseif j == 2
                    taylor_terms[i,j+1] = abs.(vec(dd(point)))'kprod/factorial(j)
                elseif j == 3
                    taylor_terms[i,j+1] = abs.(vec(ddd(point)))'kprod/factorial(j)
                elseif j == 4
                    taylor_terms[i,j+1] = abs.(vec(dddd(point)))'kprod/factorial(j)
                elseif j == 5
                    taylor_terms[i,j+1] = abs.(vec(ddddd(point)))'kprod/factorial(j)
                elseif j ==6
                    taylor_terms[i,j+1] = abs.(vec(dddddd(point)))'kprod/factorial(j)
                else
                    taylor_terms[i,j+1] = abs.(vec(ddddddd(point)))'kprod/factorial(j)
                end
            end
        end
        return taylor_terms
    end
end

function check_taylor_convergence(model::REModel,ss::Array{T,1},z::Array{T,1},degree::S,eqns::Array{S,1}) where {T<:Real,S<:Integer}

    if degree > 7
        error("The degree is too large (too many derivatives)")
    end
    if degree < 0
        error("The degree cannot be negative")
    end

    ne = model.number_equations
    nv = model.number_variables
    ns = model.number_shocks

    model_equations = model.each_eqn_function
    
    ss    = [ss;ss]
    z     = [z;z]
    point = [ss;zeros(ns)]

    taylor_terms = zeros(length(eqns),degree+1)

    taylor_terms[:,1] = model.dynamic_function(point)[eqns]
    if degree == 0
        return taylor_terms
    else
        println("Skipping the shock processes.")
        @sync @qthreads for i = 1:length(eqns)
            println("Working on equation ", eqns[i])
            kprod = 1.0

            d(x) = ForwardDiff.gradient(model_equations[eqns[i]],x,ForwardDiff.GradientConfig(model_equations[eqns[i]],x,ForwardDiff.Chunk{2}()))[1:2*nv]
            dd(x) = ForwardDiff.jacobian(d,x,ForwardDiff.JacobianConfig(d,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddd(x) = ForwardDiff.jacobian(dd,x,ForwardDiff.JacobianConfig(dd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            dddd(x) = ForwardDiff.jacobian(ddd,x,ForwardDiff.JacobianConfig(ddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddddd(x) = ForwardDiff.jacobian(dddd,x,ForwardDiff.JacobianConfig(dddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            dddddd(x) = ForwardDiff.jacobian(ddddd,x,ForwardDiff.JacobianConfig(ddddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]
            ddddddd(x) = ForwardDiff.jacobian(dddddd,x,ForwardDiff.JacobianConfig(dddddd,x,ForwardDiff.Chunk{2}()))[:,1:2*nv]

            for j = 1:degree
                kprod = kron(kprod,(z - ss))
                if j == 1
                    taylor_terms[i,j+1] = abs.(d(point))'kprod
                elseif j == 2
                    taylor_terms[i,j+1] = abs.(vec(dd(point)))'kprod/factorial(j)
                elseif j == 3
                    taylor_terms[i,j+1] = abs.(vec(ddd(point)))'kprod/factorial(j)
                elseif j == 4
                    taylor_terms[i,j+1] = abs.(vec(dddd(point)))'kprod/factorial(j)
                elseif j == 5
                    taylor_terms[i,j+1] = abs.(vec(ddddd(point)))'kprod/factorial(j)
                elseif j ==6
                    taylor_terms[i,j+1] = abs.(vec(dddddd(point)))'kprod/factorial(j)
                else
                    taylor_terms[i,j+1] = abs.(vec(ddddddd(point)))'kprod/factorial(j)
                end
            end
        end
        return taylor_terms
    end
end
