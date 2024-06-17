

function spacecraft_dynamics(P::ChaserTrajParams, x::Vector, u::Vector)
    return P.Ad * x[1:6] + P.Bd * u[1:3] ## chaser translation
end


function solve_cvx(P::ChaserTrajParams, X0c::Vector, target::Target, solver; silent_solver=false)

    N, dt = P.N, P.dt
    PCc, PCt = P.body.capture_port, target.body.capture_port

    # Variables
    X = cvx.Variable(6, N) 
    U = cvx.Variable(3, N-1) 
    sU = cvx.Variable(3, N-1) 
    ρ = cvx.Variable(N)

    # Compute expected final position of the chaser based on final predicted state of the target
    xt_f = predict_target_motion(target, P.dt * (N-1))
    xf = Vector(dcm_from_q(xt_f[1:4]) * (PCc + PCt))


    ################ Objective ################
    cost = 0.0
    unit = ones(3)

    for k=1:N

        if k != N
            # L1 norm control
            cost += P.ζ *  unit' * sU[:, k] 
        end

        # slack lifting variables 
        cost += ρ[k]

    end

    prob = cvx.minimize(cost)

    ################ Constraints ################

    ## Initial conditions
    prob.constraints += [X[:, 1] == X0c]

    ## Dynamics
    for k=1:N-1
        prob.constraints +=  [P.Ad * X[:, k] + P.Bd * U[:, k] == X[:, k+1]] 
    end

    ## State constraints
    for k=1:N
        prob.constraints += [norm(X[4:6, k]) <= P.v_max]
    end
    
    for k=1:N-1
        ## Norm control bound
        prob.constraints += [norm(U[:, k]) <= P.T_max]
        
        # Slack bounds on control (L1)
        prob.constraints += [sU[:, k] >= zeros(3)] # Positive slack
        prob.constraints += [ - sU[:, k] <= U[:, k]]
        prob.constraints += [ U[:, k] <= sU[:, k]]
    end


    
    ## Final position
    prob.constraints += [norm(X[1:3, N] - xf) <= P.goal_pos_tol] 
    ## Soft contact- near-zero relative velocity
    prob.constraints += [norm(X[4:6, N] - skew_symmetric(predict_target_motion(target, (N-1)*dt)[5:7]) * PCt) <= P.soft_vel_tol ]


    ## Docking cone constraints 
    A_dock = [1 0.0 0; 0 1 0; 0 0 0]
    c_dock = [0,0,1.0]

    if N - P.N_dock >= 0
        for k=N-P.N_dock+1:N
            wQb_t = dcm_from_q(predict_target_motion(target, (k-1)*dt)[1:4])
            Xw = X[1:3, k] - wQb_t * PCt
            prob.constraints += [ (norm((wQb_t*A_dock) * Xw )*tan((pi/2) - P.docking_angle) - (wQb_t*c_dock)' * Xw ) <= 0]
        end
    else
        for k=1:N
            wQb_t = dcm_from_q(predict_target_motion(target, (k-1)*dt)[1:4])
            Xw = X[1:3, k] - wQb_t * PCt
            prob.constraints += [ (norm((wQb_t*A_dock) * Xw )*tan((pi/2) - P.docking_angle) - (wQb_t*c_dock)' * Xw ) <= 0]
        end
    end


    ## Lifting variables - relaxation
    for k=1:N
        prob.constraints += [ 0.0 <= ρ[k] ]
        prob.constraints += [ρ[k] <= P.p_max]
        prob.constraints += [norm(X[1:3,k]) <= ρ[k]]
    end
    
    
    ## Field-of-View constraints (convex relaxation)
    β = P.ω_max * dt
    H = Matrix([I -0.5*I(3); -0.5*I(3) I])
    H_sqrt = sqrt(H)
    S = [1 -0.5*cos(β); -0.5*cos(β) 1]
    e, _ = eigen(sqrt(S))
    σmin = minimum(e)
    for k=1:N-1
        prob.constraints += [norm(H_sqrt * [X[1:3, k+1];X[1:3, k]]) <= σmin * (1/sqrt(2)) * (ρ[k+1]+ρ[k])] # SOCP
    end



    ################ Output ################
    cvx.solve!(prob, solver; silent = silent_solver)

    status = nothing

    if (prob.status != cvx.MOI.OPTIMAL) && (prob.status != cvx.MOI.ALMOST_OPTIMAL) && (prob.status != cvx.MOI.SLOW_PROGRESS) 
        #error("Convex.jl problem failed to solve.")
        print("Convex.jl problem failed to solve.")
        status = false
        return status, prob, nothing, nothing, nothing, nothing, nothing
    end
    # SLOW_PROGRESS should be checked again


    status = true
    X = Matrix(transpose(X.value))
    U = Matrix(transpose(U.value))
    sU = Matrix(transpose(sU.value))
    ρ = Matrix(ρ.value)
    J = prob.optval

    return status, prob, X, U, sU, ρ, J
end


function solve_cvx_correction(P::ChaserTrajParams, Xref::Matrix, Uref::Matrix, target::Target, solver; silent_solver=false)

    N, dt = P.N, P.dt
    PCc, PCt = P.body.capture_port, target.body.capture_port

    # Variables
    dX = cvx.Variable(6, N) 
    dU = cvx.Variable(3, N-1) 
    sU = cvx.Variable(3, N-1) 
    ρ = cvx.Variable(N)

    Sl = cvx.Variable(N)  # Penalized slack

    # Compute expected final position of the chaser based on final predicted state of the target
    xt_f = predict_target_motion(target, P.dt * (N-1))
    xf = Vector(dcm_from_q(xt_f[1:4]) * (PCc + PCt))


    ################ Objective ################
    cost = 0.0
    unit = ones(3)

    for k=1:N

        if k != N
            # L1 norm control
            cost += P.ζ *  unit' * sU[:, k] 
        end

        # slack lifting variables 
        cost += ρ[k]

        cost += P.ψ * Sl[k]

    end

    prob = cvx.minimize(cost)


    ################ Constraints ################

    ## Initial conditions
    prob.constraints += [dX[:, 1] == zeros(6)]

    ## Dynamics
    for k=1:N-1
        prob.constraints +=  [P.Ad * dX[:, k] + P.Bd * dU[:, k] == dX[:, k+1]] 
    end

    ## State constraints
    for k=1:N
        prob.constraints += [norm(Xref[k, 4:6] + dX[4:6, k]) <= P.v_max]
    end
    
    for k=1:N-1
        ## Norm control bound
        prob.constraints += [norm(Uref[k, :] + dU[:, k]) <= P.T_max]
        
        # Slack bounds on control (L1)
        prob.constraints += [sU[:, k] >= zeros(3)] # Positive slack
        prob.constraints += [ - sU[:, k] <= Uref[k, :] + dU[:, k]]
        prob.constraints += [ Uref[k, :] + dU[:, k] <= sU[:, k]]
    end


    
    ## Final position
    prob.constraints += [norm(Xref[N, 1:3] + dX[1:3, N] - xf) <= P.goal_pos_tol] 
    ## Soft contact- near-zero relative velocity
    prob.constraints += [norm(Xref[N, 4:6] + dX[4:6, N] - skew_symmetric(predict_target_motion(target, (N-1)*dt)[5:7]) * PCt) <= P.soft_vel_tol]


    ## Docking cone constraints 
    A_dock = [1 0.0 0; 0 1 0; 0 0 0]
    c_dock = [0,0,1.0]

    if N - P.N_dock >= 0
        for k=N-P.N_dock+1:N
            wQb_t = dcm_from_q(predict_target_motion(target, (k-1)*dt)[1:4])
            Xw = Xref[k, 1:3] + dX[1:3, k] - wQb_t * PCt
            prob.constraints += [ (norm((wQb_t*A_dock) * Xw )*tan((pi/2) - P.docking_angle) - (wQb_t*c_dock)' * Xw ) <= 0]
        end
    else
        for k=1:N
            wQb_t = dcm_from_q(predict_target_motion(target, (k-1)*dt)[1:4])
            Xw = Xref[k, 1:3] + dX[1:3, k] - wQb_t * PCt
            prob.constraints += [ (norm((wQb_t*A_dock) * Xw )*tan((pi/2) - P.docking_angle) - (wQb_t*c_dock)' * Xw ) <= 0]
        end
    end


    ## Lifting variables - relaxation
    for k=1:N
        prob.constraints += [ 0.0 <= ρ[k] ]
        prob.constraints += [ρ[k] <= P.p_max]
        prob.constraints += [norm(Xref[k, 1:3] + dX[1:3,k]) <= ρ[k]]


        # Slack for collision avoidance
        prob.constraints += [Sl[k] >= 0]
    end
    
    
    ## Field-of-View constraints (convex relaxation)
    β = P.ω_max * dt
    H = Matrix([I -0.5*I(3); -0.5*I(3) I])
    H_sqrt = sqrt(H)
    S = [1 -0.5*cos(β); -0.5*cos(β) 1]
    e, _ = eigen(sqrt(S))
    σmin = minimum(e)
    for k=1:N-1
        prob.constraints += [norm(H_sqrt * [Xref[k, 1:3] + dX[1:3, k];Xref[k+1, 1:3] + dX[1:3, k+1]]) <= σmin * (1/sqrt(2)) * (ρ[k+1]+ρ[k])] # SOCP
    end



    ## Collision avoidance constraints (dcol) - linearized 
    q_c = get_attitude_pointing(Xref)
    for k = 2:N
        target.body.hull.q = SVector{4, eltype(Xref)}(predict_target_motion(target, (k-1)*P.dt)[1:4])
        P.body.hull.r = SVector{3, eltype(Xref)}(Xref[k, 1:3])
        P.body.hull.q = SVector{4, eltype(Xref)}(q_c[k, :])
        αbar, x, J = dc.proximity_jacobian(P.body.hull,target.body.hull)

        jac_qc_pc = FD.jacobian(p -> get_next_fov_attitude(q_c[k-1, :], Xref[k-1, 1:3], p), Xref[k, 1:3]) # 4 x 3
        jac_col = J[4, 1:3]' + J[4, 4:7]' * jac_qc_pc
        prob.constraints += [ P.α_colmin <=  αbar + jac_col * dX[1:3, k] + Sl[k]] 
    end 


    ################ Output ################

    cvx.solve!(prob, solver; silent = silent_solver)
    status = nothing

    if (prob.status != cvx.MOI.OPTIMAL) && (prob.status != cvx.MOI.ALMOST_OPTIMAL) && (prob.status != cvx.MOI.SLOW_PROGRESS) 
        #error("Convex.jl problem failed to solve.")
        print("Convex.jl problem failed to solve.")
        status = false
        return status, prob, nothing, nothing, nothing, nothing, nothing, nothing
    end
    # SLOW_PROGRESS should be checked again

    status = true
    dX = Matrix(transpose(dX.value))
    dU = Matrix(transpose(dU.value))
    sU = Matrix(transpose(sU.value))
    ρ = Matrix(ρ.value)
    Sl= Matrix(Sl.value)
    J = prob.optval

    return status, prob, dX, dU, sU, ρ, Sl, J
end



# Full trajectory Optmisation

function solve_trajectory_optimization(P::ChaserTrajParams, X0_chaser::Vector, target::Target, solver::Function; max_iter_scp::Real=14, silent_solver::Bool=false, get_violation_dcol::Bool=true)

    Jcost = [] # costs
    dcols = [] # violations 

    ## Initial convex iteration
    println("TrajOpt: CVX iteration")
    status, prob, X, U, sU, ρ, J = solve_cvx(P, X0_chaser, target, solver(); silent_solver=silent_solver)

    if !status
        print("TrajOpt: Initial convex optimization failed")
        return status, nothing, nothing, nothing, nothing
    end

    dcol = collision_avoidance_dcol(P, X, target)
    dcol_violation = sum(min.(dcol.-1, 0.0)) # Use absolute safety criteria

    push!(dcols, dcol_violation)
    push!(Jcost, J)

    i = 1

    ## Inner loop - SCP corrections 
    while ((dcol_violation < 0.0) && i < max_iter_scp)

        println("TrajOpt: SCP iteration ", i)
        status, prob, dX, dU, sU, ρ, Sl, J = solve_cvx_correction(P, X, U, target, solver(); silent_solver=silent_solver)

        if !status
            print("TrajOpt: SCP convex optimization failed at iteration ", i)
            return status, nothing, nothing, nothing, nothing
        end

    
        X = X + dX
        U = U + dU

        dcol = collision_avoidance_dcol(P, X, target)
        dcol_violation = sum(min.(dcol.-1, 0.0)) # Use absolute safety criteria

        push!(Jcost, J)
        push!(dcols, dcol_violation)

        i = i + 1
    end



    return status, X, U, Jcost, dcols

end




"""
    collision_avoidance_dcol(P::ChaserTrajParams, X::Matrix, target::Target)

Compute the collision avoidance metric for each time step in the trajectory.

# Arguments
- `P::ChaserTrajParams`: Parameters of the chaser trajectory.
- `X::Matrix`: Matrix representing the chaser trajectory.
- `target::Target`: Target object.

# Returns
- `αviol::Vector`: Vector of collision avoidance metrics for each time step.

"""
function collision_avoidance_dcol(P::ChaserTrajParams, X::Matrix, target::Target)

    # DCOL
    q_c = get_attitude_pointing(X)
    αviol = zeros(eltype(X), P.N)

    for k = 1:P.N
        target.body.hull.q = SVector{4, eltype(X)}(predict_target_motion(target, (k-1)*P.dt)[1:4])
        P.body.hull.r = SVector{3, eltype(X)}(X[k, 1:3])
        P.body.hull.q = SVector{4, eltype(X)}(q_c[k, :])
        αbar, x, J = dc.proximity_jacobian(P.body.hull, target.body.hull)
        αviol[k] = αbar
    end

    return αviol
end




function get_attitude_pointing(X::Matrix; N::Real=size(X,1))

    q = zeros(eltype(X), N, 4)

    Pt = [0.0, 0, 1] # z-axis of reference global frame
    Pc = - X[1, 1:3] / norm(X[1, 1:3]) ## Pointing axis (z of cone) 

    # Get initial attitude
    c = cross(Pt, Pc)
    n = c / norm(c)
    θ = atan(norm(c), dot(Pc,Pt))
    q[1, :] = [cos(θ/2) ; n*sin(θ/2)]
    q[1, :] = q[1, :]/  norm(q[1, :])

    for k=2:N
        q[k, :] = get_next_fov_attitude(q[k-1, :], X[k-1, 1:3], X[k, 1:3])
    end
     
    return q
end



function get_next_fov_attitude(q0::Vector, p0::Vector, p_next::Vector)
    # Produce the next attitude based on the current and next position

    Pc = - p0 / norm(p0) ## Current pointing axis (z of cone) 
    Pc1 = -p_next / norm(p_next)

    c = cross(Pc, Pc1)

    n = c / norm(c)
    θ = atan(norm(c), dot(Pc, Pc1))
    dq = [cos(θ/2) ; n*sin(θ/2)]
    dq = dq / norm(dq)


    q_next = R(q0) * dq

    return q_next

end

function zero_order_hold_interpolation(dt_ratio::Real, U::Matrix)
    u_plan = vcat([repeat(row', dt_ratio, 1) for row in eachrow(U)]...)
    return u_plan
end

function simulate_trajectory(P::ChaserTrajParams, target::Target, X0c::Vector, X0t::Vector, Uc::Matrix, dt_sim::Real)

    dt_ratio = Int(P.dt / dt_sim)
    N_sim = (P.N-1) * dt_ratio

    Xc = zeros(eltype(X0c), N_sim, 6)
    qc = zeros(eltype(X0c), N_sim, 4)
    Xt = zeros(eltype(X0t), N_sim, 7)

    U = zero_order_hold_interpolation(dt_ratio, Uc)

    Xc[1, :] = X0c
    Xt[1, :] = X0t

    u0 = zeros(eltype(X0t), 3)

    for k=1:N_sim-1
        Xc[k+1, :] = rk4_chaser(P.body.mass, target.sma, Xc[k, :], U[k, :], dt_sim)
        Xt[k+1, :] = rk4mk_target(target.body.J, Xt[k, :], u0, dt_sim)
    end

    qc = get_attitude_pointing(Xc)

    return Xc, qc, U, Xt
end
