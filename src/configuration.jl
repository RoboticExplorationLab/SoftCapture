

struct BodyParams
    mass::Real
    J::Matrix
    capture_port::Vector # Position of the capture port w.r.t. the center of mass in the body frame
    hull::dc.AbstractPrimitive # Convex hull for collision avoidance (DifferentiableCollisions.jl)
    function BodyParams(mass::Real, J::Matrix, capture_port::Vector, hull::dc.AbstractPrimitive)
        new(mass, J, capture_port, hull)
    end
end


mutable struct StateIdx
    x
    u
    t
    nx
    nu
    N
    nz

    function StateIdx(nx, nu, N)
        # This function creates some useful indexing tools for Z 
        nz = (N) * nu + N * nx # length of Z 
        x = [(i - 1) * (nx + nu) .+ (1 : nx) for i = 1:N]
        u = [(i - 1) * (nx + nu) .+ ((nx + 1):(nx + (nu))) for i = 1:(N-1)] # zero order hold 
        t = []
        return new(x,u,t,nx,nu,N,nz)
    end
end


mutable struct ChaserTrajParams

    idx::StateIdx # Indexing 
    dt::Real # Time step
    N::Int64 # Number of knot points
    Ad::Matrix # state transition matrix
    Bd::Matrix # input matrix
    body::BodyParams # Chaser body parameters


    ζ::Real # Weight on the L1 control 
    ψ::Real # Weight on slack variable

    goal_pos_tol::Real # Tolerance on the goal position
    soft_vel_tol::Real # Tolerance on the velocity (soft capture)

    α_colmin::Real # Minimum scaling factor between the chaser and the target (must be greater than 1)
    fov_cone_angle::Real # Field of view cone angle [rad]
    docking_angle::Real # Docking cone angle [rad]
    N_dock::Int64 # Last N_dock knot points for the docking constraints 
    p_max::Real # Maximum position [m] (artifical)
    T_max::Real # Maximum thrust [N]
    v_max::Real # Maximum velocity [m/s]
    ω_max::Real # Maximum angular velocity [rad/s]
    

    function ChaserTrajParams(body::BodyParams, dt::Real, N::Int64, target_sma::Real)

        idx = StateIdx(6, 3, N)
        Ad, Bd = discretized_clohessy_wiltshire(dt, body.mass, target_sma)

        ζ = 5 # default
        ψ = 750 # default
        
        # Constraint values and tolerances 
        goal_pos_tol = 0.35
        soft_vel_tol = 0.03
        α_colmin = 1.3
        fov_cone_angle = deg2rad(60)
        docking_angle = deg2rad(30)
        N_dock = 5
        p_max = 100.0
        T_max = 100.0
        v_max = 1.5 # 1 m/s = 3.6 km/h, 1.5 m/s = 5.4 km/h
        ω_max = 0.2

        new(idx, dt, N, Ad, Bd, body, ζ, ψ, goal_pos_tol, soft_vel_tol, α_colmin, fov_cone_angle, docking_angle, N_dock, p_max, T_max, v_max, ω_max)

    end


end



mutable struct Target

    body::BodyParams # Target body parameters
    dt::Real # Time step
    N::Int64 # Number of knot points
    X0::Vector # Initial target state
    max_time_propagation::Real # Maximum time for propagation
    nominal_motion::Matrix
    sma::Real # Semi-major axis of the target orbit
    

    function Target(body::BodyParams, dt::Real,  X0::Vector, sma::Real; max_time_propagation::Real=300.0)

        N = Int64(round(max_time_propagation/dt)+1)

        X_target = zeros(eltype(X0), N, 7)
        X_target[1,:] = X0
        u0 = zeros(eltype(X_target), 3)

        for k=1:1:N-1
            X_target[k+1,:] = rk4mk_target(body.J, X_target[k,:], u0, dt)
        end

        new(body, dt, N, X0, max_time_propagation, X_target, sma)
        
    end

end



function predict_target_motion(target::Target, time_request)

    # Safeguards removed since the time is constrained in the optimization
    Z = target.nominal_motion
    ## Locate timestep
    n_find = (time_request / target.dt) + 1
    if n_find < 1 # Dirty fix
        n_find = 1
    end
    n_lb = Int(floor(n_find))
    n_ub = Int(ceil(n_find))

    if n_lb == n_ub
        return Z[n_lb, :]
    else
        ## Fit spline between the two points
        xt = n_lb:n_ub-n_lb:n_ub
        ## Angular velocity
        yt = [Z[n_lb,5:7] Z[n_ub,5:7]]
        n = length(Z[n_lb,5:7])    
        itps = [CubicSplineInterpolation(xt, yt[i,:]) for i in 1:n]
        ω = [itp[n_find] for itp in itps]
        # quaternions
        q = slerp(Z[n_lb,1:4], Z[n_ub,1:4], n_find-n_lb)
        return [q;ω]
    end
end



function repropagate_target!(target::Target, X0::Vector)

    target.X0 = X0
    xt = zeros(eltype(X0), target.N, 7)
    xt[1, :] = X0
    u0 = zeros(eltype(X0), 3)

    for i=1:target.N-1
        xt[i+1, :] = rk4mk_target(target.body.J, xt[i, :], u0, target.dt)
    end
    target.nominal_motion = xt
end


function get_predicted_target_trajectory(target::Target, N::Int64, dt::Real)
    Xt = zeros(eltype(target.X0), N, 7)
    for i=1:N
        Xt[i, :] = predict_target_motion(target, (i-1)*dt)
    end
    return Xt
end


"""
generate_relative_ic_safety_ellipse(sma::Real; A0_min::Real=25, A0_max::Real=50, B0_min::Real=50, B0_max::Real=100)

Generate initial conditions for a relative motion safety ellipse. 
A0 - along radial direction Lx
2A0 - along-track direction Ly
B0 - along cross-tack direction Lz (out-of-plane motion)

Parameters:
    - `sma::Real`: Semi-major axis of the ellipse.
    - `A0_min::Real`: Minimum value for the semi-major axis along the radial direction.
    - `A0_max::Real`: Maximum value for the semi-major axis along the radial direction.
    - `B0_min::Real`: Minimum value for the semi-major axis along the cross-tack direction.
    - `B0_max::Real`: Maximum value for the semi-major axis along the cross-tack direction.

Returns:
    - The vector of initial conditions [x0, y0, z0, ẋ0, ẏ0, ż0] for the relative motion.

"""
function generate_relative_ic_safety_ellipse(sma::Real; A0_min::Real=25, A0_max::Real=50, B0_min::Real=50, B0_max::Real=100)

    # Generate random values for A0 and B0 within the specified range
    A0 = (rand(Bool) ? 1 : -1) * (A0_min + rand() * (A0_max-A0_min))
    B0 = (rand(Bool) ? 1 : -1) * (B0_min + rand() * (B0_max-B0_min))

    # Calculate the mean motion
    n = sqrt(earth.μ / sma^3)
    
    # Generate values for x0, ẋ0, y0, ẏ0, z0, and ż0 within the specified range
    x0 = (rand(Bool) ? 1 : -1) * (A0/2 + (rand() * A0/2))
    ẋ0 = (rand(Bool) ? 1 : -1) * (n * sqrt(A0^2 - x0^2))
    y0 = 2 * ẋ0 / n 
    ẏ0 = - 2 * n * x0
    z0 = (rand(Bool) ? 1 : -1) * (B0/2 + (rand() * B0/2))
    ż0 = (rand(Bool) ? 1 : -1) * (n * sqrt(B0^2 - z0^2))

    return [x0, y0, z0, ẋ0, ẏ0, ż0]

end

"""
    sample_initial_tumble_state(max_ω_norm::Real)

This function generates a random initial tumble state within the given maximum bound.

# Arguments
- `max_ω_norm::Real`: The maximum value for the norm of the angular velocity vector.

# Returns
- `state::Vector{Float64}`: The generated initial tumble state (7x1): [attitude quaternion, angular velocity vector]

"""

function sample_initial_tumble_state(max_ω_norm::Real)
    # Attitude and angular velocity 

    r_axisangle = (rand(3).*2*π).-π
    θ_n = norm(r_axisangle)
    axis_n = r_axisangle / θ_n
    q = [cos(0.5 * θ_n); axis_n * sin(0.5 * θ_n)]
    
    # Sample random direction for the angular velocity vector
    rω = rand(3).*2 .- 1
    rω = rω / norm(rω)

    norm_sample = rand()*max_ω_norm 
    rω = rω * norm_sample
    ω = rω

    return [q;ω]
end

