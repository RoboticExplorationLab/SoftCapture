# Dynamics models


function rk4(params::NamedTuple, dynamics::Function, x::Vector, u::Vector, dt::Real)
    # Implementation of the fourth-order Runge-Kutta method for numerical integration of ordinary differential equations.
    k1 = dt * dynamics(params, x, u)
    k2 = dt * dynamics(params, x + k1 * 0.5, u)
    k3 = dt * dynamics(params, x + k2 * 0.5, u)
    k4 = dt * dynamics(params, x + k3, u)
    x = x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return x
end

function euler_rotational_dynamics(params::NamedTuple, ω::Vector, τ::Vector)
    # Computes the rotational dynamics using Euler's method.
    return params.J \ (τ - cross(ω, params.J * ω))
end

function clohessy_wiltshire(params::NamedTuple, x::Vector, u::Vector)
    """
    Computes the dynamics of a satellite in a circular orbit using the Clohessy-Wiltshire equations.

    Args:
        params (NamedTuple): Parameters of the satellite and the central body.
        x (Vector): State vector of the satellite.
        u (Vector): Control input vector.

    Returns:
        Vector: Time derivative of the state vector.
    """
    r = x[1:3] # Position
    ṙ = x[4:6] # Linear Velocity

    n = sqrt(earth.μ / params.sma^3)
    A = [0 0 0 1 0 0; 
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        3*n^2 0 0 0 2*n 0;
        0 0 0 -2*n 0 0;
        0 0 -n^2 0 0 0]
    
    B = [zeros(eltype(u), 3, 3); Matrix{eltype(x)}(I, 3, 3)]
    cw = A*[r ; ṙ] + B * (1/ earth.m) * u
    r̈ = cw[4:6]
    ẋ = [ṙ ; r̈ ]

    return ẋ
end


function clohessy_wiltshire_matrices(params::NamedTuple)
    """
    Computes the matrices A and B for the Clohessy-Wiltshire equations.

    Args:
        params (NamedTuple): Parameters of the satellite and the central body.

    Returns:
        Tuple: Matrices A and B.
    """
    n = sqrt(earth.μ / params.sma^3)
    A = [0 0 0 1 0 0; 
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        3*n^2 0 0 0 2*n 0;
        0 0 0 -2*n 0 0;
        0 0 -n^2 0 0 0]
    
    B = [zeros(3, 3); Matrix(I, 3, 3)]./params.m
    return A, B
end


function discretized_clohessy_wiltshire(params::NamedTuple, dt::Real)
    """
    Computes the discretized A and B matrices for the Clohessy-Wiltshire equations.

    Args:
        params (NamedTuple): Parameters of the satellite and the central body.
        dt (Real): Time step.

    Returns:
        Tuple: Discretized matrices Ad and Bd.
    """
    A, B = clohessy_wiltshire_matrices(params)

    nx = size(A)[1]
    nu = size(B)[2]

    M = zeros(nx+nu,nx+nu)
    M[1:nx,1:nx] = A
    M[1:nx,1+nx:nx+nu] = B

    M = exp(M*dt)

    Ad = M[1:nx,1:nx]
    Bd = M[1:nx,1+nx:nx+nu]

    return Ad, Bd
end

function discretized_clohessy_wiltshire(dt::Real, m::Real, sma::Real)
    """
    Computes the discretized A and B matrices for the Clohessy-Wiltshire equations.

    Args:
        dt (Real): Time step.
        m (Real): Mass of the spacecraft.
        sma (Real): Semi-major axis of the target's orbit.

    Returns:
        Tuple: Discretized matrices Ad and Bd.
    """
    n = sqrt(earth.μ / sma^3)
    A = [0 0 0 1 0 0; 
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        3*n^2 0 0 0 2*n 0;
        0 0 0 -2*n 0 0;
        0 0 -n^2 0 0 0]
    
    B = [zeros(3, 3); Matrix(I, 3, 3)]./m

    nx = size(A)[1]
    nu = size(B)[2]

    M = zeros(nx+nu,nx+nu)
    M[1:nx,1:nx] = A
    M[1:nx,1+nx:nx+nu] = B

    M = exp(M*dt)

    Ad = M[1:nx,1:nx]
    Bd = M[1:nx,1+nx:nx+nu]

    return Ad, Bd
end

function quat_kinematics(q::Vector, ω::Vector; H::Matrix=[zeros(eltype(q), 1, 3); Matrix{eltype(q)}(I, 3, 3)])
    """
    q: (w,x,y,z)
    """
    q̇ = 0.5 * L(q) * H * ω
    return q̇
end


function attitude_quat_rbd(params::NamedTuple, x::Vector, τ::Vector; invJ=false)
    """
    Computes the rotational dynamics of an attitude quaternion using the rigid body dynamics equations.

    Args:
        params (NamedTuple): Parameters of the satellite.
        x (Vector): State vector containing the attitude quaternion and angular velocity.
        τ (Vector): Control input vector.
        invJ (Bool, optional): Whether to use the inverse of the inertia matrix. Defaults to false.

    Returns:
        Vector: Time derivative of the state vector.
    """
    q = x[1:4] # Attitude (quaternions)
    ω = x[5:7] # Angular Velocity
    q̇ = quat_kinematics(q, ω)
    if invJ
        ωdot = params.invJ * (τ - cross(ω, params.J * ω))
    else
        ωdot = params.J \ (τ - cross(ω, params.J * ω))
    end
    
    return [q̇ ; ωdot]
end


function rk4mk_attitude(params::NamedTuple, dynamics::Function, x::Vector, u::Vector, dt::Real)
    """
    Implements the fourth-order Runge-Kutta method for numerical integration of the rotational dynamics with quaternion kinematics.

    Args:
        params (NamedTuple): Parameters of the satellite.
        dynamics (Function): Function that computes the rotational dynamics.
        x (Vector): State vector containing the attitude quaternion and angular velocity.
        u (Vector): Control input vector.
        dt (Real): Time step.

    Returns:
        Vector: Updated state vector.
    """
    H = [zeros(eltype(x), 1, 3); Matrix{eltype(x)}(I, 3, 3)]
    k1 = dt * dynamics(params, x, u)
    k2 = dt * dynamics(params, x + k1 * 0.5, u)
    k3 = dt * dynamics(params, x + k2 * 0.5, u)
    k4 = dt * dynamics(params, x + k3, u)

    f = (x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)) 

    q = x[1:4]
    ω = x[5:7]

    ## Using ExponentialAction - ForwardDiff friendly https://github.com/sethaxen/ExponentialAction.jl 
    f[1:4] = expv(1, R( H * (0.5 * dt * ω + ((dt)/6)*(k1[5:7] + k2[5:7] + k3[5:7]))), q) 

    return f 
end

function target_dynamics(J::Matrix, x::Vector, τ::Vector)
    """
    Computes the rotational dynamics of an attitude quaternion using the rigid body dynamics equations.

    Args:
        J (Matrix): Inertia matrix.
        x (Vector): State vector containing the attitude quaternion and angular velocity.
        τ (Vector): Control input vector.

    Returns:
        Vector: Time derivative of the state vector.
    """
    q = x[1:4] # Attitude (quaternions)
    ω = x[5:7] # Angular Velocity
    q̇ = quat_kinematics(q, ω)
    ωdot = J \ (τ - cross(ω, J * ω))
    return [q̇ ; ωdot]
end

function rk4mk_target(J::Matrix, x::Vector, u::Vector, dt::Real)
    """
    Implements the fourth-order Runge-Kutta method for numerical integration of the rotational dynamics with quaternion kinematics.
    This function is used for the target dynamics.

    Args:
        J (Matrix): Inertia matrix.
        x (Vector): State vector containing the attitude quaternion and angular velocity.
        u (Vector): Control input vector.
        dt (Real): Time step.

    Returns:
        Vector: Updated state vector.
    """
    H = [zeros(eltype(x), 1, 3); Matrix{eltype(x)}(I, 3, 3)]
    k1 = dt * target_dynamics(J, x, u)
    k2 = dt * target_dynamics(J, x + k1 * 0.5, u)
    k3 = dt * target_dynamics(J, x + k2 * 0.5, u)
    k4 = dt * target_dynamics(J, x + k3, u)

    f = (x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)) 

    q = x[1:4]
    ω = x[5:7]

    ## Using ExponentialAction - ForwardDiff friendly https://github.com/sethaxen/ExponentialAction.jl 
    f[1:4] = expv(1, R( H * (0.5 * dt * ω + ((dt)/6)*(k1[5:7] + k2[5:7] + k3[5:7]))), q) 

    return f 
end


function nonlinear_relative_keplerian_dynamics(mass::Real, target_sma::Real, X::Vector, F::Vector)
    """
    Computes the nonlinear relative Keplerian dynamics.

    Args:
        mass (Real): Mass.
        target_sma (Real): Semi-major axis of the target's orbit.
        params (NamedTuple): Parameters of the satellite and the central body.
        X (Vector): State vector containing the relative position and velocity.
        F (Vector): Control input vector.

    Returns:
        Vector: Time derivative of the state vector.
    """
    x = X[1]
    y = X[2]
    z = X[3]
    ẋ = X[4]
    ẏ = X[5]
    ż = X[6]

    rt = target_sma
    rtdot = sqrt(earth.μ/rt)

    rc = sqrt((rt + x)^2 + y^2 + z^2)
    θdot = sqrt(earth.μ/rt^3)

    ẍ = - (earth.μ / rc^3) * ( rt + x ) +  θdot^2 * x + 2 * θdot * (ẏ - y * (rtdot/rt)) + earth.μ / rt^2 + F[1] / mass
    ÿ = - (earth.μ / rc^3) * y + θdot^2 * y - 2 * θdot * (ẋ - x * (rtdot/rt)) + F[2] / mass
    z̈ = - (earth.μ / rc^3) * z + F[3] / mass

    return [ẋ, ẏ, ż, ẍ, ÿ, z̈]
end

function rk4_chaser(mass::Real, target_sma::Real, x::Vector, u::Vector, dt::Real)
    """
    Implements the fourth-order Runge-Kutta method for numerical integration of the rotational dynamics with quaternion kinematics.
    This function is used for the chaser dynamics.

    Args:
        mass (Real): Mass.
        target_sma (Real): Semi-major axis of the target's orbit.
        x (Vector): State vector containing the attitude quaternion and angular velocity.
        u (Vector): Control input vector.
        dt (Real): Time step.

    Returns:
        Vector: Updated state vector.
    """
    H = [zeros(eltype(x), 1, 3); Matrix{eltype(x)}(I, 3, 3)]
    k1 = dt * nonlinear_relative_keplerian_dynamics(mass, target_sma, x, u)
    k2 = dt * nonlinear_relative_keplerian_dynamics(mass, target_sma, x + k1 * 0.5, u)
    k3 = dt * nonlinear_relative_keplerian_dynamics(mass, target_sma, x + k2 * 0.5, u)
    k4 = dt * nonlinear_relative_keplerian_dynamics(mass, target_sma, x + k3, u)

    f = (x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)) 
    return f 
end


function scale_state(x0, r_scale, t_scale)
    """
    Scales the Cartesian state vector.

    Args:
        x0 (Vector): Cartesian state vector.
        r_scale (Real): Scaling factor for position.
        t_scale (Real): Scaling factor for time.

    Returns:
        Vector: Scaled state vector.
    """
    r = x0[1:3] / r_scale
    v = x0[4:6] / (r_scale / t_scale)
    return [r; v]
end

function unscale_state(x0, r_scale, t_scale)
    """
    Unscales the Cartesian state vector.

    Args:
        x0 (Vector): Scaled state vector.
        r_scale (Real): Scaling factor for position.
        t_scale (Real): Scaling factor for time.

    Returns:
        Vector: Unscaled state vector.
    """
    r = x0[1:3] * r_scale
    v = x0[4:6] * (r_scale / t_scale)
    return [r; v]
end
