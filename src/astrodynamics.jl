

const earth = (
    μ = 3.986004418e14,
    G = 6.67430e-11,
    M = 5.9722e24,
    R = 6378.0e3,
    R_GEO = 42164000.0,
    R_LEO_ub = 8413.0e3
)

# Some global constants here 

const LEO_lb = 6778.0e3 
const LEO_ub = 8413.0e3

const MEO = 20200e3
const GSO = 42164e3


mutable struct Orbit 

    region::String
    sma::Real
    ecc::Real
    i::Real # Inclination
    raan::Real # Right ascension of the ascending node 
    ω::Real # Argument of perigee
    anomaly::Real # True anomaly 
    cartesian_state::Vector

    function Orbit(region="None", sma=0.0,ecc=0.0,i=0,raan=0,ω=0.0,anomaly=0.0, cartesian_state=zeros(6))
        return new(region,sma,ecc,i,raan,ω,anomaly,cartesian_state)
    end

end

function sample_new_orbit(region="any")

    # Sampling over orbital elements with specific constraints 

    orbit = Orbit()

    c = rand()
    if region == "any"
        if c < 0.45
            orbit.region  = "LEO"
            # LEO - uniform distriution
            orbit.sma = LEO_lb + rand() * (LEO_ub - LEO_lb)
            orbit.ecc = rand() * 0.25
        elseif c >= 0.45 && c <= 0.55
            orbit.region = "MEO"
            # MEO - Gaussian tranform 
            orbit.sma = rand(Dis.Normal(MEO, 1000e3))
            orbit.ecc = rand() * 0.5
        else
            orbit.region = "GSO"
            # GSO - Gaussian here? 
            orbit.sma = rand(Dis.Normal(GSO, 100e3))
            orbit.ecc = rand() * 0.25
        end
    elseif region == "LEO"
        orbit.region  = "LEO"
        # LEO - uniform distriution
        orbit.sma = LEO_lb + rand() * (LEO_ub - LEO_lb)
        orbit.ecc = rand() * 0.25
    elseif region == "MEO"
        orbit.region = "MEO"
        # MEO - Gaussian tranform 
        orbit.sma = rand(Dis.Normal(MEO, 1000e3))
        orbit.ecc = rand() * 0.5
    elseif region == "GSO"
        orbit.region = "GSO"
        # GSO - Gaussian here? 
        orbit.sma = rand(Dis.Normal(GSO, 100e3))
        orbit.ecc = rand() * 0.25
    else
        error("Specify a correct region: LEO, MEO, GSO or any.")

    end

    

    orbit.i = rand() * 2 * π
    orbit.raan = rand() * 2 * π

    if orbit.ecc != 0.0
        orbit.ω = rand() * 2 * π
        orbit.anomaly = rand() * 2 * π
    end

    

    orbit.cartesian_state = get_CART_from_OSC([
        orbit.sma;
        orbit.ecc;
        orbit.i;
        orbit.raan;
        orbit.ω;
        orbit.anomaly
    ])


    return orbit 

end


"""
    semi_major_axis(rp::Float64, ra::Float64)::Float64

Calculate the semi-major axis of an ellipse given the periapsis and apoapsis.

# Arguments
- `rp::Float64`: Periapsis, the radius from the central body to the nearest point on the orbital path.
- `ra::Float64`: Apoapsis, the radius from the central body to the farthest point on the orbital path.

# Returns
- `Float64`: Semi-major axis, the longest semi-diameter of an ellipse.
"""
function semi_major_axis(rp::Float64, ra::Float64)::Float64
    return (rp+ra)/2
end

"""
    orbital_period(sma::Float64, mu::Float64)::Float64

Calculate the orbital period of an object in an elliptical orbit.

# Arguments
- `sma::Float64`: Semi-major axis, the longest semi-diameter of an ellipse.
- `mu::Float64`: Gravitational parameter of the central body.

# Returns
- `Float64`: Orbital period in seconds.
"""
function orbital_period(sma::Float64, mu::Float64)::Float64
    return 2.0 * π * sqrt((sma^3)/mu)
end

"""
    eccentricity(rp::Float64, ra::Float64)::Float64

Calculate the eccentricity of an ellipse given the periapsis and apoapsis.

# Arguments
- `rp::Float64`: Periapsis, the radius from the central body to the nearest point on the orbital path.
- `ra::Float64`: Apoapsis, the radius from the central body to the farthest point on the orbital path.

# Returns
- `Float64`: Eccentricity of the ellipse.
"""
function eccentricity(rp::Float64, ra::Float64)::Float64
    semi_major_axis_val = semi_major_axis(rp, ra)
    ecc = 1 - rp/semi_major_axis_val
    return ecc
end

"""
    scaled_mu(sma::Float64, T::Float64, μ::Float64)::Float64

Calculate the scaled gravitational parameter of an object in an elliptical orbit.

# Arguments
- `sma::Float64`: Semi-major axis, the longest semi-diameter of an ellipse.
- `T::Float64`: Orbital period in seconds.
- `μ::Float64`: Gravitational parameter of the central body.

# Returns
- `Float64`: Scaled gravitational parameter.
"""
function scaled_mu(sma::Float64, T::Float64, μ::Float64)::Float64
    return ((T^2)/(sma^3)) * μ
end

"""
    instantaneous_orbital_speed(μ::Float64, radius::Float64, sma::Float64)::Float64

Calculate the instantaneous orbital speed of an object in an elliptical orbit.

# Arguments
- `μ::Float64`: Gravitational parameter of the central body.
- `radius::Float64`: Distance from the central body to the object.
- `sma::Float64`: Semi-major axis, the longest semi-diameter of an ellipse.

# Returns
- `Float64`: Instantaneous orbital speed.
"""
function instantaneous_orbital_speed(μ::Float64, radius::Float64, sma::Float64)::Float64
    return sqrt(μ*((2/radius)-(1/sma)))
end

"""
    orbital_mean_motion(T::Float64)

Calculate the mean motion (angular speed) required for an object to complete one orbit.

# Arguments
- `T::Float64`: Orbital period in seconds.

# Returns
- `Float64`: Mean motion in radians per unit time.
"""
function orbital_mean_motion(T::Float64)
    return 2 * π / T
end

"""
    get_CART_from_OSC(x_oe::AbstractArray{<:Real, 1}; degrees::Bool=false)

Return the cartesian state (position and velocity, ECI) given the corresponding set of osculating orbital elements.

# Arguments
- `x_oe::AbstractArray{<:Real, 1}`: Vector of osculating orbital elements in the following order: [a, e, i, Ω, ω, M].
- `degrees::Bool=false`: Set to `true` if the input angles are in degrees.

# Returns
- `Array{Float64, 1}`: Cartesian state vector [x, y, z, xdot, ydot, zdot].
"""
function get_CART_from_OSC(x_oe::AbstractArray{<:Real, 1}; degrees::Bool=false)
    return SD.sOSCtoCART(x_oe; use_degrees=degrees)
end

"""
    get_OSC_from_CART(x_oe::AbstractArray{<:Real, 1}; degrees::Bool=false)

Return the set of osculating orbital elements given the cartesian state (position and velocity, ECI).

# Arguments
- `x_oe::AbstractArray{<:Real, 1}`: Cartesian state vector [x, y, z, xdot, ydot, zdot].
- `degrees::Bool=false`: Set to `true` if the output angles should be in degrees.

# Returns
- `Array{Float64, 1}`: Vector of osculating orbital elements in the following order: [a, e, i, Ω, ω, M].
"""
function get_OSC_from_CART(x_oe::AbstractArray{<:Real, 1}; degrees::Bool=false)
    return SD.sCARTtoOSC(x_oe; use_degrees=degrees)
end
