__precompile__(true)
module SoftCapture

using LinearAlgebra
using StaticArrays
using ExponentialAction
import SatelliteDynamics as SD

import Convex as cvx 
import ForwardDiff as FD
import DifferentiableCollisions as dc

using Random
using Distributions
using Interpolations

using PyCall
using Printf

import MeshCat as mc
using CoordinateTransformations
using Rotations
using GeometryBasics
using Colors 
using StaticArrays 
using FileIO



export uniquify, 
        mat_from_vec, 
        vec_from_mat, 
        has_nonzero_elements



export  skew_symmetric, 
        L, 
        R, 
        conj, 
        mrp_from_quat, 
        quat_from_mrp, 
        quat2angleaxis,
        axisangle_from_quat,
        quat_exp,
        dcm_from_q,
        dcm_from_mrp,
        quat_from_dcm,
        slerp,
        omega_mat,
        cayley_map,
        inv_cayley_map,
        dcm_from_axisangle

export earth, 
    LEO_lb,
    LEO_ub,
    MEO,
    GSO,
    Orbit,
    sample_new_orbit,
    semi_major_axis, 
    orbital_period, 
    eccentricity,
    scaled_mu,
    instantaneous_orbital_speed,
    orbital_mean_motion,
    get_CART_from_OSC,
    get_OSC_from_CART

export rk4,
    euler_rotational_dynamics,
    clohessy_wiltshire,
    clohessy_wiltshire_matrices,
    discretized_clohessy_wiltshire,
    quat_kinematics,
    attitude_quat_rbd,
    rk4mk_attitude,
    target_dynamics,
    rk4mk_target,
    nonlinear_relative_keplerian_dynamics,
    rk4_chaser,
    scale_state,
    unscale_state

export BodyParams,
    StateIdx,
    ChaserTrajParams,
    Target,
    predict_target_motion,
    repropagate_target!,
    get_predicted_target_trajectory,
    generate_relative_ic_safety_ellipse,
    sample_initial_tumble_state


export spacecraft_dynamics,
    solve_cvx,
    solve_cvx_correction,
    solve_trajectory_optimization,
    collision_avoidacne_dcol,
    get_attitude_pointing,
    get_next_fov_attitude,
    zero_order_hold_interpolation,
    simulate_trajectory

export set_background!,
    load_spacecraft_meshes!,
    configure_visualizer!,
    update_chaser_pointing!,
    update_chaser_velocity!,
    update_chaser_thrust!,
    update_pose_chaser!,
    update_pose_target!,
    visualize_trajectories!

export Visualizer, open, render


include("utils.jl")
include("transformations.jl")
include("astrodynamics.jl")
include("dynamics.jl")
include("configuration.jl")
include("cvx_formulation.jl")
include("visualizer.jl")


end # module