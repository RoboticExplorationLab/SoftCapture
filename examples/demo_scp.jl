
using Revise

using SoftCapture
using MeshCat

include("properties.jl")
include("plotter.jl")

##################### Initial conditions #####################

"""n_orbit = SC.Orbit()
n_orbit.region = "LEO"
n_orbit.sma = SC.semi_major_axis(SC.earth.R_LEO_ub, SC.earth.R_LEO_ub) # Semi-major axis [m]
n_orbit.cartesian_state = SC.get_CART_from_OSC([n_orbit.sma;zeros(5)])"""


sma = semi_major_axis(earth.R_LEO_ub, earth.R_LEO_ub) # Semi-major axis [m]

X0_chaser = [-6.725943874562646,19.98486574776072,16.780901068708612,0.008175460294039207,0.01100594682544459,0.011943924287753122] # [position, linear velocity]
X0_target = [1.0,0.0,0.0,0.0,-0.07845879017167497,0.03147562650127947,-0.021654991937582135] # [quaternion, angular velocity]

dt = 1.0
N = 175


RANDOM_INIT = false
if RANDOM_INIT
    X0_chaser = generate_relative_ic_safety_ellipse(sma; A0_min=10, A0_max=25, B0_min=15, B0_max=35)
    X0_target = sample_initial_tumble_state(deg2rad(10))

    # N = ... # Need to choose this one (or do a search)

end

##################### Parameters CVX #####################


chaser_body, target_body = generate_example_spacecraft_model()
target = Target(target_body, 0.2, X0_target, sma; max_time_propagation= N * dt)
Pc = ChaserTrajParams(chaser_body, dt, N, sma)

# Default values (shown for the example)
Pc.ζ = 5 #  Weight on the L1 control 
Pc.ψ = 750 # Weight on slack variable

Pc.goal_pos_tol = 0.35 # Tolerance on the goal position
Pc.soft_vel_tol = 0.03 # Tolerance on the velocity (soft capture)

Pc.α_colmin = 1.3 # Minimum scaling factor between the chaser and the target (must be greater than 1)
Pc.fov_cone_angle = deg2rad(60)# Field of view cone angle [rad]
Pc.docking_angle = deg2rad(30) # Docking cone angle [rad]
Pc.N_dock = 5 # Last N_dock knot points for the docking constraints 
Pc.p_max = 100.0 # Maximum position [m] (artifical)
Pc.T_max = 100.0 # Maximum thrust [N]
Pc.v_max = 1.5 # Maximum velocity [m/s]
Pc.ω_max = 0.2 # Maximum angular velocity [rad/s]


##################### Trajectory Optimization #####################

#status, prob, X, U, sU, ρ, J = solve_cvx(Pc, X0_chaser, target, get_MOSEK_solver())
#status, prob, dX, dU, sU, ρ, Sl, J = solve_cvx_correction(Pc, X, U, target, get_MOSEK_solver())

# status, X, U, Jcost, dcols = solve_trajectory_optimization(Pc, X0_chaser, target, get_MOSEK_solver)
status, X, U, Jcost, dcols = solve_trajectory_optimization(Pc, X0_chaser, target, get_ECOS_solver)

##################### Simulation #####################

dt_sim = 0.05
Xc, qc, Uc, Xt = simulate_trajectory(Pc, target, X0_chaser, X0_target, U, dt_sim)

##################### Visualization #####################

sim_N = size(Xc, 1)

vis = Visualizer()
set_background!(vis)
load_spacecraft_meshes!(vis)
configure_visualizer!(vis)


# Initial poses 
update_pose_chaser!(vis, X[1,:], qc[1,:], Uc[1,:])
update_pose_chaser!(vis, X[1,:], qc[1,:], zeros(3))

render(vis)
open(vis)
sleep(5.0)
print("Visualizer opened")

speed_up_sim_factor = 1.0 # change this parameters to accelerate the simulation
visualize_trajectories!(vis, Xc, qc, Uc, Xt, dt_sim, speed_up_sim_factor)


##################### Plotting #####################

plotter.plot_trajectory(Xc[:,1:3], dt_sim, "Position [m]")
plotter.plot_trajectory(Xc[:,4:6], dt_sim, "Velocity [m/s]")
plotter.plot_trajectory(Uc, dt_sim, "Thrust [N]")
