"""
Provides example of 
- geometric and inertial properties, includng the definition of convex hulls of both bodies.
- Solver functions for MOSEK and ECOS.

The algorithm requires the use of DifferentiableCollisions.jl to compute collision avoidance constraints.

"""

using StaticArrays
import DifferentiableCollisions as dc
import MathOptInterface as MOI
using Mosek
using MosekTools
using ECOS


function generate_example_spacecraft_model()

    # Generate inertial and geometrical properties for the target and chaser body

    ############### Target ###############
    
    # Inertia matrix 
    target_mass = 1000.0
    tIxx = 5.89056
    tIyy = 11.4462
    tIzz = 11.5365
    tIxy= -1.37725e-7
    tIxz= 2.97917e-8
    tIyz= 0.233516
    Jt = [tIxx tIxy tIxz; tIxy tIyy tIyz; tIxz tIyz tIzz]

    # Defines the mating point on body frame w.r.t geometrical center 
    capture_port_target = [0.0, 0.0, 1.35] 

    # Defines the convex hull for collision avoidance 
    At = SMatrix{7,3}([0.0 0 1;
                    0.0 0 -1;
                    0.0 1 0;
                    -1 3 0;
                    -3 -3.6 0;
                    1 3 0;
                    3 -3.6 0;
                    ])
    bt = SVector{7}(3*[0.4,
                    0.4,
                    0.4,
                    1.51,
                    4.5,
                    1.51,
                    4.5])
    st = dc.Polytope(At, bt)

    target = BodyParams(target_mass, Jt, capture_port_target, st)

    ############### Chaser ###############
    
    # Inertia matrix 
    chaser_mass = 1500
    cIxx = 16.072 
    cIyy = 13.3014 
    cIzz = 16.0409
    cIxy= 3.07409e-11
    cIxz= 4.25544e-8
    cIyz= 4.49124e-10
    Jc = [cIxx cIxy cIxz; cIxy cIyy cIyz; cIxz cIyz cIzz]

    # Defines the mating point on body frame w.r.t geometrical center 
    capture_port_chaser = [0.0, 0.0, 1.35]

    # Defines the convex hull for collision avoidance 
    Ac = SMatrix{8,3}([0.0 0 1;
                0.0 0 -1;
                1 0.0 0;
                -1 0.0 0;
                -2 1 0;
                -2 -1 0;
                2 1 0;
                2 -1 0;
                ])
    bc = SVector{8}(3*[0.4033,
                    0.4,
                    0.4033,
                    0.4033,
                    1.2,
                    1.2,
                    1.2,
                    1.2])
    sc = dc.Polytope(Ac, bc)

    chaser = BodyParams(chaser_mass, Jc, capture_port_chaser, sc)


    return chaser, target

end


function get_MOSEK_solver()
    # https://docs.mosek.com/latest/juliaapi/parameters.html#mosek.dparam.intpnt_co_tol_dfeas
    return () -> begin
        optimizer = Mosek.Optimizer()
        MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_PFEAS"), 1e-6) # Primal feasibility tolerance used by the interior-point optimizer for linear problems.
        MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1e-6) # Dual feasibility tolerance used by the interior-point optimizer for linear problems.
        return optimizer
    end
end

function get_ECOS_solver()
    return () -> begin
        optimizer = ECOS.Optimizer()
        MOI.set(optimizer, MOI.RawOptimizerAttribute("feastol"), 1e-6) # primal/dual infeasibility tolerance
        return optimizer
    end
end

