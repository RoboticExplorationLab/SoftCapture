
function set_background!(vis, display_earth=true)
    mc.setprop!(vis["/Background"], "top_color", Gray(0.0))
    mc.setprop!(vis["/Background"], "bottom_color", Gray(0.0))

    if display_earth
        # Plot Earth in MeshCat
        length_scale = 1.5e-6
        earth_image = mc.PngImage("assets/earth.png")
        earth_texture = mc.Texture(image=earth_image)
        earth_material = mc.MeshLambertMaterial(map=earth_texture)
        earth = HyperSphere(Point(0.0,0.0,0.0), 6.7310*3e6*length_scale)
        mc.setobject!(vis["earth"], earth, earth_material)
        mc.settransform!(vis["earth"],  Translation([-30,-30,0.0]) ∘ LinearMap(AngleAxis(0.2*pi, 0, -0.5, 1)) ∘ LinearMap(AngleAxis(pi/2, 1, 0, 0)))
    end        

end


function load_spacecraft_meshes!(vis, chaser_mesh_path = "assets/scaled_chaser.STL", target_mesh_path = "assets/scaled_target.STL", add_port = true)

    chaser_mesh = load(chaser_mesh_path)
    target_mesh = load(target_mesh_path)
    mc.setobject!(vis["chaser"], chaser_mesh, mc.MeshPhongMaterial(color = RGBA(1,0.94,0.86,1)))
    mc.setobject!(vis["target"], target_mesh, mc.MeshPhongMaterial(color = RGBA(1.0, 0.71, 0.76,1)))

    if add_port
        ## Cylinder port on spacecrafts
        cyl_target = Cylinder(Point(0,0.02*3,0.35), Point(0,0.02*3,0.43), 0.275)
        cyl_target = Cylinder(Point(0,0.02*3,0.35*3), Point(0,0.02*3,0.43*3), 3*0.28)
        mc.setobject!(vis["target/cyl"], cyl_target, mc.MeshPhongMaterial(color = RGBA(1,0,0,0.5)))

        cyl_chaser = Cylinder(Point(0,0.0,0.4), Point(0,0.0,0.45), 0.15)
        cyl_chaser = Cylinder(Point(0,0.0,0.4*3), Point(0,0.0,0.45*3), 3*0.15)
        mc.setobject!(vis["chaser/cyl"], cyl_chaser, mc.MeshPhongMaterial(color = RGBA(1,0.5,0,1)))
    end

end


function configure_visualizer!(vis;  hull_on=false, grid_on=false, axes_one=true, pointing_arrow_on=false, velocity_arrow_on=true, thrust_arrow_on=true)

    mc.setprop!(vis["/Grid"], "visible", grid_on)
    mc.setprop!(vis["/Axes"], "visible", axes_one)
    mc.setvisible!(vis[:chaser_hull], hull_on)
    mc.setvisible!(vis[:target_hull], hull_on)
    
    mc.setvisible!(vis[:chaser_pointing], pointing_arrow_on)
    mc.setvisible!(vis[:chaser_vel], velocity_arrow_on)
    mc.setvisible!(vis[:chaser_thrust], thrust_arrow_on)
end



## Poiting arrow 
function update_chaser_pointing!(vis, p0::Vector)
    origin = Point(p0[1], p0[2], p0[3])
    arrow_end = origin + (- p0 / norm(p0)).* 0.3*3
    arrow_end = Point(arrow_end[1], arrow_end[2], arrow_end[3]) 
    cone_end = origin + (- p0 / norm(p0)).* 0.45*3
    cone_end = Point(cone_end[1], cone_end[2], cone_end[3]) 

    chaser_pointing = Cylinder(origin, arrow_end, 0.10)
    mc.setobject!(vis["chaser_pointing"], chaser_pointing, mc.MeshPhongMaterial(color = RGBA(1,0.94,0.86,1)))

    cone = mc.Cone(arrow_end, cone_end, 0.2)
    mc.setobject!(vis["chaser_pointing/cone"], cone, mc.MeshPhongMaterial(color = RGBA(1,0.5,0,1)))
end

## Velocity arrow
function update_chaser_velocity!(vis, p::Vector, v::Vector)
    origin = Point(p[1], p[2], p[3])
    arrow_end = origin + (v / 1.0) * 5 
    arrow_end = Point(arrow_end[1], arrow_end[2], arrow_end[3]) 

    cone_end = origin + (v / 1.0) * 7.5
    cone_end = Point(cone_end[1], cone_end[2], cone_end[3]) 

    chaser_vel = Cylinder(origin, arrow_end, 0.05)
    mc.setobject!(vis["chaser_vel"], chaser_vel, mc.MeshPhongMaterial(color = RGBA(0,0,1,1)))

    cone = mc.Cone(arrow_end, cone_end, 0.1)
    mc.setobject!(vis["chaser_vel/cone"], cone, mc.MeshPhongMaterial(color = RGBA(0,0,1,1)))
end


## Thrust arrow
function update_chaser_thrust!(vis, p::Vector, T::Vector)
    if norm(T) > 1e-6  
        origin = Point(p[1], p[2], p[3])
        arrow_end = origin + (T / 50) * 5 
        arrow_end = Point(arrow_end[1], arrow_end[2], arrow_end[3]) 

        cone_end = origin + (T / 50) * 7.5
        cone_end = Point(cone_end[1], cone_end[2], cone_end[3]) 

        chaser_vel = Cylinder(origin, arrow_end, 0.1)
        mc.setobject!(vis["chaser_thrust"], chaser_vel, mc.MeshPhongMaterial(color = RGBA(1,0,0,1)))

        cone = mc.Cone(arrow_end, cone_end, 0.2)
        mc.setobject!(vis["chaser_thrust/cone"], cone, mc.MeshPhongMaterial(color = RGBA(1,0,0,1)))
    else
        delete!(vis["chaser_thrust"])
    end
end


function update_pose_chaser!(vis, Xch::Vector, qch::Vector, U::Vector)

    mc.settransform!(vis["chaser"], Translation(Xch[1:3]) ∘ LinearMap(QuatRotation(qch)))

    # If set to no visible, won't display anyway
    update_chaser_pointing!(vis, Xch[1:3])
    update_chaser_velocity!(vis, Xch[1:3], Xch[4:6])
    update_chaser_thrust!(vis, Xch[1:3], U[1:3])
end



function update_pose_target!(vis, Xta::Vector)
    mc.settransform!(vis["target"], LinearMap(QuatRotation(Xta[1:4])))
end

function visualize_trajectories!(vis, X, qc, U, Xt, dt, speed_up=1.0)
    sim_N = size(X, 1)
    for i in 1:sim_N
        if i != sim_N
            update_pose_chaser!(vis, X[i,:], qc[i,:], U[i,:])
        else
            update_pose_chaser!(vis, X[i,:], qc[i,:], zeros(3))
        end
        update_pose_target!(vis, Xt[i,:])
        sleep(dt / speed_up)
    end
end