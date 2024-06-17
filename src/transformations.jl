

function skew_symmetric(ω)
    """
    Returns the skew-symmetric matrix of a 3-dimensional vector v.
    """
    return @SMatrix [0 -ω[3] ω[2];
                    ω[3] 0 -ω[1];
                    -ω[2] ω[1] 0]
end


function L(q::Vector)
    """
    Left-multiply
    """
    L = zeros(eltype(q), 4, 4)
    L[1,1] = q[1]
    L[1,2:end] .= -q[2:end]
    L[2:end,1] .= q[2:end]
    L[2:end,2:end] .= q[1]*I + skew_symmetric(q[2:end])
    return L
end


function R(q::Vector)
    """
    Right-multiply
    """
    R = zeros(eltype(q), 4,4)
    R[1,1] = q[1]
    R[1,2:end] .= -q[2:end]
    R[2:end,1] .= q[2:end]
    R[2:end,2:end] .= q[1]*I - skew_symmetric(q[2:end])
    return R
end


function conj(q)
    """
    Inverse of a unit quaternion is its conjugate, i.e. same quaternion with a negated vector part 
    """
    qr = zeros(eltype(q),4)
    qr[1] = q[1]
    qr[2:end] .= -q[2:end]
    return qr
end


function mrp_from_quat(q::Vector)
    """
    Modified Rodriguez parameter inverse mapping
    """
    nq = norm(q)
    q = q/nq
    
    mrp = q[2:4] / (1 + q[1])

    # Shadow set switching 
    # Karlgaard, C. D., & Schaub, H. (2009). Nonsingular attitude filtering using modified Rodrigues parameters. 
    if norm(mrp) >= 1
        mrp = - mrp / (mrp' * mrp)
    end

    return  mrp
end

function quat_from_mrp(mrp::Vector)
    """
    Convert Modified Rodriguez parameters to quaternion
    """
    return normalize((1/(1+dot(mrp,mrp))) * vcat(1 - (dot(mrp,mrp)), 2*mrp))
end


function quat2axisangle(q)
    """
    Convert quaternion to axis-angle representation
    """
    axis = zeros(eltype(q),3)
    angle = 2 * acos(q[1])
    axis .= (q[2:end] / sqrt(1 - q[1]^2)) / norm(q[2:end] / sqrt(1 - q[1]^2))
    return axis * angle
end

function axisangle_from_quat(q)
    """
    Convert quaternion to axis-angle representation
    """
    axis = zeros(eltype(q),3)
    angle = 2 * acos(q[1])
    axis = (q[2:end] / sqrt(1 - q[1]^2)) / norm(q[2:end] / sqrt(1 - q[1]^2))
    return axis, angle
end


function quat_exp(q, eps=1e-8)
    """
    Exponential map of a quaternion
    """
    halfangle = norm(q[2:4])
    
    if halfangle < eps
        return q
    else
        c = cos(halfangle)
        s = sin(halfangle) / halfangle
        return exp(q[1])*([c, s * q[2], s * q[3], s * q[4]])
    end
end


function dcm_from_q(q::Vector)
    """
    Convert quaternion to direction cosine matrix (DCM)
    """
    q0,q1,q2,q3 = q / norm(q)

    # DCM
    return @SArray [(2*q1^2+2*q0^2-1)   2*(q1*q2 - q3*q0)   2*(q1*q3 + q2*q0);
    2*(q1*q2 + q3*q0)  (2*q2^2+2*q0^2-1)   2*(q2*q3 - q1*q0);
    2*(q1*q3 - q2*q0)   2*(q2*q3 + q1*q0)  (2*q3^2+2*q0^2-1)]
end


function dcm_from_mrp(p::Vector) 
    """
    Convert Modified Rodriguez parameters to direction cosine matrix (DCM)
    """
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)

    return @SArray[
        (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
        (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
        (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
        ]/den
end


function quat_from_dcm(R)
    """
    Convert direction cosine matrix (DCM) to quaternion
    """
    # can be + or -
    q0 = 0.5 * sqrt(R[1,1] + R[2,2] + R[3,3] + 1)
    q = [ q0, R[2,3]-R[3,2], R[3,1] - R[1,3], R[1,2] - R[2,1] ]
    q[2:4] .= q[2:4] / (4*q0)
    return q
end

function slerp(q1::Vector, q2::Vector, t::Real)
    """
    Spherical linear interpolation (SLERP) between two quaternions
    """
    # Inner product
    λ = q1' * q2

    if λ < 0
        q2 = -q2
        λ = - λ
    end

    # Quaternions nearly parallel so use linear interpolation
    if norm(1-λ) < 1e-2
        r = 1 - t
        s = t
    else
        # Calculate spherical linear interpolation factors 
        α = acos(λ)
        γ = 1 / sin(α)
        r = sin((1-t)*α)*γ
        s = sin(t*α)*γ
    end
    # Set the interpolated quaternion
    qf = r * q1 + s * q2
    # normalize
    qf = qf / norm(qf)

    return qf

end


function omega_mat(ω)
    ωx = ω[1]
    ωy = ω[2]
    ωz = ω[3]
  
    return [0 -ωx -ωy -ωz;
            ωx 0 ωz -ωy;
            ωy -ωz 0 ωx;
            ωz ωy -ωx 0]
end

function extract_ω_from_Omega(Omega)
    return Omega[2:4, 1]
end


function cayley_map(Φ)
    q = (1/sqrt(1+Φ' * Φ)) * vcat(1,Φ)
    return q / norm(q)
end

function inv_cayley_map(q)
    return q[2:4] / q[1]
end



function dcm_from_axisangle(axisangle::Vector)

    θ = norm(axisangle)
    u = axisangle / θ
    
    return @SMatrix [cos(θ)+(u[1]^2)*(1-cos(θ))   u[1]*u[2]*(1-cos(θ))-u[3]*sin(θ)   u[1]*u[3]*(1-cos(θ))+u[2]*sin(θ);
                    u[2]*u[1]*(1-cos(θ))+u[3]*sin(θ)    cos(θ)+(u[2]^2)*(1-cos(θ))      u[2]*u[3]*(1-cos(θ))-u[1]*sin(θ);
                    u[3]*u[1]*(1-cos(θ))-u[2]*sin(θ)    u[3]*u[2]*(1-cos(θ))+u[1]*sin(θ)    cos(θ)+(u[3]^2)*(1-cos(θ))]

end 