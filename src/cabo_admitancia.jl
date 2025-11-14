# ==============================================================================
# Transversal admittance Y per unit length in frequency-domain 
# of Single-Core (coaxial) and Pipe-Type cables.
# ==============================================================================

# Physical constants
const ε₀ = 8.854187817e-12  # Vacuum permittivity [F/m]

# ==============================================================================
# Calculation functions (individual formulas)
# ==============================================================================

"""
Elastance of a tubular capacitor.

# Parameters

- `radius_in`: internal radius of the capacitor \\[m\\].
- `radius_ext`: external radius of the capacitor \\[m\\].
- `epsr`: relative electric permittivity of the insulator \\[dimensionless\\].

# Returns

- `p`: The calculated elastance.
"""
function calc_tubular_elastance(
    radius_in::Real,
    radius_ext::Real,
    epsr::Real
)
    return log(radius_ext / radius_in) / (2 * π * ε₀ * epsr)
end


"""
Mutual elastance between two cables inside a tubular conductive pipe in relation 
to the pipe's internal surface.

# Parameters

- `distance_1`: distance from the center of the pipe to the center of conductor 1 \\[m\\].
- `distance_2`: distance from the center of the pipe to the center of conductor 2 \\[m\\].
- `theta`: angle between `distance_1` and `distance_2` \\[rad\\].
- `radius_1`: external radius of conductor 1 \\[m\\].
- `radius_2`: external radius of conductor 2 \\[m\\].
- `radius_in`: internal radius of the pipe \\[m\\].
- `epsr`: relative electric permittivity of the pipe's internal insulator \\[dimensionless\\].

# Returns

- `p`: The calculated mutual elastance.
"""
function calc_pipe_mutual_internal_elastance(
    distance_1::Real,
    distance_2::Real,
    theta::Real,
    radius_1::Real,
    radius_2::Real,
    radius_in::Real,
    epsr::Real
)
    r1 = radius_in
    di = distance_1
    dk = distance_2
    r_sq = r1^2 + 0.0im
    didk = di * dk
    TOL = 1e-12
    
    if abs(theta) < TOL && abs(di - dk) < TOL  # i == k
        Qik = log(r1 / radius_1 * (1 - di^2 / r_sq))
    else  # i ≠ k
        ln_1 = log((r_sq - didk * exp(1im * theta)) / r_sq)
        ln_2 = log((r_sq - didk * exp(-1im * theta)) / r_sq)
        Cn_sum = -0.5 * (ln_1 + ln_2)
        Qik = (
            log(radius_in / sqrt(di^2 + dk^2 - 2 * didk * cos(theta)))
            - Cn_sum
        )
    end
    
    return Qik / (2 * π * ε₀ * epsr)
end


"""
Self elastance of a cable inside a tubular conductive pipe in relation to 
the pipe's internal surface.

# Parameters

- `distance_1`: distance from the center of the pipe to the center of conductor 1 \\[m\\].
- `radius_1`: external radius of conductor 1 \\[m\\].
- `radius_in`: internal radius of the pipe \\[m\\].
- `epsr`: relative electric permittivity of the pipe's internal insulator \\[dimensionless\\].

# Returns

- `p`: The calculated self elastance.
"""
function calc_pipe_self_internal_elastance(
    distance_1::Real,
    radius_1::Real,
    radius_in::Real,
    epsr::Real
)
    return calc_pipe_mutual_internal_elastance(
        distance_1, distance_1, 0.0, radius_1, radius_1, radius_in, epsr
    )
end


# ==============================================================================
# Computation functions
# ==============================================================================

"""
Compute the shunt elastance matrix of a coaxial cable.

# Parameters

- `cable`: The coaxial cable object (CoaxialCable).

# Returns

- `P`: Elastance matrix.
"""
function comp_coaxial_cable_elastance(cable::CoaxialCable)
    Nc = length(cable.components)
    P = zeros(ComplexF64, Nc, Nc)

    for i in 1:Nc
        comp = cable.components[i]
        r1 = comp.radius_ext
        r2 = comp.radius_ext_insulator
        epsr = comp.epsr
        P[1:i, 1:i] .+= calc_tubular_elastance(r1, r2, epsr)
    end
    
    return P
end


"""
Compute the shunt elastance matrix of a pipe-type cable.

# Parameters

- `pipecable`: The pipe cable object.

# Returns

- `P`: Elastance matrix.
"""
function comp_pipe_cable_elastance(pipecable::PipeCable)
    Nc = count_conductors_cable(pipecable)
    P = zeros(ComplexF64, Nc, Nc)

    p_c = calc_tubular_elastance(
        pipecable.radius_ext,
        pipecable.radius_ext_insulator,
        pipecable.epsr_ext,
    )

    P[:, :] .= p_c

    N_inner = length(pipecable.cables)
    next_index_i = 1

    for i in 1:N_inner
        cable_i = pipecable.cables[i]
        x1 = cable_i.x - pipecable.x
        y1 = cable_i.y - pipecable.y
        distance_1 = hypot(x1, y1)
        radius_1 = outer_radius(cable_i)
        
        p_ii = calc_pipe_self_internal_elastance(
            distance_1, radius_1, pipecable.radius_in, pipecable.epsr_in
        )

        i1 = next_index_i
        i2 = i1 + count_conductors_cable(cable_i) - 1
        P[i1:i2, i1:i2] .+= p_ii
        next_index_i = i2 + 1

        next_index_k = next_index_i
        for k in (i+1):N_inner
            cable_k = pipecable.cables[k]
            x2 = cable_k.x - pipecable.x
            y2 = cable_k.y - pipecable.y
            distance_2 = hypot(x2, y2)
            radius_2 = outer_radius(cable_k)

            # Calculate angle between vectors
            dot_product = x1 * x2 + y1 * y2
            magnitude_product = distance_1 * distance_2
            TOL = 1e-12

            if abs(magnitude_product) < TOL
                theta = 0.0
            else
                theta = acos(dot_product / magnitude_product)
            end

            p_ik = calc_pipe_mutual_internal_elastance(
                distance_1,
                distance_2,
                theta,
                radius_1,
                radius_2,
                pipecable.radius_in,
                pipecable.epsr_in,
            )

            k1 = next_index_k
            k2 = k1 + count_conductors_cable(cable_k) - 1
            P[i1:i2, k1:k2] .+= p_ik
            P[k1:k2, i1:i2] .= P[i1:i2, k1:k2]
            next_index_k = k2 + 1
        end
    end

    return P
end


"""
Recursive function to compute elastance for any cable type.

A PipeCable and the components of a CoaxialCable receive an `_index`
attribute mapping the conductor to the matrix position.

# Parameters

- `P`: elastance matrix of the cable system that is modified in-place.
- `cable`: current cable being iterated upon.
- `start_index`: current index in the elastance matrix.

# Returns

- `Int`: index in the elastance matrix after the computation.
"""
function comp_cable_elastance_recursive!(
    P::Matrix{<:Complex{T}},
    cable::AbstractCable,
    start_index::Int
) where {T <: Real}
    n = count_conductors_cable(cable)
    end_index = start_index + n - 1
    P_k = @view P[start_index:end_index, start_index:end_index]

    if cable isa CoaxialCable
        P_k .+= comp_coaxial_cable_elastance(cable)
        for (i, comp) in enumerate(cable.components)
            comp._index = start_index + i - 1
        end

    elseif cable isa PipeCable
        P_k .+= comp_pipe_cable_elastance(cable)
        current_index = start_index
        cable._index = end_index
        for inner_cable in cable.cables
            current_index = comp_cable_elastance_recursive!(
                P, inner_cable, current_index
            )
        end

    else
        error("Unrecognized cable type `$(typeof(cable))`.")
    end

    return end_index + 1
end


"""
Compute the shunt elastance matrix of a cable system.

A PipeCable and the components of a CoaxialCable receive an `_index`
attribute mapping the conductor to the matrix position.

# Parameters

- `cable_system`: The outermost cable.

# Returns

- `P`: Elastance matrix.
"""
function comp_cable_system_elastance(cable_system::PipeCable)
    Nc = count_conductors_cable(cable_system)
    P = zeros(ComplexF64, Nc, Nc)
    current_index = 1
    comp_cable_elastance_recursive!(P, cable_system, current_index)
    return P
end


"""
Compute the shunt admittance matrix of a cable system at the given frequencies 
of interest in rad/s.

# Parameters

- `cable_system`: The outermost cable.
- `complex_frequencies`: Vector of complex angular frequencies `s = c + jω` \\[rad/s\\].

# Returns

- `Y`: Admittance matrix of shape (Nc, Nc, Nf).
"""
function comp_cable_system_admittance(
    cable_system::PipeCable,
    complex_frequencies::Vector{<:Complex{<:Real}}
)
    P = comp_cable_system_elastance(cable_system)
    P_inv = inv(P)
    P_inv = (P_inv + P_inv') / 2.0  # ensure symmetry

    Nc = size(P, 1)
    Nf = length(complex_frequencies)
    Y = zeros(ComplexF64, Nc, Nc, Nf)

    for (k, jω) in enumerate(complex_frequencies)
        Y[:, :, k] = jω * P_inv
    end

    return Y
end
