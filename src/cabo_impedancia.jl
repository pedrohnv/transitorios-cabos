# ==============================================================================
# Longitudinal impedance Z per unit length in frequency-domain 
# of Single-Core (coaxial) and Pipe-Type cables.
# ==============================================================================

using SpecialFunctions

# Physical constants
const μ₀ = 4e-7 * pi  # vacuum magnetic permeability [H/m]

# ==============================================================================
# Calculation functions (individual formulas)
# ==============================================================================

"""
Impedance of a tubular conductor considering the skin effect in its inner surface.

```math
\\frac{j\\omega \\mu_0 \\mu_{c1}}{2\\pi w_1} \\cdot \\frac{I_0(w_1) \\cdot K_1(w_2) + \\cdot I_1(w_2) \\cdot K_0(w_1)}{I_1(w_2) \\cdot K_1(w_1) - \\cdot I_1(w_1) \\cdot K_1(w_2)}
```

where ``w_1 = r_1 \\sqrt{j \\omega \\mu \\sigma}``, ``w_2 = r_2 \\sqrt{j \\omega \\mu \\sigma}``, ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

# Parameters

- `radius_in`: internal radius of the conductor \\[m\\]. Use zero if solid.
- `radius_ext`: external radius of the conductor \\[m\\].
- `rho_c`: conductivity of the conductor \\[Ω/m\\].
- `mur_pipe`: relative magnetic permeability of the conductor \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `simplified_formula`: use a simplified formula?
"""
function calc_inner_skin_effect_impedance(
    radius_in::Real,
    radius_ext::Real,
    rho_c::Real,
    mur_c::Real,
    complex_frequency::Complex{T},
    simplified_formula::Bool = false,
) where {T <: Real}
    iszero(radius_in) && return 0.0

    # Constants
    mu_c = μ₀ * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = sqrt(complex_frequency * mu_c * sigma_c)

    if simplified_formula
        cothTerm = coth(m * (radius_ext - radius_in))
        Z1 = (m / sigma_c) / (2 * pi * radius_in) * cothTerm
        Z2 = 1 / (2 * pi * radius_in * (radius_in + radius_ext) * sigma_c)
        zin = Z1 + Z2
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ext
        w_in = m * radius_in
        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besselix() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        N =
            sc *
            ( besselix(0, w_in)) *
            ( besselkx(1, w_out)) +
            ( besselkx(0, w_in)) *
            ( besselix(1, w_out))

        D =
            ( besselix(1, w_out)) *
            ( besselkx(1, w_in)) -
            sc *
            ( besselkx(1, w_out)) *
            ( besselix(1, w_in))
        # Final impedance calculation
        zin = (complex_frequency * mu_c / (2 * pi)) * (1 / w_in) * (N / D)
    end

    return zin
end


"""
Impedance of a tubular conductor considering the skin effect in its outer surface.

```math
\\frac{j\\omega \\mu_0 \\mu_{c1}}{2\\pi w_2} \\cdot \\frac{I_0(w_2) \\cdot K_1(w_1) + \\cdot I_1(w_1) \\cdot K_0(w_2)}{I_1(w_2) \\cdot K_1(w_1) - \\cdot I_1(w_1) \\cdot K_1(w_2)}
```

where ``w_1 = r_1 \\sqrt{j \\omega \\mu \\sigma}``, ``w_2 = r_2 \\sqrt{j \\omega \\mu \\sigma}``, ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

# Parameters

- `radius_in`: internal radius of the conductor \\[m\\]. Use zero if solid.
- `radius_ext`: external radius of the conductor \\[m\\].
- `rho_c`: conductivity of the conductor \\[Ω/m\\].
- `mur_pipe`: relative magnetic permeability of the conductor \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `simplified_formula`: use a simplified formula?
"""
function calc_outer_skin_effect_impedance(
    radius_in::Real,
    radius_ext::Real,
    rho_c::Real,
    mur_c::Real,
    complex_frequency::Complex{T},
    simplified_formula::Bool = false,
) where {T <: Real}
    # Constants
    mu_c = μ₀ * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = sqrt(complex_frequency * mu_c * sigma_c)

    if simplified_formula
        if radius_in < TOL
            cothTerm = coth(m * radius_ext * 0.733)
        else
            cothTerm = coth(m * (radius_ext - radius_in))
        end
        Z1 = (m / sigma_c) / (2 * pi * radius_ext) * cothTerm

        if radius_in < TOL
            Z2 = 0.3179 / (sigma_c * pi * radius_ext^2)
        else
            Z2 = 1 / (sigma_c * 2 * pi * radius_ext * (radius_in + radius_ext))
        end
        zin = Z1 + Z2
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ext
        w_in = m * radius_in
        if radius_in < TOL
            N = (besselix(0, w_out))
            D = (besselix(1, w_out))
        else
            s_in = exp(abs(real(w_in)) - w_out)
            s_out = exp(abs(real(w_out)) - w_in)
            sc = s_in / s_out  # Should be applied to all besseli() involving w_in

            # Bessel function terms with uncertainty handling using the macro
            N =
                (besselix(0, w_out)) *
                (besselkx(1, w_in)) +
                sc *
                (besselkx(0, w_out)) *
                (besselix(1, w_in))

            D =
                (besselix(1, w_out)) *
                (besselkx(1, w_in)) -
                sc *
                (besselkx(1, w_out)) *
                (besselix(1, w_in))
        end
        # Final impedance calculation
        zin = (complex_frequency * mu_c / (2 * pi)) * (1 / w_out) * (N / D)
    end
    return zin
end


"""
Mutual impedance between the inner and outer surfaces of tubular conductor considering the skin effect.

```math
\\frac{1}{2\\pi r_1 r_2 \\sigma} \\cdot \\frac{1}{I_1(w_2) \\cdot K_1(w_1) - \\cdot I_1(w_1) \\cdot K_1(w_2)}
```

where ``w_1 = r_1 \\sqrt{j \\omega \\mu \\sigma}``, ``w_2 = r_2 \\sqrt{j \\omega \\mu \\sigma}``, ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

# Parameters

- `radius_in`: internal radius of the conductor \\[m\\]. Use zero if solid.
- `radius_ext`: external radius of the conductor \\[m\\].
- `rho_c`: conductivity of the conductor \\[Ω/m\\].
- `mur_pipe`: relative magnetic permeability of the conductor \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `simplified_formula`: use a simplified formula?
"""
function calc_mutual_skin_effect_impedance(
    radius_in::Real,
    radius_ext::Real,
    rho_c::Real,
    mur_c::Real,
    complex_frequency::Complex{T},
    simplified_formula::Bool = false,
) where {T <: Real}
    if iszero(radius_in)
        return calc_outer_skin_effect_impedance(
            radius_in,
            radius_ext,
            rho_c,
            mur_c,
            complex_frequency,
        )
    end

    # Constants
    mu_c = μ₀ * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = sqrt(complex_frequency * mu_c * sigma_c)

    if simplified_formula
        cschTerm = csch(m * (radius_ext - radius_in))
        zm = m / (sigma_c * pi * (radius_in + radius_ext)) * cschTerm
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ext
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besselix() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        D =
            ( besselix(1, w_out)) *
            ( besselkx(1, w_in)) -
            sc *
            ( besselkx(1, w_out)) *
            ( besselix(1, w_in))

        # Final mutual impedance calculation
        zm = 1 / (2 * pi * radius_ext * radius_in * sigma_c * D * s_out)
    end

    return zm
end


"""
Mutual impedance between the inner and outer surfaces of a tubular capacitor considering the skin effect.

```math
\\frac{j \\omega \\mu}{2\\pi} \\cdot \\ln \\frac{r_2}{r_1}
```

where ``w_1 = r_1 \\sqrt{j \\omega \\mu \\sigma}``, ``w_2 = r_2 \\sqrt{j \\omega \\mu \\sigma}``, ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

# Parameters

- `radius_in`: internal radius of the conductor \\[m\\]. Use zero if solid.
- `radius_ext`: external radius of the conductor \\[m\\].
- `rho_c`: conductivity of the conductor \\[Ω/m\\].
- `mur_pipe`: relative magnetic permeability of the conductor \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `simplified_formula`: use a simplified formula?
"""
function calc_tubular_capacitor_impedance(
    radius_in::Real,
    radius_ext::Real,
    mur_ins::Real,
    complex_frequency::Complex{T},
) where {T <: Real}
    if iszero(radius_in)
        return calc_outer_skin_effect_impedance(
            radius_in,
            radius_ext,
            rho_c,
            mur_c,
            complex_frequency,
        )
    end
    return complex_frequency * μ₀ * mur_ins * log(radius_ext / radius_in) / (2 * pi)
end


"""
Mutual impedance between two cables inside a tubular conductive pipe in relation to the pipe's internal surface.

```math
\\frac{j\\omega \\mu_0}{2\\pi} \\cdot \\left[ \\frac{\\mu_p K_0(w_1)}{w_1 K_1(w_1)} + Q_{ik} + 2 \\mu_p \\sum_{n=1}^{\\infty} \\frac{C_n}{n (1+\\mu_p) + w_1 K_{n-1}(w_1)/K_n(w_1)} \\right]
```

where

```math
Q_{ik} = \\ln \\left[ \\frac{r_{p1}}{\\sqrt{d_i^2 + d_k^2 - 2 d_i d_k \\cos\\theta_{ik}}} \\right] - \\sum_{n=1}^{\\infty} \\frac{C_n}{n}
```

```math
C_n = \\left( \\frac{d_i d_k}{r_{p1}^2} \\right)^n \\cdot \\cos(n \\theta_{ik})
```

```math
w_1 = r_{p1} \\sqrt{j\\omega \\mu_0 \\mu_p \\sigma_p}
```

and ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

Note that ``\\sum_{n=1}^{\\infty} \\frac{C_n}{n}`` has a closed formula solution:

```math
\\sum_{n=1}^{\\infty} \\frac{C_n}{n} = \\frac{-1}{2} \\left[ \\ln\\left( \\frac{r_{p1}^2 - d_i  d_k  \\exp(j\\theta_{ik})}{r_{p1}^2} \\right) + \\ln\\left( \\frac{r_{p1}^2 - d_i  d_k \\exp(-j\\theta_{ik})}{r_{p1}^2} \\right) \\right]
```

# Parameters

- `distance_1`: distance from the center of the pipe to the center of conductor 1 \\[m\\].
- `distance_2`: distance from the center of the pipe to the center of conductor 2 \\[m\\].
- `theta`: angle between `distance_1` and `distance_2` \\[rad\\].
- `radius_1`: external radius of conductor 1 \\[m\\].
- `radius_2`: external radius of conductor 2 \\[m\\].
- `radius_in`: internal radius of the pipe \\[m\\].
- `sigma_c`: conductivity of the pipe \\[S/m\\].
- `mur_pipe`: relative magnetic permeability of the pipe \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `maxiter`: maximum number of terms to consider in the infinite series.
- `tol`: tolerance to check for convergence of the infinite series.
"""
function calc_pipe_mutual_internal_impedance(
    distance_1::Real,
    distance_2::Real,
    theta::Real,
    radius_1::Real,
    radius_2::Real,
    radius_in::Real,
    rho_c::Real,
    mur_pipe::Real,
    complex_frequency::Complex{T},
    maxiter::Int = 100,
    tol::Real = 1e-12,
) where {T <: Real}
    sigma_c = 1.0 / rho_c

    r1 = radius_in
    di = distance_1
    dk = distance_2
    
    r² = r1^2 + 0.0im
    didk = di * dk
    w = r1 * sqrt(complex_frequency * μ₀ * mur_pipe * sigma_c)

    K0 =  besselkx(0, w)
    K1 =  besselkx(1, w)
    term1 = mur_pipe * K0 / (w * K1)
    
    if iszero(theta) && isapprox(di, dk)  # i == k
        Qik = log(r1 / radius_1 * (1 - di^2 / r²))
    else  # i ≠ k
        ln_1 = log( (r² - didk * exp(1im * theta)) / (r²) )
        ln_2 = log( (r² - didk * exp(-1im * theta)) / (r²) )
        Cn_sum = -0.5 * (ln_1 + ln_2)
        Qik = log(radius_in / sqrt(di^2 + dk^2 - 2 * didk * cos(theta))) - Cn_sum
    end

    Cn = n -> (didk / r²)^n * cos(n * theta)
    Dn = (n, K_prev, K_n) -> (n * (1 + mur_pipe) + w * K_prev / K_n)
    K_prev = K0
    K_n = K1
    inf_series = Cn(1) / Dn(1, K_prev, K_n)
    converged = false
    overflowing = false
    for n in 2:maxiter
        if !overflowing
            try
                # TODO is there a recursive relation we can use to avoid calculating Kₙ(w) every iteration? See for a start:
                # M. Onoe (1955) Formulae and Tables, The Modified Quotients of Cylinder Functions. Technical report Technical Report UDC 517.564.3:518.25, Vol. 4, Report of the Institute of Industrial Science, University of Tokyo, Institute of Industrial Science, Chiba City, Japan. 
                # https://www.ams.org/journals/mcom/1956-10-053/S0025-5718-1956-0076107-X/S0025-5718-1956-0076107-X.pdf
                K_prev = K_n
                K_n =  besselkx(n, w)
                term_n = Cn(n) / Dn(n, K_prev, K_n)
            catch err
                if err.id == 4  # overflow, then `Kₙ(w)/Kₙ₊₁(w) = 0`
                    overflowing = true
                else
                    rethrow()
                end
            end
        end

        if overflowing
            term_n = Cn(n) / (n * (1 + mur_pipe))
        end

        inf_series += term_n
        if n > 5 && abs(term_n) < tol * abs(inf_series)
            converged = true
            break
        elseif n == maxiter && !converged
            @warn "calc_pipe_mutual_internal_impedance did not converge to specified tolerance."
            @show abs(term_n) / abs(inf_series)
        end
    end
    return complex_frequency * μ₀ / (2π) * (term1 + Qik + inf_series * 2 * mur_pipe)
end


"""
Self impedance of a cable inside a tubular conductive pipe in relation to the pipe's internal surface.

```math
\\frac{j\\omega \\mu_0}{2\\pi} \\cdot \\left[ \\frac{\\mu_p K_0(w_1)}{w_1 K_1(w_1)} + Q_{ii} + 2 \\mu_p \\sum_{n=1}^{\\infty} \\frac{C_n}{n (1+\\mu_p) + w_1 K_{n-1}(w_1)/K_n(w_1)} \\right]
```

where

```math
Q_{ii} = \\ln \\left[ \\frac{r_{p1}}{r_i} \\cdot \\left( 1 - \\left(\\frac{d_i}{r_{p1}}\\right)^2 \\right) \\right]
```

```math
C_n = \\left( \\frac{d_i d_k}{r_{p1}^2} \\right)^n \\cdot \\cos(n \\theta_{ik})
```

```math
w_1 = r_{p1} \\sqrt{j\\omega \\mu_0 \\mu_p \\sigma_p}
```

and ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of the first and second kind, respectively, of order ``\\nu``.

# Parameters

- `distance_1`: distance from the center of the pipe to the center of conductor 1 \\[m\\].
- `radius_1`: external radius of conductor 1 \\[m\\].
- `radius_in`: internal radius of the pipe \\[m\\].
- `sigma_c`: conductivity of the pipe \\[S/m\\].
- `mur_pipe`: relative magnetic permeability of the pipe \\[dimensionless\\].
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `maxiter`: maximum number of terms to consider in the infinite series.
- `tol`: tolerance to check for convergence of the infinite series.
"""
function calc_pipe_self_internal_impedance(
    distance_1::Real,
    radius_1::Real,
    radius_in::Real,
    rho_c::Real,
    mur_pipe::Real,
    complex_frequency::Complex{T},
    maxiter::Int = 100,
    tol::Real = 1e-12,
) where {T <: Real}
    return calc_pipe_mutual_internal_impedance(
        distance_1,
        distance_1,
        zero(T),
        radius_1,
        radius_1,
        radius_in,
        rho_c,
        mur_pipe,
        complex_frequency,
        maxiter,
        tol,
    )
end


# ==============================================================================
# Computation functions
# ==============================================================================

"""Compute the series impedance matrix of a coaxial cable at the given complex frequency in rad/s."""
function comp_coaxial_cable_impedance(
    cable::CoaxialCable,
    complex_frequency::Complex{T},
) where {T <: Real}
    # number of conductors
    Nc = length(cable.design_data.components)
    Z = zeros(Complex{T}, Nc, Nc)

    # There is a recursive relation:
    # -ΔV[i + 1]/Δx = z_loop * I_loop[i + 1] - z_mutual_io * I_loop[i]
    for i = 1:Nc
        comp = cable.design_data.components[i]
        rho_c = comp.conductor_props.rho
        mur_c = comp.conductor_props.mu_r
        radius_in = comp.conductor_group.radius_in
        radius_ext = comp.conductor_group.radius_ext
        radius_ext_insulator = comp.insulator_group.radius_ext
        mur_d = comp.insulator_props.mu_r
        
        z_inner = calc_inner_skin_effect_impedance(
            radius_in,
            radius_ext,
            rho_c,
            mur_c,
            complex_frequency,
        )

        z_outer = calc_outer_skin_effect_impedance(
            radius_in,
            radius_ext,
            rho_c,
            mur_c,
            complex_frequency,
        )

        if i > 1
            z_mutual_io = calc_mutual_skin_effect_impedance(
                radius_in,
                radius_ext,
                rho_c,
                mur_c,
                complex_frequency,
            )
        else
            z_mutual_io = 0
        end

        z_mutual_d = calc_tubular_capacitor_impedance(
            radius_ext,
            radius_ext_insulator,
            mur_d,
            complex_frequency,
        )

        if i > 1
            Z[1:(i - 1), 1:(i - 1)] .+= z_outer + z_mutual_d + z_inner - 2 * z_mutual_io
            Z[i, 1:(i - 1)] .+= z_outer + z_mutual_d - z_mutual_io
            Z[1:(i - 1), i] .+= z_outer + z_mutual_d - z_mutual_io
        end
        Z[i, i] += z_outer + z_mutual_d
    end
    return Z
end


"""Compute the series impedance matrix of a pipe-type cable at the given
complex frequency in rad/s.

This function considers only the outermost layer of the inner conductors.    
It does NOT compute the internal impedance of the inner conductors.
"""
function comp_pipe_cable_impedance(
    pipecable::PipeCable,
    complex_frequency::Complex{T},
) where {T <: Real}
    Nc = count_conductors_cable(pipecable)
    Z = zeros(Complex{T}, Nc, Nc)

    z_inner = calc_inner_skin_effect_impedance(
        pipecable.radius_in,
        pipecable.radius_ext,
        pipecable.rho_c,
        pipecable.mur_c,
        complex_frequency,
    )

    z_outer = calc_outer_skin_effect_impedance(
        pipecable.radius_in,
        pipecable.radius_ext,
        pipecable.rho_c,
        pipecable.mur_c,
        complex_frequency,
    )

    z_mutual_io = calc_mutual_skin_effect_impedance(
        pipecable.radius_in,
        pipecable.radius_ext,
        pipecable.rho_c,
        pipecable.mur_c,
        complex_frequency,
    )

    z_mutual_d = calc_tubular_capacitor_impedance(
        pipecable.radius_ext,
        pipecable.radius_ext_insulator,
        pipecable.mur_d_ext,
        complex_frequency,
    )

    Z[1:Nc-1, 1:Nc-1] .= z_outer + z_mutual_d - 2 * z_mutual_io + z_inner
    Z[1:Nc-1, end] .= z_outer + z_mutual_d - z_mutual_io
    Z[end, 1:Nc-1] .= z_outer + z_mutual_d - z_mutual_io
    Z[end, end] = z_outer + z_mutual_d

    N_inner = length(pipecable.cables)
    next_index_i = 1
    for i in 1:N_inner
        cable_i = pipecable.cables[i]
        x1 = cable_i.x - pipecable.x
        y1 = cable_i.y - pipecable.y
        distance_1 = sqrt(x1^2 + y1^2)

        if cable_i isa PipeCable
            radius_1 = cable_i.radius_ext_insulator
        else
            radius_1 = cable_i.components[end].radius_ext_insulator
        end

        zm_ii = calc_pipe_self_internal_impedance(
            distance_1,
            radius_1,
            pipecable.radius_in,
            pipecable.rho_c,
            pipecable.mur_c,
            complex_frequency,
        )

        i1 = next_index_i
        i2 = i1 + count_conductors_cable(cable_i) - 1
        Z[i1:i2, i1:i2] .+= zm_ii
        next_index_i = i2 + 1
        
        next_index_k = next_index_i
        for k in (i+1):N_inner
            cable_k = pipecable.cables[k]
            x2 = cable_k.x - pipecable.x
            y2 = cable_k.y - pipecable.y
            distance_2 = sqrt(x2^2 + y2^2)

            if cable_k isa PipeCable
                radius_2 = cable_k.radius_ext_insulator
            else
                radius_2 = cable_k.components[end].radius_ext_insulator
            end

            theta = acos((x1 * x2 + y1 * y2) / (distance_1 * distance_2))

            zm_ik = calc_pipe_mutual_internal_impedance(
                distance_1,
                distance_2,
                theta,
                radius_1,
                radius_2,
                pipecable.radius_in,
                pipecable.rho_c,
                pipecable.mur_c,
                complex_frequency,
            )

            k1 = next_index_k
            k2 = k1 + count_conductors_cable(cable_k) - 1
            Z[i1:i2, k1:k2] .+= zm_ik
            Z[k1:k2, i1:i2] .= Z[i1:i2, k1:k2]
            next_index_k = k2 + 1
        end
    end
    return Z
end


"""Recursive function to compute impedance for any cable type.

# Parameters
- `Z`: impedance matrix of the cable system that is modified in-place.
- `cable`: current cable being iterated upon.
- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `start_index`: current index in the impedance matrix.

# Returns
- `next_index`: index in the impedance matrix after the computation.
"""
function comp_cable_impedance_recursive!(
    Z::Matrix{Complex{T}},
    cable::AbstractCable,
    complex_frequency::Complex{T},
    start_index::Int,
) where {T <: Real}
    n = count_conductors_cable(cable)
    end_index = start_index + n - 1
    if cable isa CoaxialCable
        Z_coaxial = comp_coaxial_cable_impedance(cable, complex_frequency)
        Z[start_index:end_index, start_index:end_index] += Z_coaxial

    elseif cable isa PipeCable
        Z_pipe = comp_pipe_cable_impedance(cable, complex_frequency)
        Z[start_index:end_index, start_index:end_index] += Z_pipe

        current_index = start_index
        for inner_cable in cable.cables
            current_index = comp_cable_impedance_recursive!(
                Z, inner_cable, complex_frequency, current_index)
        end

    else
        @error "Unrecognized cable type `$(typeof(cable))`."
    end

    return start_index + n
end


"""Compute the series impedance matrix of a cable system at the given complex frequency in rad/s."""
function comp_cable_system_impedance(
    cable_system::PipeCable,
    complex_frequency::Complex{T},
    sigma_mar::Real = 5.0,
    epsr_mar::Real = 81.0,
) where {T <: Real}
    Nc = count_conductors_cable(cable_system)
    Z = zeros(Complex{T}, Nc, Nc)
    current_index = 1
    for cable in cable_system.cables
        current_index = comp_cable_impedance_recursive!(Z, cable, complex_frequency, current_index)
    end

    if !isnothing(sigma_mar)
        Z[:, :, k] += cZmar(jω, rca, sigma_mar, epsr_mar)
    end

    return Z
end


"""Impedância de retorno pelo mar.

# Parâmetros

- `freq_s`: Frequência angular complexa no formato `c + jw` [rad/s].
- `rca`: Raio da capa externa da armadura [m].
- `sigma`: Permissividade elétrica do mar [S/m]. O padrão é `5.0`.
    Use `nothing` para ignorar o mar.
- `epsr`: Permissividade elétrica do mar. O padrão é `81.0`.

# Retorna

- `Z`: impedância do mar [Ω/m].
"""
function cZmar(
    freq_s::Complex{T},
    rca::Real,
    sigma::Real = 5.0,
    epsr::Real = 81.0,
) where {T <: Real}
    if isnothing(sigma) || isnothing(epsr)
        return 0.0
    end
    jw = freq_s
    rho = 1.0 / sigma
    eta = sqrt(jw * mu_0 * (sigma + jw * epsilon_0 * epsr))
    Zmar = eta * rho / (2 * pi * rca) * besselkx(0, eta * rca) / besselkx(1, eta * rca)
    return Zmar
end
