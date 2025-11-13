# -*- coding: utf-8 -*-
"""Longitudinal impedance Z per unit length in frequency-domain of Single-Core (coaxial) and Pipe-Type cables."""

from typing import Optional
import numpy as np
from scipy.special import ive, kve
import warnings
from scipy.constants import pi, mu_0, epsilon_0
from .cabos import (
    CoaxialCable,
    PipeCable,
    count_conductors_cable,
)

# Constantes
TOL = 1e-6

# %% Calculation functions (individual formulas)


def calc_inner_skin_effect_impedance(
    radius_in: float,
    radius_ext: float,
    rho_c: float,
    mur_c: float,
    complex_frequency: complex,
) -> complex:
    """
    Impedance of a tubular conductor considering the skin effect in its inner surface.

    Parameters
    ----------
    radius_in : float
        internal radius of the conductor [m]. Use zero if solid.
    radius_ext : float
        external radius of the conductor [m].
    rho_c : float
        resistivity of the conductor [Ω/m].
    mur_c : float
        relative magnetic permeability of the conductor [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    simplified_formula : bool, optional
        use a simplified formula? Default is False.

    Returns
    -------
    complex
        The calculated impedance.
    """
    if abs(radius_in) < TOL:
        return 0.0 + 0.0j

    # Constants
    mu_c = mu_0 * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = np.sqrt(complex_frequency * mu_c * sigma_c)

    w_out = m * radius_ext
    w_in = m * radius_in
    s_in = np.exp(abs(np.real(w_in)) - w_out)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", category=RuntimeWarning)
            s_out = np.exp(abs(np.real(w_out)) - w_in)
            sc = s_in / s_out  # Should be applied to all ive() involving w_in
    except RuntimeWarning:  # s_out = Inf
        sc = 0.0

    N = sc * ive(0, w_in) * kve(1, w_out) + kve(0, w_in) * ive(1, w_out)
    D = ive(1, w_out) * kve(1, w_in) - sc * kve(1, w_out) * ive(1, w_in)
    zin = (complex_frequency * mu_c / (2 * pi)) * (1 / w_in) * (N / D)
    return zin


def calc_outer_skin_effect_impedance(
    radius_in: float,
    radius_ext: float,
    rho_c: float,
    mur_c: float,
    complex_frequency: complex,
    simplified_formula: bool = False,
) -> complex:
    """
    Impedance of a tubular conductor considering the skin effect in its outer surface.

    Parameters
    ----------
    radius_in : float
        internal radius of the conductor [m]. Use zero if solid.
    radius_ext : float
        external radius of the conductor [m].
    rho_c : float
        resistivity of the conductor [Ω/m].
    mur_c : float
        relative magnetic permeability of the conductor [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    simplified_formula : bool, optional
        use a simplified formula? Default is False.

    Returns
    -------
    complex
        The calculated impedance.
    """
    # Constants
    mu_c = mu_0 * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = np.sqrt(complex_frequency * mu_c * sigma_c)

    w_out = m * radius_ext
    w_in = m * radius_in

    if radius_in < TOL:
        N = ive(0, w_out)
        D = ive(1, w_out)
    else:
        s_in = np.exp(abs(np.real(w_in)) - w_out)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", category=RuntimeWarning)
                s_out = np.exp(abs(np.real(w_out)) - w_in)
                sc = (
                    s_in / s_out
                )  # Should be applied to all ive() involving w_in
        except RuntimeWarning:  # s_out = Inf
            sc = 0.0
        N = ive(0, w_out) * kve(1, w_in) + sc * kve(0, w_out) * ive(1, w_in)
        D = ive(1, w_out) * kve(1, w_in) - sc * kve(1, w_out) * ive(1, w_in)

    zin = (complex_frequency * mu_c / (2 * pi)) * (1 / w_out) * (N / D)
    return zin


def calc_mutual_skin_effect_impedance(
    radius_in: float,
    radius_ext: float,
    rho_c: float,
    mur_c: float,
    complex_frequency: complex,
    simplified_formula: bool = False,
) -> complex:
    """
    Mutual impedance between the inner and outer surfaces of tubular conductor considering the skin effect.

    Parameters
    ----------
    radius_in : float
        internal radius of the conductor [m]. Use zero if solid.
    radius_ext : float
        external radius of the conductor [m].
    rho_c : float
        resistivity of the conductor [Ω/m].
    mur_c : float
        relative magnetic permeability of the conductor [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    simplified_formula : bool, optional
        use a simplified formula? Default is False.

    Returns
    -------
    complex
        The calculated mutual impedance.
    """
    if abs(radius_in) < TOL:
        return calc_outer_skin_effect_impedance(
            radius_in, radius_ext, rho_c, mur_c, complex_frequency
        )

    # Constants
    mu_c = mu_0 * mur_c
    sigma_c = 1.0 / rho_c

    # Calculate the reciprocal of the skin depth
    m = np.sqrt(complex_frequency * mu_c * sigma_c)

    w_out = m * radius_ext
    w_in = m * radius_in
    s_in = np.exp(abs(np.real(w_in)) - w_out)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", category=RuntimeWarning)
            s_out = np.exp(abs(np.real(w_out)) - w_in)
            sc = s_in / s_out  # Should be applied to all ive() involving w_in
            D = ive(1, w_out) * kve(1, w_in) - sc * kve(1, w_out) * ive(
                1, w_in
            )
            zm = 1 / (2 * pi * radius_ext * radius_in * sigma_c * D * s_out)
    except RuntimeWarning:  # s_out = Inf
        zm = 0.0
    return zm


def calc_tubular_capacitor_impedance(
    radius_in: float,
    radius_ext: float,
    mur_ins: float,
    complex_frequency: complex,
) -> complex:
    """
    Mutual impedance between the inner and outer surfaces of a tubular capacitor considering the skin effect.

    Parameters
    ----------
    radius_in : float
        internal radius of the capacitor [m]. Use zero if solid.
    radius_ext : float
        external radius of the capacitor [m].
    mur_ins : float
        relative magnetic permeability of the capacitor [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].

    Returns
    -------
    complex
        The calculated impedance.
    """
    return (
        complex_frequency
        * mu_0
        * mur_ins
        * np.log(radius_ext / radius_in)
        / (2 * pi)
    )


def calc_pipe_mutual_internal_impedance(
    distance_1: float,
    distance_2: float,
    theta: float,
    radius_1: float,
    radius_2: float,
    radius_in: float,
    rho_c: float,
    mur_c: float,
    complex_frequency: complex,
    maxiter: int = 1000,
    tol: float = 1e-12,
) -> complex:
    """
    Mutual impedance between two cables inside a tubular conductive pipe in relation to the pipe's internal surface.

    Parameters
    ----------
    distance_1 : float
        distance from the center of the pipe to the center of the first conductor [m].
    distance_2 : float
        distance from the center of the pipe to the center of the second conductor [m].
    theta : float
        angle between distance_1 and distance_2 [rad].
    radius_1 : float
        external radius of the insulator of the first conductor [m].
    radius_2 : float
        external radius of the insulator of the second conductor [m].
    radius_in : float
        internal radius of the pipe [m].
    rho_c : float
        resistivity of the pipe [Ω/m].
    mur_c : float
        relative magnetic permeability of the pipe [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    maxiter : int, optional
        maximum number of terms to consider in the infinite series. Default is 1000.
    tol : float, optional
        tolerance to check for convergence of the infinite series. Default is 1e-12.

    Returns
    -------
    complex
        The calculated mutual impedance.
    """
    sigma_c = 1.0 / rho_c

    r1 = radius_in
    di = distance_1
    dk = distance_2

    r_sq = r1**2 + 0.0j
    didk = di * dk

    w = r1 * np.sqrt(complex_frequency * mu_0 * mur_c * sigma_c)

    K0 = kve(0, w)
    K1 = kve(1, w)
    # inner impedance of a pipe of infinite thickness
    inner_term = mur_c * K0 / (w * K1)

    # Qik is the impedance through the insulator
    if abs(theta) < TOL and abs(di - dk) < TOL:  # i == k
        Qik = np.log(r1 / radius_1 * (1 - di**2 / r_sq))
    else:  # i ≠ k
        ln_1 = np.log((r_sq - didk * np.exp(1j * theta)) / r_sq)
        ln_2 = np.log((r_sq - didk * np.exp(-1j * theta)) / r_sq)
        Cn_sum = -0.5 * (ln_1 + ln_2)
        Qik = (
            np.log(
                radius_in / np.sqrt(di**2 + dk**2 - 2 * didk * np.cos(theta))
            )
            - Cn_sum
        )

    if abs(di) < TOL or abs(dk) < TOL:  # results in `Cn = 0` ∀n
        return complex_frequency * mu_0 / (2 * pi) * (inner_term + Qik)

    # mutual impedance between both inner cables
    def Cn(n):
        return (didk / r_sq) ** n * np.cos(n * theta)

    def Dn(n, K_prev, K_n):
        return n * (1 + mur_c) + w * K_prev / K_n

    K_prev = K0
    K_n = K1
    inf_series = Cn(1) / Dn(1, K_prev, K_n)
    converged = False
    overflowing = False

    for n in range(2, maxiter + 1):
        if not overflowing:
            K_prev = K_n
            K_n = kve(n, w)
            term_n = Cn(n) / Dn(n, K_prev, K_n)
            if np.isnan(K_n):
                overflowing = True

        if overflowing:
            term_n = Cn(n) / (n * (1 + mur_c))

        inf_series += term_n
        if n > 5 and abs(term_n) < tol * abs(inf_series):
            converged = True
            break
        elif n == maxiter and not converged:
            warnings.warn(
                f"calc_pipe_mutual_internal_impedance did not converge to specified tolerance. Last term: {abs(term_n) / abs(inf_series)}"
            )

    return (
        complex_frequency
        * mu_0
        / (2 * pi)
        * (inner_term + Qik + 2 * mur_c * inf_series)
    )


def calc_pipe_self_internal_impedance(
    distance_1: float,
    radius_1: float,
    radius_in: float,
    rho_c: float,
    mur_c: float,
    complex_frequency: complex,
    maxiter: int = 1000,
    tol: float = 1e-12,
) -> complex:
    """
    Self impedance of a cable inside a tubular conductive pipe in relation to the pipe's internal surface.

    Parameters
    ----------
    distance_1 : float
        distance from the center of the pipe to the center of the first conductor [m].
    radius_1 : float
        external radius of the insulator of the first conductor [m].
    radius_in : float
        internal radius of the pipe [m].
    rho_c : float
        resistivity of the pipe [Ω/m].
    mur_c : float
        relative magnetic permeability of the pipe [dimensionless].
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    maxiter : int, optional
        maximum number of terms to consider in the infinite series. Default is 1000.
    tol : float, optional
        tolerance to check for convergence of the infinite series. Default is 1e-12.

    Returns
    -------
    complex
        The calculated self impedance.
    """
    if abs(distance_1) < TOL and abs(radius_1 - radius_in) < TOL:
        return 0.0 + 0.0j

    return calc_pipe_mutual_internal_impedance(
        distance_1,
        distance_1,
        0.0,
        radius_1,
        radius_1,
        radius_in,
        rho_c,
        mur_c,
        complex_frequency,
        maxiter=maxiter,
        tol=tol,
    )


# %% Computation functions


def comp_coaxial_cable_impedance(cable, complex_frequency: complex):
    """
    Compute the series impedance matrix of a coaxial cable at the given complex frequency in rad/s.

    Parameters
    ----------
    cable : CoaxialCable
        The coaxial cable object.
    complex_frequency : complex
        Complex angular frequency s = c + jω [rad/s].

    Returns
    -------
    numpy.ndarray
        Impedance matrix.
    """
    # number of conductors
    Nc = len(cable.components)
    Z = np.zeros((Nc, Nc), dtype=complex)

    for i in range(Nc):
        comp = cable.components[i]
        rho_c = comp.rho_c
        mur_c = comp.mur_c
        radius_in = comp.radius_in
        radius_ext = comp.radius_ext
        radius_ext_insulator = comp.radius_ext_insulator
        mur_d = comp.mur_c

        z_inner = calc_inner_skin_effect_impedance(
            radius_in, radius_ext, rho_c, mur_c, complex_frequency
        )

        z_outer = calc_outer_skin_effect_impedance(
            radius_in, radius_ext, rho_c, mur_c, complex_frequency
        )

        if i > 0:
            z_mutual_io = calc_mutual_skin_effect_impedance(
                radius_in, radius_ext, rho_c, mur_c, complex_frequency
            )
        else:
            z_mutual_io = 0.0 + 0.0j

        z_mutual_d = calc_tubular_capacitor_impedance(
            radius_ext, radius_ext_insulator, mur_d, complex_frequency
        )

        if i > 0:
            Z[0:i, 0:i] += z_outer + z_mutual_d - 2 * z_mutual_io + z_inner
            Z[i, 0:i] += z_outer + z_mutual_d - z_mutual_io
            Z[0:i, i] += z_outer + z_mutual_d - z_mutual_io

        Z[i, i] += z_outer + z_mutual_d

    return Z


def comp_pipe_cable_impedance(pipecable, complex_frequency: complex):
    """
    Compute the series impedance matrix of a pipe-type cable at the given complex frequency in rad/s.

    Parameters
    ----------
    pipecable : PipeCable
        The pipe cable object.
    complex_frequency : complex
        Complex angular frequency s = c + jω [rad/s].

    Returns
    -------
    numpy.ndarray
        Impedance matrix.
    """
    maxiter = 10000
    tol = 1e-12

    Nc = count_conductors_cable(pipecable)
    Z = np.zeros((Nc, Nc), dtype=complex)

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

    z_inner = calc_inner_skin_effect_impedance(
        pipecable.radius_in,
        pipecable.radius_ext,
        pipecable.rho_c,
        pipecable.mur_c,
        complex_frequency,
    )

    Z[0 : Nc - 1, 0 : Nc - 1] = (
        z_outer + z_mutual_d - 2 * z_mutual_io + z_inner
    )
    Z[0 : Nc - 1, Nc - 1] = z_outer + z_mutual_d - z_mutual_io
    Z[Nc - 1, 0 : Nc - 1] = z_outer + z_mutual_d - z_mutual_io
    Z[Nc - 1, Nc - 1] = z_outer + z_mutual_d

    N_inner = len(pipecable.cables)
    next_index_i = 0

    for i in range(N_inner):
        cable_i = pipecable.cables[i]
        x1 = cable_i.x - pipecable.x
        y1 = cable_i.y - pipecable.y
        distance_1 = np.hypot(x1, y1)
        radius_1 = cable_i.outer_radius()

        zm_ii = calc_pipe_self_internal_impedance(
            distance_1,
            radius_1,
            pipecable.radius_in,
            pipecable.rho_c,
            pipecable.mur_c,
            complex_frequency,
            maxiter=maxiter,
            tol=tol,
        )

        i1 = next_index_i
        i2 = i1 + count_conductors_cable(cable_i) - 1
        Z[i1 : i2 + 1, i1 : i2 + 1] += zm_ii
        next_index_i = i2 + 1

        next_index_k = next_index_i
        for k in range(i + 1, N_inner):
            cable_k = pipecable.cables[k]
            x2 = cable_k.x - pipecable.x
            y2 = cable_k.y - pipecable.y
            distance_2 = np.hypot(x2, y2)
            radius_2 = cable_k.outer_radius()

            # Calculate angle between vectors
            dot_product = x1 * x2 + y1 * y2
            magnitude_product = distance_1 * distance_2
            if abs(magnitude_product) < TOL:
                theta = 0.0
            else:
                theta = np.arccos(dot_product / magnitude_product)

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
                maxiter=maxiter,
                tol=tol,
            )

            k1 = next_index_k
            k2 = k1 + count_conductors_cable(cable_k) - 1
            Z[i1 : i2 + 1, k1 : k2 + 1] += zm_ik
            Z[k1 : k2 + 1, i1 : i2 + 1] = Z[i1 : i2 + 1, k1 : k2 + 1]
            next_index_k = k2 + 1

    return Z


def comp_cable_impedance_recursive(Z, cable, complex_frequency, start_index):
    """
    Recursive function to compute impedance for any cable type.

    A PipeCable and the components of a CoaxialCable receive an `_index`
    attribute mapping the conductor to the matrix position.

    Parameters
    ----------
    Z : numpy.ndarray
        impedance matrix of the cable system that is modified in-place.
    cable : AbstractCable
        current cable being iterated upon.
    complex_frequency : complex
        complex angular frequency s = c + jω [rad/s].
    start_index : int
        current index in the impedance matrix.

    Returns
    -------
    int
        index in the impedance matrix after the computation.
    """
    n = count_conductors_cable(cable)
    end_index = start_index + n
    Z_k = Z[start_index:end_index, start_index:end_index]

    if isinstance(cable, CoaxialCable):
        Z_k += comp_coaxial_cable_impedance(cable, complex_frequency)
        for i, comp in enumerate(cable.components, start_index):
            comp._index = i

    elif isinstance(cable, PipeCable):
        Z_k += comp_pipe_cable_impedance(cable, complex_frequency)
        current_index = start_index
        cable._index = end_index - 1
        for i, inner_cable in enumerate(cable.cables, start_index):
            current_index = comp_cable_impedance_recursive(
                Z, inner_cable, complex_frequency, current_index
            )

    else:
        raise ValueError(f"Unrecognized cable type `{type(cable)}`.")

    return end_index


def comp_cable_system_impedance(
    cable_system, complex_frequencies, sigma_mar=5.0, epsr_mar=81
):
    """
    Compute the series impedance matrix of a cable system at the given frequencies of interest in rad/s.

    A PipeCable and the components of a CoaxialCable receive an `_index`
    attribute mapping the conductor to the matrix position.

    Parameters
    ----------
    cable_system : PipeCable
        The outermost cable.
    complex_frequencies : list of complex
        List of complex angular frequencies s = c + jω [rad/s].
    sigma_mar: Opcional, float
        Permissividade elétrica do mar [S/m]. O padrão é `5.0`.
        Use None para ignorar o mar.
    epsr_mar: Opcional, float
        Permissividade elétrica do mar. O padrão é `81.0`.

    Returns
    -------
    numpy.ndarray
        Impedance matrix of shape (Nc, Nc, Nf).
    """
    Nc = count_conductors_cable(cable_system)
    Z = np.zeros((Nc, Nc, len(complex_frequencies)), dtype=complex)
    rca = cable_system.radius_ext_insulator
    for k, jω in enumerate(complex_frequencies):
        Z_k = Z[:, :, k]
        current_index = 0
        current_index = comp_cable_impedance_recursive(
            Z_k, cable_system, jω, current_index
        )
        if sigma_mar is not None:
            Z[:, :, k] += cZmar(jω, rca, sigma_mar, epsr_mar)
    return Z


def cZmar(
    freq_s: complex,
    rca: float,
    sigma: Optional[float] = 5.0,
    epsr: Optional[float] = 81.0,
) -> complex:
    """Impedância de retorno pelo mar do cabo pipe-Type (PT).

    Parâmetros
    ----------
    freq_s: complex
        Frequência angular complexa no formato `c + jw` [rad/s].
    rca: float
        Raio da capa externa da armadura [m].
    sigma: Opcional, float
        Permissividade elétrica do mar [S/m]. O padrão é `5.0`.
        Use None para ignorar o mar (retorna 0).
    epsr: Opcional, float
        Permissividade elétrica do mar. O padrão é `81.0`.

    Retorna
    -------
        Z: impedância do mar [Ω/m].
    """
    if sigma is None or epsr is None:
        return 0.0

    jw = freq_s
    rho = 1.0 / sigma
    eta = np.sqrt(jw * mu_0 * (sigma + jw * epsilon_0 * epsr))
    Zmar = eta * rho / (2 * pi * rca) * kve(0, eta * rca) / kve(1, eta * rca)
    return Zmar
