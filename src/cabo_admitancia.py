# -*- coding: utf-8 -*-
"""Transversal admittance Y per unit length in frequency-domain of Single-Core (coaxial) and Pipe-Type cables."""

import numpy as np
from scipy.linalg import inv
import math
from typing import List
from scipy.constants import pi, epsilon_0
from .cabos import (
    CoaxialCable,
    PipeCable,
    count_conductors_cable,
)


# %% Calculation functions (individual formulas)


def calc_tubular_elastance(
    radius_in: float, radius_ext: float, epsr: float
) -> complex:
    """
    Elastance of a tubular capacitor.

    Parameters
    ----------
    radius_in : float
        internal radius of the capacitor [m].
    radius_ext : float
        external radius of the capacitor [m].
    epsr : float
        relative electric permittivity of the insulator [dimensionless].

    Returns
    -------
    complex
        The calculated elastance.
    """
    return np.log(radius_ext / radius_in) / (2 * pi * epsilon_0 * epsr)


def calc_pipe_mutual_internal_elastance(
    distance_1: float,
    distance_2: float,
    theta: float,
    radius_1: float,
    radius_2: float,
    radius_in: float,
    epsr: float,
) -> complex:
    """
    Mutual elastance between two cables inside a tubular conductive pipe in relation to the pipe's internal surface.

    Parameters
    ----------
    distance_1 : float
        distance from the center of the pipe to the center of conductor 1 [m].
    distance_2 : float
        distance from the center of the pipe to the center of conductor 2 [m].
    theta : float
        angle between distance_1 and distance_2 [rad].
    radius_1 : float
        external radius of conductor 1 [m].
    radius_2 : float
        external radius of conductor 2 [m].
    radius_in : float
        internal radius of the pipe [m].
    epsr : float
        relative electric permittivity of the pipe's internal insulator [dimensionless].

    Returns
    -------
    complex
        The calculated mutual elastance.
    """
    r1 = radius_in
    di = distance_1
    dk = distance_2

    r_sq = r1**2 + 0.0j
    didk = di * dk

    TOL = 1e-12

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

    return Qik / (2 * pi * epsilon_0 * epsr)


def calc_pipe_self_internal_elastance(
    distance_1: float, radius_1: float, radius_in: float, epsr: float
) -> complex:
    """
    Self elastance of a cable inside a tubular conductive pipe in relation to the pipe's internal surface.

    Parameters
    ----------
    distance_1 : float
        distance from the center of the pipe to the center of conductor 1 [m].
    radius_1 : float
        external radius of conductor 1 [m].
    radius_in : float
        internal radius of the pipe [m].
    epsr : float
        relative electric permittivity of the pipe's internal insulator [dimensionless].

    Returns
    -------
    complex
        The calculated self elastance.
    """
    return calc_pipe_mutual_internal_elastance(
        distance_1, distance_1, 0.0, radius_1, radius_1, radius_in, epsr
    )


# %% Computation functions


def comp_coaxial_cable_elastance(cable) -> np.ndarray:
    """
    Compute the shunt elastance matrix of a coaxial cable.

    Parameters
    ----------
    cable : CoaxialCable
        The coaxial cable object.

    Returns
    -------
    numpy.ndarray
        Elastance matrix.
    """
    # number of conductors
    Nc = len(cable.components)
    P = np.zeros((Nc, Nc), dtype=complex)

    for i in range(Nc):
        comp = cable.components[i]
        r1 = comp.radius_ext
        r2 = comp.radius_ext_insulator
        epsr = comp.epsr
        P[0 : i + 1, 0 : i + 1] += calc_tubular_elastance(r1, r2, epsr)

    return P


def comp_pipe_cable_elastance(pipecable) -> np.ndarray:
    """
    Compute the shunt elastance matrix of a pipe-type cable.

    Parameters
    ----------
    pipecable : PipeCable
        The pipe cable object.

    Returns
    -------
    numpy.ndarray
        Elastance matrix.
    """
    Nc = count_conductors_cable(pipecable)
    P = np.zeros((Nc, Nc), dtype=complex)

    p_c = calc_tubular_elastance(
        pipecable.radius_ext,
        pipecable.radius_ext_insulator,
        pipecable.epsr_ext,
    )

    P[:, :] = p_c

    N_inner = len(pipecable.cables)
    next_index_i = 0

    for i in range(N_inner):
        cable_i = pipecable.cables[i]
        x1 = cable_i.x - pipecable.x
        y1 = cable_i.y - pipecable.y
        distance_1 = math.hypot(x1, y1)
        radius_1 = cable_i.outer_radius()

        p_ii = calc_pipe_self_internal_elastance(
            distance_1, radius_1, pipecable.radius_in, pipecable.epsr_in
        )

        i1 = next_index_i
        i2 = i1 + count_conductors_cable(cable_i) - 1
        P[i1 : i2 + 1, i1 : i2 + 1] += p_ii
        next_index_i = i2 + 1

        next_index_k = next_index_i
        for k in range(i + 1, N_inner):
            cable_k = pipecable.cables[k]
            x2 = cable_k.x - pipecable.x
            y2 = cable_k.y - pipecable.y
            distance_2 = math.hypot(x2, y2)
            radius_2 = cable_k.outer_radius()

            # Calculate angle between vectors
            dot_product = x1 * x2 + y1 * y2
            magnitude_product = distance_1 * distance_2
            TOL = 1e-12
            if abs(magnitude_product) < TOL:
                theta = 0.0
            else:
                theta = math.acos(dot_product / magnitude_product)

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
            P[i1 : i2 + 1, k1 : k2 + 1] += p_ik
            P[k1 : k2 + 1, i1 : i2 + 1] = P[i1 : i2 + 1, k1 : k2 + 1]
            next_index_k = k2 + 1

    return P


def comp_cable_elastance_recursive(
    P: np.ndarray, cable, start_index: int
) -> int:
    """
    Recursive function to compute elastance for any cable type.

    A PipeCable and the components of a CoaxialCable receive an `_index`
    attribute mapping the conductor to the matrix position.

    Parameters
    ----------
    P : numpy.ndarray
        elastance matrix of the cable system that is modified in-place.
    cable : AbstractCable
        current cable being iterated upon.
    start_index : int
        current index in the elastance matrix.

    Returns
    -------
    int
        index in the elastance matrix after the computation.
    """
    n = count_conductors_cable(cable)
    end_index = start_index + n
    P_k = P[start_index:end_index, start_index:end_index]

    if isinstance(cable, CoaxialCable):
        P_k += comp_coaxial_cable_elastance(cable)
        for i, comp in enumerate(cable.components, start_index):
            comp._index = i

    elif isinstance(cable, PipeCable):
        P_k += comp_pipe_cable_elastance(cable)
        current_index = start_index
        cable._index = end_index - 1
        for inner_cable in cable.cables:
            current_index = comp_cable_elastance_recursive(
                P, inner_cable, current_index
            )
    else:
        raise ValueError(f"Unrecognized cable type `{type(cable)}`.")

    return end_index


def comp_cable_system_elastance(cable_system) -> np.ndarray:
    """
    Compute the shunt elastance matrix of a cable system.

    A PipeCable and the components of a CoaxialCable receive an `_index`
    attribute mapping the conductor to the matrix position.

    Parameters
    ----------
    cable_system : PipeCable
        The outermost cable.

    Returns
    -------
    numpy.ndarray
        Elastance matrix.
    """
    Nc = count_conductors_cable(cable_system)
    P = np.zeros((Nc, Nc), dtype=complex)
    current_index = 0
    comp_cable_elastance_recursive(P, cable_system, current_index)
    return P


def comp_cable_system_admittance(
    cable_system, complex_frequencies: List[complex]
) -> np.ndarray:
    """
    Compute the shunt admittance matrix of a cable system at the given frequencies of interest in rad/s.

    Parameters
    ----------
    cable_system : PipeCable
        The outermost cable.
    complex_frequencies : list of complex
        List of complex angular frequencies s = c + jω [rad/s].

    Returns
    -------
    numpy.ndarray
        Admittance matrix of shape (Nc, Nc, Nf).
    """
    P = comp_cable_system_elastance(cable_system)
    P_inv = inv(P)
    P_inv = (P_inv + P_inv.T) / 2.0  # ensure symmetry
    Nc = P.shape[0]
    Y = np.zeros((Nc, Nc, len(complex_frequencies)), dtype=complex)

    for k, jω in enumerate(complex_frequencies):
        Y[:, :, k] = jω * P_inv

    return Y
