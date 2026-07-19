"""Euler four-circle diffractometer kinematics (Bruker/Huber vertical-omega convention)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from xrd_ewald import incident_wavevector_z, scattered_wavevector_y, xray_wavenumber


@dataclass(frozen=True)
class FourCircleAngles:
    omega: float
    chi: float
    phi: float
    two_theta: float

    def as_tuple(self):
        return (self.omega, self.chi, self.phi, self.two_theta)

    def distance_to(self, other: FourCircleAngles) -> float:
        diffs = [
            _angle_diff(self.omega, other.omega),
            _angle_diff(self.chi, other.chi),
            _angle_diff(self.phi, other.phi),
            _angle_diff(self.two_theta, other.two_theta),
        ]
        return float(np.linalg.norm(diffs))


def _angle_diff(a_deg, b_deg):
    return ((a_deg - b_deg + 180.0) % 360.0) - 180.0


def rotation_x(angle_deg):
    angle = np.deg2rad(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def rotation_y(angle_deg):
    angle = np.deg2rad(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def rotation_z(angle_deg):
    angle = np.deg2rad(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def sample_rotation_matrix(omega_deg, chi_deg, phi_deg, r_mount=None):
    """R_sample = R_omega · R_chi · R_phi · R_mount."""
    r_mount = np.eye(3) if r_mount is None else np.asarray(r_mount, dtype=float)
    return rotation_y(omega_deg) @ rotation_x(chi_deg) @ rotation_z(phi_deg) @ r_mount


def forward_kinematics(
    omega_deg,
    chi_deg,
    phi_deg,
    two_theta_deg,
    wavelength_angstrom,
    r_mount=None,
):
    """
    Return (ki_lab, ks_lab, Q_lab) for specular theta-2theta in the XZ plane.

    Beams are fixed in the laboratory frame; sample rotation (omega, chi, phi)
    enters through ki_in_crystal_frame() for Ewald weighting on the crystal grid.
    """
    del omega_deg, chi_deg, phi_deg, r_mount
    theta_b = np.deg2rad(two_theta_deg) / 2.0
    k0 = xray_wavenumber(wavelength_angstrom)
    ki = k0 * np.array([np.cos(theta_b), 0.0, -np.sin(theta_b)])
    ks = k0 * np.array([np.cos(theta_b), 0.0, np.sin(theta_b)])
    q_lab = ks - ki
    return ki, ks, q_lab


def scattering_vector_lab(
    omega_deg,
    chi_deg,
    phi_deg,
    two_theta_deg,
    wavelength_angstrom,
    r_mount=None,
):
    """Laboratory-frame scattering vector Q = k_s - k_i."""
    _ki, _ks, q_lab = forward_kinematics(
        omega_deg, chi_deg, phi_deg, two_theta_deg, wavelength_angstrom, r_mount
    )
    return q_lab


def ki_in_crystal_frame(
    omega_deg,
    chi_deg,
    phi_deg,
    two_theta_deg,
    wavelength_angstrom,
    r_mount=None,
):
    """Transform the fixed laboratory incident beam into crystal coordinates."""
    ki_lab, _, _ = forward_kinematics(
        omega_deg, chi_deg, phi_deg, two_theta_deg, wavelength_angstrom, r_mount
    )
    r_sample = sample_rotation_matrix(omega_deg, chi_deg, phi_deg, r_mount)
    return r_sample.T @ ki_lab


def q_lab_to_film_coords(q_lab, r_mount=None):
    """
    Project Q_lab onto film in-plane (X) and out-of-plane (Z) axes.

    Assumes the oriented film has surface normal along laboratory +Z and
    in-plane reference along +X after mount rotation.
    """
    q_lab = np.asarray(q_lab, dtype=float).reshape(3)
    return float(q_lab[0]), float(q_lab[2])


def crystal_g_vector(hkl, reciprocal_matrix):
    """Reciprocal-space vector for Miller indices in Cartesian coordinates."""
    hkl_vec = np.asarray(hkl, dtype=float)
    return hkl_vec @ np.asarray(reciprocal_matrix, dtype=float)


def g_vector_in_lab(hkl, reciprocal_matrix, omega_deg, chi_deg, phi_deg, r_mount=None):
    """Transform crystal g(hkl) into the laboratory frame."""
    g_crystal = crystal_g_vector(hkl, reciprocal_matrix)
    r_sample = sample_rotation_matrix(omega_deg, chi_deg, phi_deg, r_mount)
    return r_sample @ g_crystal


def bragg_error(
    omega_deg,
    chi_deg,
    phi_deg,
    two_theta_deg,
    g_crystal,
    wavelength_angstrom,
    r_mount=None,
):
    """Combined Bragg and direction mismatch for inverse kinematics."""
    g_crystal = np.asarray(g_crystal, dtype=float)
    g_norm = np.linalg.norm(g_crystal)
    if g_norm < 1e-12:
        return float("inf")

    k0 = xray_wavenumber(wavelength_angstrom)
    ki, ks, q_target = forward_kinematics(
        omega_deg, chi_deg, phi_deg, two_theta_deg, wavelength_angstrom, r_mount
    )
    g_lab = sample_rotation_matrix(omega_deg, chi_deg, phi_deg, r_mount) @ g_crystal
    closure_err = float(np.linalg.norm((ki + g_lab) - ks))
    q_norm = np.linalg.norm(q_target)
    if q_norm < 1e-12:
        return closure_err + abs(np.linalg.norm(ki + g_lab) - k0)

    sphere_err = abs(np.linalg.norm(ki + g_lab) - k0)
    direction_err = np.linalg.norm(g_lab / g_norm - q_target / q_norm)
    magnitude_err = abs(g_norm - q_norm) / g_norm
    return closure_err + sphere_err + direction_err + magnitude_err


def is_symmetric_reflection(g_lab, z_hat=None, tolerance=0.99):
    z_hat = np.array([0.0, 0.0, 1.0]) if z_hat is None else np.asarray(z_hat, dtype=float)
    g_norm = np.linalg.norm(g_lab)
    if g_norm < 1e-12:
        return False
    return abs(np.dot(g_lab / g_norm, z_hat)) >= tolerance


def recommend_scan_mode(g_lab):
    if is_symmetric_reflection(g_lab):
        return "theta_2theta"
    g_xy = np.linalg.norm(g_lab[:2])
    g_norm = np.linalg.norm(g_lab)
    if g_norm > 0 and g_xy / g_norm > 0.99:
        return "phi_scan"
    return "four_circle"


def solve_bragg_angles(
    hkl,
    reciprocal_matrix,
    wavelength_angstrom,
    r_mount=None,
    current_angles: FourCircleAngles | None = None,
    two_theta_range=(5.0, 150.0),
    omega_range=(-30.0, 30.0),
    chi_range=(0.0, 90.0),
    phi_range=(0.0, 360.0),
    coarse_step=5.0,
    full_search=False,
):
    """
    Find four-circle settings that satisfy the Bragg condition for (hkl).

    Returns sorted list of FourCircleAngles (lowest error first).
    When full_search is False, symmetric reflections return immediately and
    non-symmetric reflections skip the expensive grid scan.
    """
    g_crystal = crystal_g_vector(hkl, reciprocal_matrix)
    g_norm = np.linalg.norm(g_crystal)
    if g_norm < 1e-12:
        return []

    wavelength = float(wavelength_angstrom)
    d_spacing = 2.0 * np.pi / g_norm
    if wavelength > 2.0 * d_spacing:
        return []

    theta_b = math.degrees(math.asin(wavelength / (2.0 * d_spacing)))
    two_theta_bragg = 2.0 * theta_b

    g_lab_zero = sample_rotation_matrix(0.0, 0.0, 0.0, r_mount) @ g_crystal
    solutions: list[tuple[float, FourCircleAngles]] = []
    if is_symmetric_reflection(g_lab_zero):
        candidate = FourCircleAngles(0.0, 0.0, 0.0, two_theta_bragg)
        err = bragg_error(0.0, 0.0, 0.0, two_theta_bragg, g_crystal, wavelength, r_mount)
        solutions.append((err, candidate))
        if err < 0.05 and not full_search:
            return [candidate]

    if not full_search:
        solutions.sort(key=lambda item: item[0])
        return [angles for _err, angles in solutions]

    omega_vals = _frange(*omega_range, coarse_step)
    chi_vals = _frange(*chi_range, coarse_step)
    phi_vals = _frange(*phi_range, coarse_step)
    two_theta_vals = _frange(
        max(two_theta_range[0], two_theta_bragg - 10.0),
        min(two_theta_range[1], two_theta_bragg + 10.0),
        max(coarse_step / 2.0, 1.0),
    )

    for omega in omega_vals:
        for chi in chi_vals:
            for phi in phi_vals:
                for two_theta in two_theta_vals:
                    err = bragg_error(
                        omega, chi, phi, two_theta, g_crystal, wavelength, r_mount
                    )
                    if err < 0.05:
                        solutions.append(
                            (err, FourCircleAngles(omega, chi, phi, two_theta))
                        )

    if not solutions:
        return []

    solutions.sort(key=lambda item: item[0])
    unique: list[FourCircleAngles] = []
    for _err, angles in solutions:
        if all(existing.distance_to(angles) > 1.0 for existing in unique):
            unique.append(angles)
    return unique


def pick_recommended_solution(
    solutions: Iterable[FourCircleAngles],
    current_angles: FourCircleAngles | None = None,
) -> FourCircleAngles | None:
    solutions = list(solutions)
    if not solutions:
        return None
    if current_angles is None:
        return solutions[0]
    return min(solutions, key=lambda angles: angles.distance_to(current_angles))


def hkl_tracker_angles(
    hkl,
    reciprocal_matrix,
    wavelength_angstrom,
    r_mount=None,
    current_angles: FourCircleAngles | None = None,
):
    """Return the best four-circle setting to track a single reflection."""
    solutions = solve_bragg_angles(
        hkl,
        reciprocal_matrix,
        wavelength_angstrom,
        r_mount=r_mount,
        current_angles=current_angles,
    )
    return pick_recommended_solution(solutions, current_angles)


def _frange(start, stop, step):
    values = []
    value = start
    while value <= stop + 1e-9:
        values.append(float(value))
        value += step
    return values
