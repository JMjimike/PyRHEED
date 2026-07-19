"""Probe-specific geometry wrappers for RHEED and XRD simulations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

import ewald
import xrd_ewald
from four_circle_geometry import (
    FourCircleAngles,
    bragg_error,
    forward_kinematics,
    g_vector_in_lab,
    ki_in_crystal_frame,
    recommend_scan_mode,
    sample_rotation_matrix,
    solve_bragg_angles,
)


class ProbeGeometry(Protocol):
    def k0(self) -> float: ...

    def incident_k(self) -> np.ndarray: ...

    def ewald_weight(self, Kx, Ky, Kz) -> np.ndarray: ...

    def intensity_corrections(self, Kx, Ky, Kz) -> np.ndarray: ...


@dataclass
class RHEEDGeometry:
    energy_kev: float
    fwhm_ev: float = 0.0
    theta_deg: float = 2.0

    def k0(self):
        return ewald.electron_wavenumber(self.energy_kev)

    def incident_k(self):
        return ewald.incident_wavevector_y(self.k0(), self.theta_deg)

    def ewald_weight(self, Kx, Ky, Kz):
        ki = self.incident_k()
        k0 = self.k0()
        delta = ewald.ewald_deviation(Kx, Ky, Kz, ki, k0)
        sigma = ewald.gaussian_energy_sigma_delta(k0, self.energy_kev, self.fwhm_ev)
        if sigma <= 0.0:
            return (np.abs(delta) < 1e-12).astype(float)
        return np.exp(-0.5 * (delta / sigma) ** 2)

    def intensity_corrections(self, Kx, Ky, Kz):
        return np.ones_like(Kx, dtype=float)

    def apply_broadening(self, intensity, Kx, Ky, Kz):
        return ewald.apply_ewald_broadening(
            intensity,
            Kx,
            Ky,
            Kz,
            self.energy_kev,
            self.fwhm_ev,
            self.theta_deg,
        )


@dataclass
class XRDGeometry:
    wavelength_angstrom: float = 1.54184
    fwhm_deg: float = 0.0
    omega_deg: float = 0.0
    chi_deg: float = 0.0
    phi_deg: float = 0.0
    two_theta_deg: float = 0.0
    r_mount: np.ndarray | None = None
    apply_lorentz: bool = True

    def k0(self):
        return xrd_ewald.xray_wavenumber(self.wavelength_angstrom)

    def incident_k_lab(self):
        ki, _, _ = forward_kinematics(
            self.omega_deg,
            self.chi_deg,
            self.phi_deg,
            self.two_theta_deg,
            self.wavelength_angstrom,
            self.r_mount,
        )
        return ki

    def incident_k(self):
        return ki_in_crystal_frame(
            self.omega_deg,
            self.chi_deg,
            self.phi_deg,
            self.two_theta_deg,
            self.wavelength_angstrom,
            self.r_mount,
        )

    def ewald_weight(self, Kx, Ky, Kz):
        ki = self.incident_k()
        k0 = self.k0()
        delta = xrd_ewald.ewald_deviation(Kx, Ky, Kz, ki, k0)
        sigma = xrd_ewald.angular_sigma_delta(k0, self.fwhm_deg)
        if sigma <= 0.0:
            return (np.abs(delta) < 1e-9).astype(float)
        return np.exp(-0.5 * (delta / sigma) ** 2)

    def intensity_corrections(self, Kx, Ky, Kz):
        if not self.apply_lorentz:
            return np.ones_like(Kx, dtype=float)
        two_theta = np.full_like(Kx, self.two_theta_deg, dtype=float)
        vectorized = np.vectorize(xrd_ewald.lorentz_polarization_factor)
        return vectorized(two_theta)

    def apply_broadening(self, intensity, Kx, Ky, Kz):
        weighted = np.real(intensity).astype(float) * self.ewald_weight(Kx, Ky, Kz)
        return weighted * self.intensity_corrections(Kx, Ky, Kz)


@dataclass
class FourCircleGeometry(XRDGeometry):
    """Four-circle XRD geometry with explicit sample/detector angles."""

    def current_angles(self) -> FourCircleAngles:
        return FourCircleAngles(
            self.omega_deg, self.chi_deg, self.phi_deg, self.two_theta_deg
        )

    def forward_kin(self):
        return forward_kinematics(
            self.omega_deg,
            self.chi_deg,
            self.phi_deg,
            self.two_theta_deg,
            self.wavelength_angstrom,
            self.r_mount,
        )

    def solve_bragg(self, hkl, reciprocal_matrix):
        return solve_bragg_angles(
            hkl,
            reciprocal_matrix,
            self.wavelength_angstrom,
            r_mount=self.r_mount,
            current_angles=self.current_angles(),
        )

    def bragg_error(self, hkl, reciprocal_matrix):
        from four_circle_geometry import crystal_g_vector

        g_crystal = crystal_g_vector(hkl, reciprocal_matrix)
        return bragg_error(
            self.omega_deg,
            self.chi_deg,
            self.phi_deg,
            self.two_theta_deg,
            g_crystal,
            self.wavelength_angstrom,
            self.r_mount,
        )


def auto_k_range(wavelength_angstrom, two_theta_max_deg=80.0, margin=1.05):
    """Suggest symmetric K grid limits for XRD."""
    q_max = xrd_ewald.max_q_magnitude(wavelength_angstrom, two_theta_max_deg)
    limit = q_max * margin
    return (-limit, limit)
