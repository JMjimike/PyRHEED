"""Automatic reciprocal-space grid for XRD scans (Ewald-aware bounding box)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from four_circle_geometry import (
    FourCircleAngles,
    crystal_g_vector,
    sample_rotation_matrix,
    scattering_vector_lab,
)
from xrd_constraints import ScanConfig, ValidationResult, validate_scan_config
from xrd_ewald import angular_sigma_delta, max_q_magnitude, resolve_wavelength, xray_wavenumber
from xrd_scan import generate_trajectory


@dataclass
class KGridSpec:
    kx_min: float
    kx_max: float
    ky_min: float
    ky_max: float
    kz_min: float
    kz_max: float
    n_para: int
    n_perp: int
    summary: str = ""


def q_crystal_at_angles(angles: FourCircleAngles, wavelength, r_mount):
    """Scattering vector Q in crystal Cartesian coordinates."""
    q_lab = scattering_vector_lab(
        angles.omega,
        angles.chi,
        angles.phi,
        angles.two_theta,
        wavelength,
        r_mount,
    )
    r_sample = sample_rotation_matrix(
        angles.omega, angles.chi, angles.phi, r_mount
    )
    return r_sample.T @ q_lab


def _pad_from_fwhm(wavelength, fwhm_deg):
    k0 = xray_wavenumber(wavelength)
    sigma_k = angular_sigma_delta(k0, max(fwhm_deg, 0.01))
    return max(3.0 * sigma_k, 0.08)


def _steps_for_span(span, dk_target, min_steps=8, max_steps=64):
    if span <= 0:
        return min_steps
    return int(max(min_steps, min(max_steps, math.ceil(span / dk_target) + 1)))


def auto_k_grid_for_scan(
    config: ScanConfig,
    validation: ValidationResult | None = None,
) -> KGridSpec:
    """
    Build a K-grid box covering the Ewald trajectory plus angular broadening.

    Uses scan angles, optional (hkl), 2θ clip, and Angular FWHM — not manual K sliders.
    """
    if validation is None:
        validation = validate_scan_config(config)
    if not validation.ok:
        raise ValueError("; ".join(validation.errors))

    wavelength = resolve_wavelength(config.wavelength)
    k0 = xray_wavenumber(wavelength)
    tt_min, tt_max = config.two_theta_clip
    q_inst = max_q_magnitude(wavelength, tt_max)
    pad = _pad_from_fwhm(wavelength, config.fwhm_deg)
    dk_target = max(pad / 2.5, q_inst / 80.0)

    q_points: list[np.ndarray] = []
    try:
        for angles in generate_trajectory(config):
            q_points.append(q_crystal_at_angles(angles, wavelength, config.r_mount))
    except ValueError:
        pass

    if config.hkl and config.reciprocal_matrix is not None:
        g = crystal_g_vector(config.hkl, config.reciprocal_matrix)
        if float(np.linalg.norm(g)) > 1e-12:
            q_points.append(np.asarray(g, dtype=float))

    if not q_points:
        # Fallback: symmetric box from instrument limit
        limit = q_inst * 1.05
        return KGridSpec(
            -limit,
            limit,
            -limit,
            limit,
            0.0,
            limit,
            _steps_for_span(2 * limit, dk_target),
            _steps_for_span(limit, dk_target),
            summary="Auto K grid: instrument limit q_max={:.3f} A^-1.".format(q_inst),
        )

    stack = np.vstack(q_points)
    min_in_plane = max(pad, 0.25)
    kx_lo = float(np.min(stack[:, 0])) - pad
    kx_hi = float(np.max(stack[:, 0])) + pad
    ky_lo = float(np.min(stack[:, 1])) - pad
    ky_hi = float(np.max(stack[:, 1])) + pad
    kz_lo = float(np.min(stack[:, 2])) - pad
    kz_hi = float(np.max(stack[:, 2])) + pad

    if kx_hi - kx_lo < min_in_plane:
        mid = 0.5 * (kx_hi + kx_lo)
        kx_lo, kx_hi = mid - 0.5 * min_in_plane, mid + 0.5 * min_in_plane
    if ky_hi - ky_lo < min_in_plane:
        mid = 0.5 * (ky_hi + ky_lo)
        ky_lo, ky_hi = mid - 0.5 * min_in_plane, mid + 0.5 * min_in_plane
    if kz_hi - kz_lo < pad:
        mid = 0.5 * (kz_hi + kz_lo)
        kz_lo, kz_hi = mid - 0.5 * pad, mid + 0.5 * pad

    # Clip to instrument reach (sphere radius ~ 2 k0)
    q_cap = 2.0 * k0 * 1.02
    kx_lo = max(-q_cap, kx_lo)
    kx_hi = min(q_cap, kx_hi)
    ky_lo = max(-q_cap, ky_lo)
    ky_hi = min(q_cap, ky_hi)
    kz_lo = max(-q_cap, kz_lo)
    kz_hi = min(q_cap, kz_hi)

    if config.mode == "theta_2theta" and kz_lo < 0:
        kz_lo = 0.0

    n_para = _steps_for_span(max(kx_hi - kx_lo, ky_hi - ky_lo), dk_target)
    n_perp = _steps_for_span(kz_hi - kz_lo, dk_target, min_steps=12)

    summary = (
        "Auto K grid from scan + FWHM={fwhm:.3f} deg: "
        "Kx=[{kx0:.2f},{kx1:.2f}], Ky=[{ky0:.2f},{ky1:.2f}], "
        "Kz=[{kz0:.2f},{kz1:.2f}], steps={np}x{nz} (dk~{dk:.3f} A^-1)."
    ).format(
        fwhm=config.fwhm_deg,
        kx0=kx_lo,
        kx1=kx_hi,
        ky0=ky_lo,
        ky1=ky_hi,
        kz0=kz_lo,
        kz1=kz_hi,
        np=n_para,
        nz=n_perp,
        dk=dk_target,
    )
    return KGridSpec(
        kx_lo, kx_hi, ky_lo, ky_hi, kz_lo, kz_hi, n_para, n_perp, summary
    )
