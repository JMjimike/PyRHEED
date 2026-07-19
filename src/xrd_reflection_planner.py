"""Reflection planner: map Miller indices to experimental XRD settings."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from four_circle_geometry import (
    FourCircleAngles,
    crystal_g_vector,
    g_vector_in_lab,
    is_symmetric_reflection,
    pick_recommended_solution,
    recommend_scan_mode,
    solve_bragg_angles,
)
from xrd_constraints import ScanRange, suggest_scan_ranges
from xrd_ewald import resolve_wavelength, xray_wavenumber


@dataclass
class ReflectionPlan:
    hkl: tuple[int, int, int]
    d_hkl: float
    theta_b_deg: float
    two_theta_deg: float
    omega_deg: float
    chi_deg: float
    phi_deg: float
    intensity: float | None
    scan_mode: str
    accessible: bool
    alt_solutions: list[FourCircleAngles] = field(default_factory=list)
    suggested_scan: dict = field(default_factory=dict)
    note: str = ""

    def to_row(self):
        return {
            "h": self.hkl[0],
            "k": self.hkl[1],
            "l": self.hkl[2],
            "d": self.d_hkl,
            "theta_B": self.theta_b_deg,
            "two_theta": self.two_theta_deg,
            "omega": self.omega_deg,
            "chi": self.chi_deg,
            "phi": self.phi_deg,
            "intensity": self.intensity,
            "scan_mode": self.scan_mode,
            "accessible": self.accessible,
            "note": self.note,
        }


def bragg_angles_from_d(d_hkl, wavelength_angstrom):
    wavelength = float(wavelength_angstrom)
    if d_hkl <= 0 or wavelength > 2.0 * d_hkl:
        return None, None
    theta_b = math.degrees(math.asin(wavelength / (2.0 * d_hkl)))
    return theta_b, 2.0 * theta_b


def plan_reflection(
    hkl,
    reciprocal_matrix,
    wavelength="CuKa",
    r_mount=None,
    current_angles: FourCircleAngles | None = None,
    intensity: float | None = None,
    full_search=False,
):
    """Build a ReflectionPlan for a single (hkl) reflection."""
    hkl = tuple(int(v) for v in hkl)
    wavelength_a = resolve_wavelength(wavelength)
    g_crystal = crystal_g_vector(hkl, reciprocal_matrix)
    g_norm = float(np.linalg.norm(g_crystal))
    if g_norm < 1e-12:
        return ReflectionPlan(
            hkl=hkl,
            d_hkl=0.0,
            theta_b_deg=0.0,
            two_theta_deg=0.0,
            omega_deg=0.0,
            chi_deg=0.0,
            phi_deg=0.0,
            intensity=intensity,
            scan_mode="none",
            accessible=False,
            note="Zero-length reciprocal vector.",
        )

    d_hkl = 2.0 * np.pi / g_norm
    theta_b, two_theta = bragg_angles_from_d(d_hkl, wavelength_a)
    if theta_b is None:
        return ReflectionPlan(
            hkl=hkl,
            d_hkl=d_hkl,
            theta_b_deg=0.0,
            two_theta_deg=0.0,
            omega_deg=0.0,
            chi_deg=0.0,
            phi_deg=0.0,
            intensity=intensity,
            scan_mode="none",
            accessible=False,
            note="Wavelength exceeds 2d; reflection not reachable.",
        )

    solutions = solve_bragg_angles(
        hkl,
        reciprocal_matrix,
        wavelength_a,
        r_mount=r_mount,
        current_angles=current_angles,
        full_search=full_search,
    )
    g_lab_zero = g_vector_in_lab(hkl, reciprocal_matrix, 0.0, 0.0, 0.0, r_mount)
    if is_symmetric_reflection(g_lab_zero) and theta_b is not None:
        recommended = FourCircleAngles(0.0, 0.0, 0.0, two_theta)
        solutions = [recommended] + [
            sol for sol in solutions if sol.as_tuple() != recommended.as_tuple()
        ]
    else:
        recommended = pick_recommended_solution(solutions, current_angles)
    if recommended is None:
        g_lab = g_vector_in_lab(hkl, reciprocal_matrix, 0.0, 0.0, 0.0, r_mount)
        return ReflectionPlan(
            hkl=hkl,
            d_hkl=d_hkl,
            theta_b_deg=theta_b,
            two_theta_deg=two_theta,
            omega_deg=0.0,
            chi_deg=0.0,
            phi_deg=0.0,
            intensity=intensity,
            scan_mode=recommend_scan_mode(g_lab),
            accessible=False,
            note="Bragg condition not satisfied with current geometry search.",
        )

    g_lab = g_vector_in_lab(
        hkl,
        reciprocal_matrix,
        recommended.omega,
        recommended.chi,
        recommended.phi,
        r_mount,
    )
    alts = [sol for sol in solutions if sol.as_tuple() != recommended.as_tuple()]
    scan_mode = recommend_scan_mode(g_lab)
    primary, secondary, ui_mode = suggest_scan_ranges(
        scan_mode, recommended.two_theta, two_theta_clip=(0.0, 179.0)
    )
    suggested = {
        "scan_mode": ui_mode,
        "primary": (primary.start, primary.stop, primary.step),
    }
    if secondary is not None:
        suggested["secondary"] = (secondary.start, secondary.stop, secondary.step)
    return ReflectionPlan(
        hkl=hkl,
        d_hkl=d_hkl,
        theta_b_deg=theta_b,
        two_theta_deg=recommended.two_theta,
        omega_deg=recommended.omega,
        chi_deg=recommended.chi,
        phi_deg=recommended.phi,
        intensity=intensity,
        scan_mode=scan_mode,
        accessible=True,
        alt_solutions=alts,
        suggested_scan=suggested,
        note="Recommended setting closest to current angles."
        if current_angles
        else "Primary Bragg solution.",
    )


def enumerate_reflections(
    structure,
    wavelength="CuKa",
    h_max=5,
    k_max=5,
    l_max=5,
    r_mount=None,
    intensity_fn=None,
    full_search=False,
    two_theta_min=0.0,
    two_theta_max=180.0,
):
    """Yield ReflectionPlan objects for all integer hkl within limits."""
    reciprocal = structure.lattice.reciprocal_lattice.matrix
    wavelength_a = resolve_wavelength(wavelength)
    seen_hkl = set()
    for h in range(-h_max, h_max + 1):
        for k in range(-k_max, k_max + 1):
            for l in range(-l_max, l_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                hkl_key = (abs(h), abs(k), abs(l))
                if hkl_key in seen_hkl:
                    continue
                seen_hkl.add(hkl_key)
                rep_h, rep_k, rep_l = h, k, l
                if rep_l < 0 or (rep_l == 0 and rep_k < 0) or (
                    rep_l == 0 and rep_k == 0 and rep_h < 0
                ):
                    rep_h, rep_k, rep_l = -rep_h, -rep_k, -rep_l

                g_norm = float(np.linalg.norm(crystal_g_vector((h, k, l), reciprocal)))
                if g_norm < 1e-12:
                    continue
                d_hkl = 2.0 * np.pi / g_norm
                theta_b, two_theta = bragg_angles_from_d(d_hkl, wavelength_a)
                if theta_b is None:
                    continue
                if not (two_theta_min <= two_theta <= two_theta_max):
                    continue

                intensity = intensity_fn(rep_h, rep_k, rep_l) if intensity_fn else None
                yield plan_reflection(
                    (rep_h, rep_k, rep_l),
                    reciprocal,
                    wavelength=wavelength,
                    r_mount=r_mount,
                    intensity=intensity,
                    full_search=full_search,
                )


def accessible_peak_table(plans, min_intensity=0.0):
    rows = []
    for plan in plans:
        if not plan.accessible:
            continue
        if plan.intensity is not None and plan.intensity < min_intensity:
            continue
        rows.append(plan.to_row())
    rows.sort(key=lambda row: row["two_theta"])
    return rows
