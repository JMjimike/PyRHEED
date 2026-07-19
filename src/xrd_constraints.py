"""Wavelength and geometry validation for XRD scan configurations."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from four_circle_geometry import FourCircleAngles, crystal_g_vector
from xrd_ewald import resolve_wavelength


@dataclass
class ScanRange:
    start: float
    stop: float
    step: float

    def values(self):
        if self.step <= 0:
            return []
        start, stop, step = self.start, self.stop, self.step
        if start > stop:
            start, stop = stop, start
        values = []
        value = start
        while value <= stop + 1e-9:
            values.append(float(value))
            value += step
        return values


@dataclass
class ScanConfig:
    mode: str
    wavelength: float
    fixed: FourCircleAngles
    primary: ScanRange
    secondary: ScanRange | None = None
    coupling: str = "none"  # "none" | "omega_2theta"
    two_theta_clip: tuple[float, float] = (0.0, 180.0)
    hkl: tuple[int, int, int] | None = None
    reciprocal_matrix: np.ndarray | None = None
    r_mount: np.ndarray | None = None
    fwhm_deg: float = 0.0


@dataclass
class ValidationResult:
    ok: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    primary: ScanRange | None = None
    secondary: ScanRange | None = None


def _clip_range(scan_range: ScanRange, low: float, high: float) -> ScanRange:
    start = max(low, min(scan_range.start, scan_range.stop))
    stop = min(high, max(scan_range.start, scan_range.stop))
    if stop < start:
        start, stop = stop, start
    return ScanRange(start, stop, scan_range.step)


def max_accessible_two_theta(wavelength_angstrom, two_theta_clip=(0.0, 180.0)):
    """Upper 2θ limit: instrument clip capped by elastic |Q| <= 2k0 (θ_B <= 90°)."""
    instrument_max = float(two_theta_clip[1])
    # Bragg kinematics: |Q| = 2 k0 sin(θ); maximum at θ = 90° gives 2θ = 180°.
    physical_max = 179.99
    return min(instrument_max, physical_max)


def bragg_two_theta(hkl, reciprocal_matrix, wavelength_angstrom):
    if hkl is None or reciprocal_matrix is None:
        return None
    g = crystal_g_vector(hkl, reciprocal_matrix)
    g_norm = float(np.linalg.norm(g))
    if g_norm < 1e-12:
        return None
    d_hkl = 2.0 * np.pi / g_norm
    wavelength = float(wavelength_angstrom)
    if wavelength > 2.0 * d_hkl:
        return None
    theta_b = math.degrees(math.asin(wavelength / (2.0 * d_hkl)))
    return 2.0 * theta_b


def validate_scan_config(config: ScanConfig) -> ValidationResult:
    result = ValidationResult(ok=True)
    wavelength = resolve_wavelength(config.wavelength)
    tt_min, tt_max = config.two_theta_clip
    tt_limit = max_accessible_two_theta(wavelength, config.two_theta_clip)

    if config.primary.step <= 0:
        result.errors.append("Primary scan step must be positive.")
    if config.mode in ("chi_omega", "rsm") and config.coupling != "omega_2theta":
        if config.secondary is None or config.secondary.step <= 0:
            result.errors.append("Secondary scan step must be positive for 2D scans.")

    if config.hkl and config.reciprocal_matrix is not None:
        bragg_tt = bragg_two_theta(
            config.hkl, config.reciprocal_matrix, wavelength
        )
        if bragg_tt is None:
            result.errors.append(
                "Target (hkl) is not reachable at the selected wavelength (lambda > 2d)."
            )
        elif config.mode in ("omega", "chi", "phi", "chi_omega", "rsm"):
            delta = abs(config.fixed.two_theta - bragg_tt)
            if delta > 10.0:
                result.warnings.append(
                    "Fixed 2theta differs from Bragg angle by {:.1f} deg.".format(delta)
                )

    primary = config.primary
    if config.mode == "theta_2theta":
        primary = _clip_range(primary, tt_min, min(tt_max, tt_limit))
    elif config.mode == "omega_2theta":
        primary = _clip_range(primary, -180.0, 180.0)
    elif config.mode == "omega":
        primary = _clip_range(primary, -180.0, 180.0)
    elif config.mode == "chi":
        primary = _clip_range(primary, 0.0, 180.0)
    elif config.mode == "phi":
        primary = _clip_range(primary, 0.0, 360.0)

    secondary = config.secondary
    if config.mode == "chi_omega" and secondary is not None:
        secondary = _clip_range(secondary, -180.0, 180.0)
        primary = _clip_range(primary, 0.0, 180.0)
    elif config.mode == "rsm" and config.coupling != "omega_2theta" and secondary is not None:
        secondary = _clip_range(secondary, tt_min, min(tt_max, tt_limit))
        primary = _clip_range(primary, -180.0, 180.0)

    if not primary.values():
        result.errors.append("Primary scan range is empty after clipping.")

    if config.mode in ("chi_omega", "rsm") and config.coupling != "omega_2theta":
        if secondary is None or not secondary.values():
            result.errors.append("Secondary scan range is empty after clipping.")

    for angle in (config.fixed.omega, config.fixed.chi, config.fixed.phi, config.fixed.two_theta):
        if not math.isfinite(angle):
            result.errors.append("Fixed angles must be finite.")

    result.primary = primary
    result.secondary = secondary
    result.ok = not result.errors
    return result


def suggest_scan_ranges(
    plan_scan_mode: str,
    two_theta_bragg: float,
    two_theta_clip=(0.0, 179.0),
):
    """Return suggested primary/secondary ranges from a reflection plan."""
    tt_min, tt_max = two_theta_clip
    width_1d = min(35.0, (tt_max - tt_min) * 0.5)
    width_rsm = 4.0
    if plan_scan_mode in ("theta_2theta", "phi_scan"):
        primary = ScanRange(
            max(tt_min, two_theta_bragg - width_1d),
            min(tt_max, two_theta_bragg + width_1d),
            0.05,
        )
        return primary, None, "theta_2theta"
    if plan_scan_mode == "phi_scan":
        primary = ScanRange(0.0, 360.0, 1.0)
        return primary, None, "phi"
    primary = ScanRange(-width_rsm, width_rsm, 0.02)
    secondary = ScanRange(
        max(tt_min, two_theta_bragg - width_rsm),
        min(tt_max, two_theta_bragg + width_rsm),
        0.02,
    )
    return primary, secondary, "rsm"
