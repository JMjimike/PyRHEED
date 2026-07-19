"""XRD continuous scan engine: 1D spectra and 2D maps (chi-omega, RSM)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from four_circle_geometry import (
    FourCircleAngles,
    q_lab_to_film_coords,
    scattering_vector_lab,
)
from probe_geometry import XRDGeometry
from xrd_constraints import ScanConfig, validate_scan_config


@dataclass
class ScanResult:
    mode: str
    dimension: int
    x_axis: np.ndarray
    intensity: np.ndarray
    y_axis: np.ndarray | None = None
    grid: np.ndarray | None = None
    qx_axis: np.ndarray | None = None
    qz_axis: np.ndarray | None = None
    angles: list[FourCircleAngles] = field(default_factory=list)
    x_label: str = ""
    y_label: str = ""
    axis_display: str = "scan"


def _omega_2theta_offset(fixed: FourCircleAngles) -> float:
    return fixed.two_theta - 2.0 * fixed.omega


def _apply_coupled(omega: float, offset: float) -> float:
    return 2.0 * omega + offset


def generate_trajectory(config: ScanConfig) -> list[FourCircleAngles]:
    validation = validate_scan_config(config)
    if not validation.ok:
        raise ValueError("; ".join(validation.errors))
    primary = validation.primary or config.primary
    secondary = validation.secondary
    fixed = config.fixed
    mode = config.mode
    offset = _omega_2theta_offset(fixed)
    trajectory: list[FourCircleAngles] = []

    if mode == "theta_2theta":
        for two_theta in primary.values():
            trajectory.append(
                FourCircleAngles(fixed.omega, fixed.chi, fixed.phi, two_theta)
            )
    elif mode == "omega_2theta":
        for omega in primary.values():
            trajectory.append(
                FourCircleAngles(
                    omega, fixed.chi, fixed.phi, _apply_coupled(omega, offset)
                )
            )
    elif mode == "omega":
        for omega in primary.values():
            trajectory.append(
                FourCircleAngles(omega, fixed.chi, fixed.phi, fixed.two_theta)
            )
    elif mode == "chi":
        for chi in primary.values():
            trajectory.append(
                FourCircleAngles(fixed.omega, chi, fixed.phi, fixed.two_theta)
            )
    elif mode == "phi":
        for phi in primary.values():
            trajectory.append(
                FourCircleAngles(fixed.omega, fixed.chi, phi, fixed.two_theta)
            )
    elif mode == "chi_omega":
        if secondary is None:
            raise ValueError("chi_omega requires secondary scan range.")
        for chi in primary.values():
            for omega in secondary.values():
                trajectory.append(
                    FourCircleAngles(omega, chi, fixed.phi, fixed.two_theta)
                )
    elif mode == "rsm":
        if config.coupling == "omega_2theta":
            for omega in primary.values():
                trajectory.append(
                    FourCircleAngles(
                        omega,
                        fixed.chi,
                        fixed.phi,
                        _apply_coupled(omega, offset),
                    )
                )
        else:
            if secondary is None:
                raise ValueError("Independent RSM requires secondary 2theta range.")
            for omega in primary.values():
                for two_theta in secondary.values():
                    trajectory.append(
                        FourCircleAngles(
                            omega, fixed.chi, fixed.phi, two_theta
                        )
                    )
    else:
        raise ValueError("Unknown scan mode: {!r}".format(mode))
    return trajectory


def integrate_point(intensity, Kx, Ky, Kz, geometry: XRDGeometry) -> float:
    weighted = geometry.apply_broadening(intensity, Kx, Ky, Kz)
    return float(np.sum(weighted))


def _build_geometry(config: ScanConfig, angles: FourCircleAngles) -> XRDGeometry:
    return XRDGeometry(
        wavelength_angstrom=float(config.wavelength),
        fwhm_deg=float(getattr(config, "fwhm_deg", 0.0)),
        omega_deg=angles.omega,
        chi_deg=angles.chi,
        phi_deg=angles.phi,
        two_theta_deg=angles.two_theta,
        r_mount=config.r_mount,
        apply_lorentz=True,
    )


def _q_for_angles(config: ScanConfig, angles: FourCircleAngles):
    q_lab = scattering_vector_lab(
        angles.omega,
        angles.chi,
        angles.phi,
        angles.two_theta,
        config.wavelength,
        config.r_mount,
    )
    return q_lab_to_film_coords(q_lab, config.r_mount)


def run_scan(
    intensity,
    Kx,
    Ky,
    Kz,
    config: ScanConfig,
    progress_callback=None,
    abort_flag=None,
) -> ScanResult:
    """Integrate Ewald-weighted intensity along a scan trajectory."""
    validation = validate_scan_config(config)
    if not validation.ok:
        raise ValueError("; ".join(validation.errors))

    trajectory = generate_trajectory(config)
    if not trajectory:
        raise ValueError("Scan trajectory is empty.")

    values: list[float] = []
    for index, angles in enumerate(trajectory):
        if abort_flag is not None and abort_flag():
            raise RuntimeError("Scan aborted.")
        geometry = _build_geometry(config, angles)
        values.append(integrate_point(intensity, Kx, Ky, Kz, geometry))
        if progress_callback is not None:
            progress_callback(index + 1, len(trajectory))

    mode = config.mode
    primary = validation.primary or config.primary
    secondary = validation.secondary or config.secondary

    if mode == "chi_omega":
        x_vals = np.array(primary.values(), dtype=float)
        y_vals = np.array(secondary.values(), dtype=float)
        grid = np.array(values, dtype=float).reshape(len(x_vals), len(y_vals))
        qx = np.zeros_like(grid)
        qz = np.zeros_like(grid)
        idx = 0
        for i in range(len(x_vals)):
            for j in range(len(y_vals)):
                qx[i, j], qz[i, j] = _q_for_angles(config, trajectory[idx])
                idx += 1
        return ScanResult(
            mode=mode,
            dimension=2,
            x_axis=x_vals,
            y_axis=y_vals,
            grid=grid,
            qx_axis=qx,
            qz_axis=qz,
            intensity=np.array([], dtype=float),
            angles=trajectory,
            x_label="Chi (deg)",
            y_label="Omega (deg)",
            axis_display="chi_omega",
        )

    if mode == "rsm":
        if config.coupling == "omega_2theta":
            x_vals = np.array(primary.values(), dtype=float)
            intensity_arr = np.array(values, dtype=float)
            qx_list = []
            qz_list = []
            for angles in trajectory:
                qx, qz = _q_for_angles(config, angles)
                qx_list.append(qx)
                qz_list.append(qz)
            return ScanResult(
                mode=mode,
                dimension=1,
                x_axis=x_vals,
                intensity=intensity_arr,
                qx_axis=np.array(qx_list, dtype=float),
                qz_axis=np.array(qz_list, dtype=float),
                angles=trajectory,
                x_label="Omega (deg)",
                axis_display="scan",
            )

        x_vals = np.array(primary.values(), dtype=float)
        y_vals = np.array(secondary.values(), dtype=float)
        grid = np.array(values, dtype=float).reshape(len(x_vals), len(y_vals))
        qx = np.zeros_like(grid)
        qz = np.zeros_like(grid)
        idx = 0
        for i in range(len(x_vals)):
            for j in range(len(y_vals)):
                qx[i, j], qz[i, j] = _q_for_angles(config, trajectory[idx])
                idx += 1
        return ScanResult(
            mode=mode,
            dimension=2,
            x_axis=x_vals,
            y_axis=y_vals,
            grid=grid,
            qx_axis=qx,
            qz_axis=qz,
            intensity=np.array([], dtype=float),
            angles=trajectory,
            x_label="Omega (deg)",
            y_label="2theta (deg)",
            axis_display="omega_2theta",
        )

    if mode == "omega_2theta":
        x_vals = np.array(
            [_apply_coupled(om, _omega_2theta_offset(config.fixed)) for om in primary.values()],
            dtype=float,
        )
    else:
        x_vals = np.array(primary.values(), dtype=float)

    labels = {
        "theta_2theta": "2theta (deg)",
        "omega_2theta": "2theta (deg)",
        "omega": "Omega (deg)",
        "chi": "Chi (deg)",
        "phi": "Phi (deg)",
    }
    return ScanResult(
        mode=mode,
        dimension=1,
        x_axis=x_vals,
        intensity=np.array(values, dtype=float),
        angles=trajectory,
        x_label=labels.get(mode, "Scan"),
        axis_display="scan",
    )


def _column_from_label(label: str) -> str:
    key = label.lower().replace("(deg)", "deg").replace("(a^-1)", "inv_a")
    key = key.replace(" ", "_").strip("_")
    key = key.replace("2theta", "two_theta")
    return key or "scan_angle"


def scan_plot_kwargs(result: ScanResult) -> dict:
    title = "XRD {} scan".format(result.mode.replace("_", "-"))
    if result.dimension == 1:
        return {
            "window_title": title,
            "spectrum_title": title,
            "spectrum_x_label": result.x_label or "Scan angle (deg)",
            "spectrum_y_label": "Intensity (arb. units)",
        }
    return {
        "window_title": title,
        "spectrum_title": title,
        "spectrum_x_label": result.x_label,
        "spectrum_y_label": result.y_label,
    }


def export_scan_csv(result: ScanResult, path: str):
    """Write scan result to CSV."""
    import pandas as pd

    x_col = _column_from_label(result.x_label)
    if result.dimension == 1:
        rows = []
        for x, inten, angles in zip(result.x_axis, result.intensity, result.angles):
            row = {
                x_col: x,
                "intensity": inten,
                "omega": angles.omega,
                "chi": angles.chi,
                "phi": angles.phi,
                "two_theta": angles.two_theta,
            }
            if result.qx_axis is not None and result.qz_axis is not None:
                idx = len(rows)
                row["qx"] = float(result.qx_axis[idx])
                row["qz"] = float(result.qz_axis[idx])
            rows.append(row)
        pd.DataFrame(rows).to_csv(path, index=False)
        return

    y_col = _column_from_label(result.y_label)
    rows = []
    for i, x in enumerate(result.x_axis):
        for j, y in enumerate(result.y_axis):
            angle_index = i * len(result.y_axis) + j
            angles = result.angles[angle_index]
            row = {
                x_col: x,
                y_col: y,
                "intensity": result.grid[i, j],
                "omega": angles.omega,
                "chi": angles.chi,
                "phi": angles.phi,
                "two_theta": angles.two_theta,
            }
            if result.qx_axis is not None and result.qz_axis is not None:
                if result.qx_axis.ndim == 2:
                    row["qx"] = result.qx_axis[i, j]
                    row["qz"] = result.qz_axis[i, j]
                else:
                    row["qx"] = result.qx_axis[angle_index]
                    row["qz"] = result.qz_axis[angle_index]
            rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
