"""X-ray Ewald-sphere utilities for single-crystal XRD."""

import numpy as np

_GAUSSIAN_FWHM_TO_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))

# Common characteristic wavelengths in angstroms (pymatgen XRD convention).
WAVELENGTHS = {
    "CuKa": 1.54184,
    "CuKa1": 1.54056,
    "CuKa2": 1.54439,
    "MoKa": 0.71073,
    "MoKa1": 0.70930,
    "AgKa": 0.560885,
    "CrKa": 2.29100,
    "FeKa": 1.93735,
    "CoKa": 1.79026,
}


def resolve_wavelength(wavelength):
    """Return wavelength in angstroms from a float or preset keyword."""
    if isinstance(wavelength, str):
        if wavelength not in WAVELENGTHS:
            raise KeyError("Unknown wavelength preset: {!r}".format(wavelength))
        return WAVELENGTHS[wavelength]
    return float(wavelength)


def xray_wavenumber(wavelength_angstrom):
    """X-ray wavenumber k0 = 2π/λ in inverse angstroms."""
    return 2.0 * np.pi / float(wavelength_angstrom)


def incident_wavevector_z(k0):
    """Incident wavevector along -Z (laboratory frame)."""
    return np.array([0.0, 0.0, -k0], dtype=float)


def scattered_wavevector_y(two_theta_deg, k0):
    """Scattered wavevector for detector arm 2θ about +Y (XZ scattering plane)."""
    two_theta = np.deg2rad(two_theta_deg)
    return k0 * np.array([np.sin(two_theta), 0.0, np.cos(two_theta)], dtype=float)


def ewald_deviation(Kx, Ky, Kz, ki, k0):
    """Deviation from the elastic Ewald sphere: |ki + Q| - k0."""
    ks_mag = np.sqrt(
        (Kx + ki[0]) ** 2 + (Ky + ki[1]) ** 2 + (Kz + ki[2]) ** 2
    )
    return ks_mag - k0


def angular_sigma_delta(k0, fwhm_deg):
    """k-space broadening sigma from Gaussian angular FWHM in degrees."""
    if fwhm_deg <= 0.0:
        return 0.0
    sigma_rad = np.deg2rad(fwhm_deg) / _GAUSSIAN_FWHM_TO_SIGMA
    return k0 * sigma_rad


def apply_xrd_ewald_broadening(
    intensity,
    Kx,
    Ky,
    Kz,
    wavelength_angstrom,
    fwhm_deg=0.0,
    ki=None,
):
    """Apply Ewald-sphere constraint with optional angular broadening."""
    k0 = xray_wavenumber(wavelength_angstrom)
    if ki is None:
        ki = incident_wavevector_z(k0)
    delta = ewald_deviation(Kx, Ky, Kz, ki, k0)
    sigma_delta = angular_sigma_delta(k0, fwhm_deg)
    if sigma_delta <= 0.0:
        weight = (np.abs(delta) < 1e-9).astype(float)
    else:
        weight = np.exp(-0.5 * (delta / sigma_delta) ** 2)
    return np.real(intensity).astype(float) * weight


def project_xrd_ewald_to_xz(
    intensity,
    Kx,
    Ky,
    Kz,
    wavelength_angstrom,
    fwhm_deg=0.0,
    ki=None,
):
    """Ewald-broadened intensity summed over Ky onto a Kx-Kz map."""
    weighted = apply_xrd_ewald_broadening(
        intensity, Kx, Ky, Kz, wavelength_angstrom, fwhm_deg, ki=ki
    )
    return np.sum(weighted, axis=0)


def max_q_magnitude(wavelength_angstrom, two_theta_max_deg):
    """Estimate maximum |Q| reachable at a given 2θ limit."""
    k0 = xray_wavenumber(wavelength_angstrom)
    theta_max = np.deg2rad(two_theta_max_deg) / 2.0
    return 2.0 * k0 * np.sin(theta_max)


def lorentz_polarization_factor(two_theta_deg):
    """Lorentz-polarization factor for unpolarized reflection geometry."""
    theta = np.deg2rad(two_theta_deg) / 2.0
    sin_t = np.sin(theta)
    cos_t = np.cos(theta)
    if sin_t == 0.0 or cos_t == 0.0:
        return 0.0
    cos2t = np.cos(2.0 * theta)
    return (1.0 + cos2t**2) / (sin_t**2 * cos_t)
