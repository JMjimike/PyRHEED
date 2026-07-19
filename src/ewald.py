import numpy as np

_GAUSSIAN_FWHM_TO_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


def electron_wavenumber(energy_kev):
    """Electron wavenumber k0 in inverse Angstroms."""
    return 0.512 * np.sqrt(energy_kev * 1000.0)


def incident_wavevector_y(k0, theta_deg):
    """Incident wavevector for beam from +y with grazing angle theta (surface normal = z)."""
    theta = np.deg2rad(theta_deg)
    return np.array([0.0, -k0 * np.cos(theta), k0 * np.sin(theta)])


def ewald_deviation(Kx, Ky, Kz, ki, k0):
    """Deviation from the elastic Ewald sphere: |ki + Q| - k0."""
    ks_mag = np.sqrt(
        Kx * Kx + (Ky + ki[1]) ** 2 + (Kz + ki[2]) ** 2
    )
    return ks_mag - k0


def gaussian_energy_sigma_delta(k0, energy_kev, fwhm_ev):
    """k-space broadening sigma from Gaussian energy spread (FWHM in eV)."""
    if fwhm_ev <= 0.0:
        return 0.0
    energy_ev = energy_kev * 1000.0
    sigma_e = fwhm_ev / _GAUSSIAN_FWHM_TO_SIGMA
    return k0 * sigma_e / (2.0 * energy_ev)


def apply_ewald_broadening(
    intensity, Kx, Ky, Kz, energy_kev, fwhm_ev, theta_deg
):
    """Apply Ewald-sphere constraint with Gaussian energy broadening to intensity."""
    k0 = electron_wavenumber(energy_kev)
    ki = incident_wavevector_y(k0, theta_deg)
    delta = ewald_deviation(Kx, Ky, Kz, ki, k0)
    sigma_delta = gaussian_energy_sigma_delta(k0, energy_kev, fwhm_ev)
    if sigma_delta <= 0.0:
        weight = (np.abs(delta) < 1e-12).astype(float)
    else:
        weight = np.exp(-0.5 * (delta / sigma_delta) ** 2)
    return np.asarray(intensity, dtype=float) * weight


def project_ewald_to_xz(
    intensity, Kx, Ky, Kz, energy_kev, fwhm_ev, theta_deg
):
    """Ewald-broadened intensity summed over Ky onto a single Kx-Kz map."""
    weighted = apply_ewald_broadening(
        intensity, Kx, Ky, Kz, energy_kev, fwhm_ev, theta_deg
    )
    return np.sum(weighted, axis=0)
