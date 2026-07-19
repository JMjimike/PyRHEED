"""Unified atomic form factor tables for RHEED (electron) and XRD (X-ray)."""

from __future__ import annotations

import json
import os
from typing import Protocol, runtime_checkable

import numpy as np

from electron_form_factors import (
    PengFormFactorTable,
    _parse_charge_species,
    peng_form_factor,
)

__all__ = [
    "FormFactorTable",
    "PengFormFactorTable",
    "XRayFormFactorTable",
    "xray_form_factor",
]


@runtime_checkable
class FormFactorTable(Protocol):
    """Protocol for probe-specific atomic scattering factors on a K grid."""

    def fe(self, Kx, Ky, Kz, specie) -> np.ndarray:
        """Return scattering factor array matching the shape of Kx."""

    def consume_warning(self, specie):
        """Optional one-time warning for fallback coefficients; may return None."""


def xray_form_factor(Q2, atomic_number, coeffs):
    """Cromer-Mann X-ray form factor: f(s) = Z - 41.78214*s² * Σ a_i exp(-b_i s²)."""
    s2 = np.asarray(Q2, dtype=float) / (16.0 * np.pi**2)
    coeffs = np.asarray(coeffs, dtype=float)
    gaussian = np.zeros_like(s2, dtype=float)
    for index in range(coeffs.shape[0]):
        gaussian += coeffs[index, 0] * np.exp(-coeffs[index, 1] * s2)
    return atomic_number - 41.78214 * s2 * gaussian


class XRayFormFactorTable:
    """X-ray atomic scattering factors from pymatgen coefficient tables."""

    def __init__(self, scattering_params):
        self.scattering_params = scattering_params
        self._warned_fallback = set()
        self.debye_waller_factors: dict[str, float] = {}

    @classmethod
    def load(cls, pymatgen_diffraction_dir=None):
        if pymatgen_diffraction_dir is None:
            import pymatgen.analysis.diffraction as diffraction_pkg

            pymatgen_diffraction_dir = os.path.dirname(diffraction_pkg.__file__)
        json_path = os.path.join(
            pymatgen_diffraction_dir, "atomic_scattering_params.json"
        )
        with open(json_path, "rb") as handle:
            scattering_params = json.load(handle)
        return cls(scattering_params)

    def set_debye_waller_factors(self, factors: dict[str, float] | None):
        self.debye_waller_factors = factors or {}

    def resolve_element(self, specie):
        element, _sign, _magnitude = _parse_charge_species(specie)
        if element not in self.scattering_params:
            raise KeyError(
                "No X-ray scattering coefficients for {!r}.".format(specie)
            )
        return element

    def fe(self, Kx, Ky, Kz, specie):
        Q2 = Kx * Kx + Ky * Ky + Kz * Kz
        element = self.resolve_element(specie)
        from pymatgen.core import Element

        atomic_number = Element(element).Z
        coeffs = np.asarray(self.scattering_params[element], dtype=float)
        f = xray_form_factor(Q2, atomic_number, coeffs)
        b_factor = self.debye_waller_factors.get(element, 0.0)
        if b_factor:
            s2 = np.asarray(Q2, dtype=float) / (16.0 * np.pi**2)
            f = f * np.exp(-b_factor * s2)
        return f

    def consume_warning(self, specie):
        return None
