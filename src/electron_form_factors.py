import json
import os
import re

import numpy as np


def _base_element_symbol(specie):
    """Strip legacy valence-shell suffixes (e.g. Sival -> Si)."""
    if specie.endswith("val"):
        return specie[:-3]
    if len(specie) > 2 and specie[-1] == "v" and specie[0].isupper():
        return specie[:-1]
    return specie


def _parse_charge_species(specie):
    """Return (element, charge_sign, charge_magnitude) or neutral (element, None, 0)."""
    cleaned = _base_element_symbol(specie)
    match = re.match(r"^([A-Z][a-z]?)(\d*)([+-])$", cleaned)
    if not match:
        element = re.match(r"^([A-Z][a-z]?)", cleaned)
        if element:
            return element.group(1), None, 0
        return cleaned, None, 0
    element, count, sign = match.groups()
    magnitude = int(count) if count else 1
    return element, sign, magnitude


def cif_species_to_peng_ionic_key(specie):
    """Map CIF-style species (Si4+, O2-) to abTEM peng_ionic keys (Si++++, O--)."""
    element, sign, magnitude = _parse_charge_species(specie)
    if sign is None:
        return None
    repeat = sign * magnitude
    return element + repeat


def peng_form_factor(Q2, coeffs):
    """Peng 5-term Gaussian: fe(s) = sum ai * exp(-bi * s^2), s = |Q|/(4*pi)."""
    a = np.asarray(coeffs[0], dtype=float)
    b = np.asarray(coeffs[1], dtype=float)
    s2 = np.asarray(Q2, dtype=float) / (16.0 * np.pi**2)
    fe = np.zeros_like(s2, dtype=float)
    for index in range(len(a)):
        fe += a[index] * np.exp(-b[index] * s2)
    return fe


class PengFormFactorTable:
    """Peng 1996/1998 electron atomic scattering factors (abTEM coefficient tables)."""

    def __init__(self, neutral, ionic):
        self.neutral = neutral
        self.ionic = ionic
        self._warned_fallback = set()

    @classmethod
    def load(cls, files_dir):
        neutral_path = os.path.join(files_dir, "peng_high.json")
        ionic_path = os.path.join(files_dir, "peng_ionic.json")
        with open(neutral_path, encoding="utf-8") as handle:
            neutral = json.load(handle)
        with open(ionic_path, encoding="utf-8") as handle:
            ionic = json.load(handle)
        return cls(neutral, ionic)

    def resolve_coefficients(self, specie):
        """
        Return (coeffs, warning_message).
        coeffs shape (2, 5); warning_message is None or a fallback notice.
        """
        element, sign, _magnitude = _parse_charge_species(specie)
        if sign is not None:
            ionic_key = cif_species_to_peng_ionic_key(specie)
            if ionic_key in self.ionic:
                return np.asarray(self.ionic[ionic_key], dtype=float), None
            if element in self.neutral:
                warning = (
                    "No Peng ionic coefficients for {!r}; "
                    "using neutral {!r} scattering factors.".format(specie, element)
                )
                return np.asarray(self.neutral[element], dtype=float), warning
            return None, "No Peng scattering coefficients for {!r}.".format(specie)

        lookup = element
        if lookup in self.neutral:
            return np.asarray(self.neutral[lookup], dtype=float), None
        return None, "No Peng scattering coefficients for {!r}.".format(specie)

    def fe(self, Kx, Ky, Kz, specie):
        """Electron form factor array on the K grid."""
        Q2 = (
            np.multiply(Kx, Kx)
            + np.multiply(Ky, Ky)
            + np.multiply(Kz, Kz)
        )
        coeffs, warning = self.resolve_coefficients(specie)
        if coeffs is None:
            raise KeyError(warning)
        return peng_form_factor(Q2, coeffs)

    def consume_warning(self, specie):
        """Return a one-time fallback warning for specie, or None."""
        coeffs, warning = self.resolve_coefficients(specie)
        if warning is None or specie in self._warned_fallback:
            return None
        if coeffs is not None and "using neutral" in warning:
            self._warned_fallback.add(specie)
            return warning
        return None
