"""Verify Peng electron form factors against abTEM coefficient tables."""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from electron_form_factors import (  # noqa: E402
    PengFormFactorTable,
    cif_species_to_peng_ionic_key,
    peng_form_factor,
)


def main():
    table = PengFormFactorTable.load(os.path.join(ROOT, "src", "files"))
    coeffs_si = table.neutral["Si"]
    fe_zero = float(peng_form_factor(0.0, coeffs_si))
    fe_q2 = float(peng_form_factor(4.0, coeffs_si))
    assert abs(fe_zero - 5.8268) < 0.01, fe_zero
    assert abs(fe_q2 - 3.2637) < 0.02, fe_q2

    assert cif_species_to_peng_ionic_key("Si4+") == "Si++++"
    assert cif_species_to_peng_ionic_key("O2-") == "O--"
    assert cif_species_to_peng_ionic_key("Fe3+") == "Fe+++"

    fe_si_neutral = float(table.fe(0.0, 0.0, 0.0, "Si"))
    fe_si_ion = float(table.fe(0.0, 0.0, 0.0, "Si4+"))
    assert fe_si_neutral != fe_si_ion

    coeffs_unknown, warning = table.resolve_coefficients("Zz9+")
    assert coeffs_unknown is None
    assert warning is not None

    coeffs_fallback, warning = table.resolve_coefficients("Si9+")
    assert warning is not None and "using neutral" in warning
    assert np.allclose(coeffs_fallback, table.neutral["Si"])

    print("Peng form factor verification passed.")
    print("Si f_e(s=0) = {:.4f}".format(fe_zero))
    print("Si f_e(|Q|=2) = {:.4f}".format(fe_q2))
    print("Si4+ f_e(s=0) = {:.4f}".format(fe_si_ion))


if __name__ == "__main__":
    main()
