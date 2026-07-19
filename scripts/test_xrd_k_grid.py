"""Tests for automatic XRD K-grid sizing."""

from __future__ import annotations

import os
import sys
import unittest

SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC_DIR))

from four_circle_geometry import FourCircleAngles
from structure_io import get_mount_matrix, load_ordered_structure, orient_crystal_axes
from xrd_constraints import ScanConfig, ScanRange, validate_scan_config
from xrd_k_grid import auto_k_grid_for_scan

DYBCO_CIF = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "cif_samples",
        "orthorhombic_O7",
        "DyBa2Cu3O7_mp-622105.cif",
    )
)


class TestAutoKGrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = orient_crystal_axes(load_ordered_structure(DYBCO_CIF), "c", "a")
        cls.reciprocal = cls.structure.lattice.reciprocal_lattice.matrix
        cls.r_mount = get_mount_matrix(cls.structure)
        cls.wavelength = 1.54184

    def test_theta_2theta_grid_covers_005(self):
        config = ScanConfig(
            mode="theta_2theta",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0, 0, 0, 38.4),
            primary=ScanRange(10.0, 80.0, 2.0),
            two_theta_clip=(10.0, 80.0),
            hkl=(0, 0, 5),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            fwhm_deg=0.1,
        )
        validation = validate_scan_config(config)
        spec = auto_k_grid_for_scan(config, validation)
        self.assertGreaterEqual(spec.kz_max, 2.5)
        self.assertGreaterEqual(spec.n_perp, 12)
        self.assertGreater(spec.n_para, 0)
        self.assertIn("Auto K grid", spec.summary)

    def test_fwhm_increases_padding(self):
        base = ScanConfig(
            mode="omega",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0, 0, 0, 38.4),
            primary=ScanRange(-1.0, 1.0, 0.5),
            two_theta_clip=(10.0, 80.0),
            hkl=(0, 0, 5),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            fwhm_deg=0.02,
        )
        wide = ScanConfig(
            mode="omega",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0, 0, 0, 38.4),
            primary=ScanRange(-1.0, 1.0, 0.5),
            two_theta_clip=(10.0, 80.0),
            hkl=(0, 0, 5),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            fwhm_deg=0.5,
        )
        spec_narrow = auto_k_grid_for_scan(base, validate_scan_config(base))
        spec_wide = auto_k_grid_for_scan(wide, validate_scan_config(wide))
        span_n = (spec_narrow.kz_max - spec_narrow.kz_min)
        span_w = (spec_wide.kz_max - spec_wide.kz_min)
        self.assertGreaterEqual(span_w, span_n)


if __name__ == "__main__":
    unittest.main()
