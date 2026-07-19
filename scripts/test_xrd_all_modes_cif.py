"""Run every XRD scan mode on DyBCO CIF (end-to-end, no GUI)."""

from __future__ import annotations

import os
import sys
import time
import unittest

import numpy as np

SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC_DIR))

from four_circle_geometry import FourCircleAngles
from process import DiffractionPattern
from probe_geometry import auto_k_range
from scattering_factors import XRayFormFactorTable
from structure_io import get_mount_matrix, load_ordered_structure, orient_crystal_axes
from xrd_constraints import ScanConfig, ScanRange, validate_scan_config
from xrd_reflection_planner import plan_reflection
from xrd_scan import run_scan

DYBCO_CIF = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "cif_samples",
        "orthorhombic_O7",
        "DyBa2Cu3O7_mp-622105.cif",
    )
)


def build_atom_dict(structure):
    atoms = {}
    for site in structure:
        species = site.specie.symbol if site.is_ordered else list(site.species)[0].symbol
        atoms[tuple(site.coords)] = species
    return atoms


class TestAllScanModesOnDyBCO(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = orient_crystal_axes(load_ordered_structure(DYBCO_CIF), "c", "a")
        cls.reciprocal = cls.structure.lattice.reciprocal_lattice.matrix
        cls.r_mount = get_mount_matrix(cls.structure)
        cls.wavelength = 1.54184
        cls.xray = XRayFormFactorTable.load()
        cls.atoms = build_atom_dict(cls.structure)
        cls.plan_005 = plan_reflection(
            (0, 0, 5),
            cls.reciprocal,
            wavelength=cls.wavelength,
            r_mount=cls.r_mount,
            full_search=True,
        )
        k_limit = auto_k_range(cls.wavelength, 80.0)[1]
        grid = np.linspace(-0.5, 0.5, 8)
        z_grid = np.linspace(0.0, k_limit, 16)
        Kx, Ky, Kz = np.meshgrid(grid, grid, z_grid)
        worker = DiffractionPattern(Kx, Ky, Kz, cls.xray, cls.atoms, useCUDA=False)
        worker.run()
        cls.intensity = worker.intensity
        cls.Kx, cls.Ky, cls.Kz = Kx, Ky, Kz
        cls.fixed = FourCircleAngles(
            cls.plan_005.omega_deg,
            cls.plan_005.chi_deg,
            cls.plan_005.phi_deg,
            cls.plan_005.two_theta_deg,
        )

    def _run_mode(self, mode, **kwargs):
        defaults = dict(
            mode=mode,
            wavelength=self.wavelength,
            fixed=self.fixed,
            primary=ScanRange(10.0, 80.0, 2.0),
            two_theta_clip=(10.0, 80.0),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            hkl=(0, 0, 5),
            fwhm_deg=0.1,
        )
        defaults.update(kwargs)
        config = ScanConfig(**defaults)
        validation = validate_scan_config(config)
        self.assertTrue(validation.ok, validation.errors)
        start = time.time()
        result = run_scan(self.intensity, self.Kx, self.Ky, self.Kz, config)
        elapsed = time.time() - start
        if result.dimension == 1:
            peak = float(np.max(result.intensity))
            self.assertGreater(peak, 0.0, msg="{} produced zero signal".format(mode))
            print(
                "{:14s} 1D  n={:4d}  peak={:.3g}  t={:.2f}s".format(
                    mode, len(result.x_axis), peak, elapsed
                )
            )
        else:
            peak = float(np.max(result.grid))
            self.assertGreater(peak, 0.0, msg="{} produced zero signal".format(mode))
            print(
                "{:14s} 2D  shape={}  peak={:.3g}  t={:.2f}s".format(
                    mode, result.grid.shape, peak, elapsed
                )
            )
        return result

    def test_theta_2theta(self):
        self._run_mode(
            "theta_2theta",
            fixed=FourCircleAngles(0, 0, 0, 38.4),
            primary=ScanRange(10.0, 80.0, 2.0),
        )

    def test_omega_2theta(self):
        self._run_mode(
            "omega_2theta",
            primary=ScanRange(-2.0, 2.0, 0.5),
            coupling="omega_2theta",
        )

    def test_omega(self):
        self._run_mode("omega", primary=ScanRange(-1.0, 1.0, 0.5))

    def test_chi(self):
        self._run_mode("chi", primary=ScanRange(-2.0, 2.0, 1.0))

    def test_phi(self):
        self._run_mode("phi", primary=ScanRange(0.0, 360.0, 45.0))

    def test_chi_omega(self):
        self._run_mode(
            "chi_omega",
            primary=ScanRange(0.0, 4.0, 2.0),
            secondary=ScanRange(-1.0, 1.0, 1.0),
        )

    def test_rsm_coupled(self):
        self._run_mode(
            "rsm",
            primary=ScanRange(-2.0, 2.0, 1.0),
            coupling="omega_2theta",
        )

    def test_rsm_independent(self):
        self._run_mode(
            "rsm",
            primary=ScanRange(-1.0, 1.0, 1.0),
            secondary=ScanRange(36.0, 40.0, 1.0),
            coupling="none",
        )


if __name__ == "__main__":
    print("DyBCO CIF:", DYBCO_CIF)
    unittest.main(verbosity=2)
