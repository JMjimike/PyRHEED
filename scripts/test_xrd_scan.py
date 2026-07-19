"""Tests for XRD continuous scan engine."""

from __future__ import annotations

import os
import sys
import unittest

import numpy as np

SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC_DIR))

from four_circle_geometry import FourCircleAngles, ki_in_crystal_frame
from process import DiffractionPattern
from probe_geometry import XRDGeometry
from scattering_factors import XRayFormFactorTable
from structure_io import get_mount_matrix, load_ordered_structure, orient_crystal_axes
from xrd_constraints import ScanConfig, ScanRange, validate_scan_config
from xrd_reflection_planner import plan_reflection
from xrd_scan import generate_trajectory, run_scan

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
        coord = tuple(site.coords)
        atoms[coord] = species
    return atoms


class TestXRDScan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = orient_crystal_axes(load_ordered_structure(DYBCO_CIF), "c", "a")
        cls.reciprocal = cls.structure.lattice.reciprocal_lattice.matrix
        cls.r_mount = get_mount_matrix(cls.structure)
        cls.wavelength = 1.54184
        cls.xray = XRayFormFactorTable.load()
        cls.atoms = build_atom_dict(cls.structure)

    def _small_intensity_grid(self):
        grid = np.linspace(-0.5, 0.5, 6)
        z_grid = np.linspace(0.0, 3.5, 20)
        Kx, Ky, Kz = np.meshgrid(grid, grid, z_grid)
        worker = DiffractionPattern(
            Kx, Ky, Kz, self.xray, self.atoms, useCUDA=False
        )
        worker.run()
        return worker.intensity, Kx, Ky, Kz

    def _scan_config(self, **kwargs):
        defaults = dict(
            mode="theta_2theta",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0.0, 0.0, 0.0, 38.4),
            primary=ScanRange(35.0, 42.0, 0.2),
            two_theta_clip=(10.0, 80.0),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            fwhm_deg=0.1,
        )
        defaults.update(kwargs)
        return ScanConfig(**defaults)

    def test_ki_crystal_frame_differs_with_omega(self):
        ki0 = ki_in_crystal_frame(0.0, 0.0, 0.0, 38.4, self.wavelength, self.r_mount)
        ki1 = ki_in_crystal_frame(1.0, 0.0, 0.0, 38.4, self.wavelength, self.r_mount)
        self.assertGreater(np.linalg.norm(ki0 - ki1), 1e-6)

    def test_validate_theta_2theta_scan(self):
        config = ScanConfig(
            mode="theta_2theta",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0.0, 0.0, 0.0, 38.4),
            primary=ScanRange(10.0, 80.0, 5.0),
            two_theta_clip=(10.0, 80.0),
            reciprocal_matrix=self.reciprocal,
            r_mount=self.r_mount,
            hkl=(0, 0, 5),
        )
        result = validate_scan_config(config)
        self.assertTrue(result.ok)
        self.assertGreater(len(result.primary.values()), 0)

    def test_theta_2theta_scan_finds_005(self):
        intensity, Kx, Ky, Kz = self._small_intensity_grid()
        config = self._scan_config()
        result = run_scan(intensity, Kx, Ky, Kz, config)
        self.assertEqual(result.dimension, 1)
        self.assertGreater(float(np.max(result.intensity)), 0.0)
        peak_idx = int(np.argmax(result.intensity))
        self.assertAlmostEqual(result.x_axis[peak_idx], 38.4, delta=2.5)

    def test_rescan_differs_with_fixed_omega(self):
        intensity, Kx, Ky, Kz = self._small_intensity_grid()
        config_a = self._scan_config(
            mode="theta_2theta",
            fixed=FourCircleAngles(0.0, 0.0, 0.0, 38.4),
            primary=ScanRange(36.0, 40.0, 0.5),
        )
        config_b = self._scan_config(
            mode="theta_2theta",
            fixed=FourCircleAngles(2.0, 0.0, 0.0, 38.4),
            primary=ScanRange(36.0, 40.0, 0.5),
        )
        res_a = run_scan(intensity, Kx, Ky, Kz, config_a)
        res_b = run_scan(intensity, Kx, Ky, Kz, config_b)
        self.assertFalse(np.allclose(res_a.intensity, res_b.intensity))

    def test_chi_omega_2d_shape(self):
        intensity, Kx, Ky, Kz = self._small_intensity_grid()
        config = ScanConfig(
            mode="chi_omega",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0.0, 0.0, 0.0, 38.4),
            primary=ScanRange(0.0, 2.0, 1.0),
            secondary=ScanRange(-1.0, 1.0, 1.0),
            two_theta_clip=(10.0, 80.0),
            r_mount=self.r_mount,
            fwhm_deg=0.1,
        )
        result = run_scan(intensity, Kx, Ky, Kz, config)
        self.assertEqual(result.dimension, 2)
        self.assertEqual(result.grid.shape, (3, 3))

    def test_rsm_coupled_trajectory(self):
        config = ScanConfig(
            mode="rsm",
            wavelength=self.wavelength,
            fixed=FourCircleAngles(0.0, 0.0, 0.0, 38.4),
            primary=ScanRange(-1.0, 1.0, 0.5),
            coupling="omega_2theta",
            two_theta_clip=(10.0, 80.0),
            r_mount=self.r_mount,
        )
        traj = generate_trajectory(config)
        self.assertEqual(len(traj), 5)
        self.assertAlmostEqual(traj[0].two_theta, 36.4, places=3)

    def test_plan_reflection_005_still_works(self):
        plan = plan_reflection(
            (0, 0, 5),
            self.reciprocal,
            wavelength=self.wavelength,
            r_mount=self.r_mount,
            full_search=True,
        )
        self.assertTrue(plan.accessible)
        self.assertAlmostEqual(plan.two_theta_deg, 38.4, delta=0.5)


if __name__ == "__main__":
    unittest.main()
