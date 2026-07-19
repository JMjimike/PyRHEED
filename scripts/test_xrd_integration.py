"""End-to-end XRD feature tests on real CIF structures (no GUI)."""

from __future__ import annotations

import os
import sys
import tempfile
import time
import unittest

import numpy as np
import pandas as pd

SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC_DIR))

from process import DiffractionPattern
from probe_geometry import FourCircleGeometry, XRDGeometry, auto_k_range
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from scattering_factors import XRayFormFactorTable
from structure_io import (
    get_mount_matrix,
    load_ordered_structure,
    orient_crystal_axes,
)
from xrd_ewald import project_xrd_ewald_to_xz
from xrd_reflection_planner import (
    accessible_peak_table,
    enumerate_reflections,
    plan_reflection,
)

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


class TestXRDIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = orient_crystal_axes(
            load_ordered_structure(DYBCO_CIF), "c", "a"
        )
        cls.reciprocal = cls.structure.lattice.reciprocal_lattice.matrix
        cls.r_mount = get_mount_matrix(cls.structure)
        cls.wavelength = 1.54184
        cls.xray = XRayFormFactorTable.load()

    def test_enumerate_reflections_fast(self):
        start = time.time()
        plans = list(
            enumerate_reflections(
                self.structure,
                wavelength=self.wavelength,
                h_max=4,
                k_max=4,
                l_max=8,
                r_mount=self.r_mount,
                full_search=False,
                two_theta_min=10.0,
                two_theta_max=80.0,
            )
        )
        elapsed = time.time() - start
        rows = accessible_peak_table(plans)
        self.assertLess(elapsed, 30.0, "enumerate_reflections too slow")
        self.assertGreater(len(rows), 0, "expected accessible peaks for DyBCO c-axis")
        df = pd.DataFrame(rows)
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as handle:
            path = handle.name
        df.to_csv(path, index=False)
        self.assertGreater(os.path.getsize(path), 0)
        os.remove(path)

    def test_plan_reflection_005(self):
        plan = plan_reflection(
            (0, 0, 5),
            self.reciprocal,
            wavelength=self.wavelength,
            r_mount=self.r_mount,
            full_search=True,
        )
        self.assertTrue(plan.accessible)
        self.assertAlmostEqual(plan.two_theta_deg, 38.4, delta=0.5)

    def test_diffraction_and_ewald_projection(self):
        k_limit = auto_k_range(self.wavelength, 80.0)[1]
        grid = np.linspace(-k_limit, k_limit, 8)
        z_grid = np.linspace(0.0, k_limit, 8)
        Kx, Ky, Kz = np.meshgrid(grid, grid, z_grid)
        atoms = build_atom_dict(self.structure)
        worker = DiffractionPattern(
            Kx,
            Ky,
            Kz,
            self.xray,
            atoms,
            constant_atomic_structure_factor=False,
            useCUDA=False,
        )
        worker.run()
        intensity = worker.intensity
        geometry = XRDGeometry(
            wavelength_angstrom=self.wavelength,
            two_theta_deg=38.4,
            r_mount=self.r_mount,
        )
        weighted = geometry.apply_broadening(intensity, Kx, Ky, Kz)
        projection = np.sum(weighted, axis=0)
        self.assertEqual(projection.shape, (8, 8))
        xz = project_xrd_ewald_to_xz(
            intensity, Kx, Ky, Kz, self.wavelength, fwhm_deg=0.0
        )
        self.assertEqual(xz.shape, (8, 8))

    def test_four_circle_geometry(self):
        geometry = FourCircleGeometry(
            wavelength_angstrom=self.wavelength,
            omega_deg=0.0,
            chi_deg=0.0,
            phi_deg=0.0,
            two_theta_deg=38.4,
            r_mount=self.r_mount,
        )
        solutions = geometry.solve_bragg((0, 0, 5), self.reciprocal)
        self.assertTrue(solutions)
        ki, ks, q = geometry.forward_kin()
        self.assertAlmostEqual(np.linalg.norm(ki), geometry.k0(), places=4)

    def test_powder_crosscheck(self):
        calculator = XRDCalculator(wavelength=self.wavelength)
        pattern = calculator.get_pattern(
            self.structure, scaled=False, two_theta_range=(10, 80)
        )
        rows = accessible_peak_table(
            list(
                enumerate_reflections(
                    self.structure,
                    wavelength=self.wavelength,
                    h_max=4,
                    k_max=4,
                    l_max=8,
                    r_mount=self.r_mount,
                    two_theta_min=10.0,
                    two_theta_max=80.0,
                )
            )
        )
        planner_005 = next(row for row in rows if row["h"] == 0 and row["k"] == 0 and row["l"] == 5)
        powder_005 = [
            x
            for x, hkls in zip(pattern.x, pattern.hkls)
            if any(item["hkl"] == (0, 0, 5) for item in hkls)
        ]
        self.assertTrue(powder_005)
        self.assertAlmostEqual(planner_005["two_theta"], powder_005[0], delta=0.3)


if __name__ == "__main__":
    unittest.main()
