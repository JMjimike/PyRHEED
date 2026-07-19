"""Validation tests for XRD extensions to PyRHEED."""

from __future__ import annotations

import math
import os
import sys
import unittest

import numpy as np

SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC_DIR))

from electron_form_factors import PengFormFactorTable
from four_circle_geometry import (
    bragg_error,
    forward_kinematics,
    solve_bragg_angles,
)
from process import DiffractionPattern
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from scattering_factors import XRayFormFactorTable
from structure_io import load_ordered_structure, orient_crystal_axes
from xrd_reflection_planner import plan_reflection


CIF_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "cif_samples",
        "orthorhombic_O7",
        "DyBa2Cu3O7_mp-622105.cif",
    )
)


class TestXRDExtensions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = orient_crystal_axes(
            load_ordered_structure(CIF_PATH), "c", "a"
        )
        cls.reciprocal = cls.structure.lattice.reciprocal_lattice.matrix
        cls.wavelength = 1.54184

    def test_xray_form_factor_table_loads(self):
        table = XRayFormFactorTable.load()
        k = np.array([0.0, 1.0, 2.0])
        f = table.fe(k, k, k, "Cu")
        self.assertEqual(f.shape, k.shape)
        self.assertTrue(np.all(f > 0))

    def test_diffraction_pattern_accepts_form_factor_table(self):
        atoms = {(0.0, 0.0, 0.0): "Cu"}
        k = np.array([[[0.0]]])
        peng = PengFormFactorTable.load(os.path.join(SRC_DIR, "files"))
        xray = XRayFormFactorTable.load()
        worker_peng = DiffractionPattern(k, k, k, peng, atoms, True, False)
        worker_xray = DiffractionPattern(
            k, k, k, form_factor_table=xray, atoms=atoms, constant_atomic_structure_factor=True, useCUDA=False
        )
        self.assertIs(worker_peng.form_factor_table, peng)
        self.assertIs(worker_xray.form_factor_table, xray)

    def test_dybco_005_two_theta_cu_kalpha(self):
        plan = plan_reflection(
            (0, 0, 5),
            self.reciprocal,
            wavelength=self.wavelength,
        )
        self.assertTrue(plan.accessible)
        self.assertAlmostEqual(plan.two_theta_deg, 38.4, delta=0.5)
        self.assertAlmostEqual(plan.omega_deg, 0.0, delta=0.1)
        self.assertAlmostEqual(plan.chi_deg, 0.0, delta=0.1)
        self.assertAlmostEqual(plan.phi_deg, 0.0, delta=0.1)
        self.assertEqual(plan.scan_mode, "theta_2theta")

    def test_xrd_calculator_agrees_with_planner_00l(self):
        calculator = XRDCalculator(wavelength=self.wavelength)
        pattern = calculator.get_pattern(self.structure, scaled=False, two_theta_range=(10, 90))
        plan = plan_reflection((0, 0, 5), self.reciprocal, wavelength=self.wavelength)
        matched = [
            two_theta
            for two_theta, hkls in zip(pattern.x, pattern.hkls)
            if any(item["hkl"] == (0, 0, 5) for item in hkls)
        ]
        self.assertTrue(matched)
        self.assertAlmostEqual(matched[0], plan.two_theta_deg, delta=0.2)

    def test_four_circle_inverse_forward_bragg(self):
        solutions = solve_bragg_angles(
            (0, 0, 5),
            self.reciprocal,
            self.wavelength,
        )
        self.assertTrue(solutions)
        angles = solutions[0]
        g_crystal = np.array([0, 0, 5]) @ self.reciprocal
        err = bragg_error(
            angles.omega,
            angles.chi,
            angles.phi,
            angles.two_theta,
            g_crystal,
            self.wavelength,
        )
        self.assertLess(err, 0.2)
        _ki, _ks, q_lab = forward_kinematics(
            angles.omega,
            angles.chi,
            angles.phi,
            angles.two_theta,
            self.wavelength,
        )
        self.assertGreater(np.linalg.norm(q_lab), 0.0)


if __name__ == "__main__":
    unittest.main()
