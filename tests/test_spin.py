"""Cross-checks for the spinor representations of spatial operations."""
import math
import unittest

from taesc_prep.spin import cartesian_rotation, flatten_su2, proper_part, rotation_to_su2


class SpinRepresentationTests(unittest.TestCase):
    def test_identity_is_the_unit_matrix(self):
        u = rotation_to_su2([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual(flatten_su2(u),
                         [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])

    def test_rotations_about_z_match_the_reference(self):
        # ABACUS ModuleSymmetry::SpinRotation::so3_to_su2 values.
        cases = {
            90: [0.7071067811865476, 0.7071067811865475, 0.0, 0.0,
                 0.0, 0.0, 0.7071067811865476, -0.7071067811865475],
            180: [6.123233995736766e-17, -1.0, 0.0, 0.0,
                  0.0, 0.0, 6.123233995736766e-17, 1.0],
            60: [0.8660254037844387, 0.5, 0.0, 0.0,
                 0.0, 0.0, 0.8660254037844387, -0.4999999999999999],
        }
        c = math.cos(math.radians(1))
        for degrees, expected in cases.items():
            angle = math.radians(degrees)
            rotation = [[math.cos(angle), -math.sin(angle), 0.0],
                        [math.sin(angle), math.cos(angle), 0.0],
                        [0.0, 0.0, 1.0]]
            got = flatten_su2(rotation_to_su2(rotation))
            for a, b in zip(got, expected):
                self.assertAlmostEqual(a, b, places=12)

    def test_improper_operations_use_the_proper_part(self):
        # -I has det = -1; flipping it back gives identity in SU(2).
        u = rotation_to_su2([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        self.assertEqual(flatten_su2(u),
                         [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        # The whole matrix is negated, as ABACUS's proper_part does.
        self.assertEqual(proper_part([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                         [[1, 0, 0], [0, -1, 0], [0, 0, -1]])

    def test_theta_pi_uses_the_symmetric_axis(self):
        # 180 degrees about x: R = diag(1, -1, -1).
        u = flatten_su2(rotation_to_su2([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
        self.assertAlmostEqual(u[0], 6.123233995736766e-17, places=15)
        self.assertAlmostEqual(u[3], -1.0, places=12)
        self.assertAlmostEqual(u[5], -1.0, places=12)

    def test_fractional_to_cartesian_uses_the_lattice(self):
        # latvec^-1 R latvec, the conversion ABACUS applies before so3_to_su2.
        lattice = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
        fractional = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        self.assertEqual(cartesian_rotation(lattice, fractional),
                         [[0.0, -1.5, 0.0], [2 / 3, 0.0, 0.0], [0.0, 0.0, 1.0]])

    def test_cubic_lattice_passes_the_rotation_through(self):
        identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        fractional = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        self.assertEqual(cartesian_rotation(identity, fractional), fractional)


if __name__ == "__main__":
    unittest.main()
