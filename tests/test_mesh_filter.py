"""The sampled-mesh filter that symrot_R.txt needs before it drives k-stars."""
import unittest

from taesc_prep.abacus import preserves_mesh, reciprocal_rotation


class MeshFilterTests(unittest.TestCase):
    def test_identity_preserves_every_mesh(self):
        identity = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        for mesh in ((1, 1, 3), (1, 1, 2), (2, 2, 2), (2, 3, 4)):
            self.assertTrue(preserves_mesh(identity, mesh), mesh)

    def test_cubic_only_operation(self):
        # maps a cubic mesh onto itself, but not (1,1,3)
        operation = [1, 0, -1, 1, 0, 0, 1, -1, 0]
        self.assertTrue(preserves_mesh(operation, (2, 2, 2)))
        self.assertFalse(preserves_mesh(operation, (1, 1, 3)))
        self.assertFalse(preserves_mesh(operation, (2, 2, 1)))

    def test_rotation_about_z_survives_a_z_mesh(self):
        rotation = [0, -1, 0, 1, 0, 0, 0, 0, 1]
        self.assertTrue(preserves_mesh(rotation, (1, 1, 3)))
        self.assertTrue(preserves_mesh(rotation, (2, 2, 2)))

    def test_antunitary_flips_the_mesh_test(self):
        # Theta * inversion reverses k, so a mesh is preserved only if it is
        # symmetric about the origin.
        inversion = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
        self.assertTrue(preserves_mesh(inversion, (2, 2, 2)))
        self.assertTrue(preserves_mesh(inversion, (2, 2, 2), antiunitary=True))

    def test_reciprocal_of_a_known_operation(self):
        # [1,0,-1, 1,0,0, 1,-1,0] has det = -1; its reciprocal is
        # [0,0,-1, 1,1,1, 0,-1,0] (verified against the ABACUS-side helper).
        self.assertEqual(reciprocal_rotation([1, 0, -1, 1, 0, 0, 1, -1, 0]),
                         [0, 0, -1, 1, 1, 1, 0, -1, 0])

    def test_non_unimodular_rotation_is_rejected(self):
        with self.assertRaises(ValueError):
            reciprocal_rotation([2, 0, 0, 0, 1, 0, 0, 0, 1])


if __name__ == "__main__":
    unittest.main()


class ShiftedGridTests(unittest.TestCase):
    """A shifted grid samples points that a zero-offset test accepts wrongly."""

    INVERSION = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
    IDENTITY = [1, 0, 0, 0, 1, 0, 0, 0, 1]

    def test_gamma_grid_keeps_inversion(self):
        self.assertTrue(preserves_mesh(self.INVERSION, (2, 1, 1), (0.0, 0.0, 0.0)))

    def test_quarter_shifted_grid_drops_inversion(self):
        # Points {1/8, 5/8} map to {7/8, 3/8}, outside the sampled set.
        self.assertFalse(preserves_mesh(self.INVERSION, (2, 1, 1), (0.25, 0.0, 0.0)))

    def test_quarter_shifted_grid_keeps_identity(self):
        self.assertTrue(preserves_mesh(self.IDENTITY, (2, 1, 1), (0.25, 0.0, 0.0)))

    def test_half_shifted_grid_keeps_inversion(self):
        # Points {1/4, 3/4} are symmetric about the origin.
        self.assertTrue(preserves_mesh(self.INVERSION, (2, 1, 1), (0.5, 0.0, 0.0)))
