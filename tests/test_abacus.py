"""Analytic geometry tests and ABACUS-format parser/CLI tests."""

import contextlib
import io
import math
from pathlib import Path
import tempfile
import unittest

from taesc_prep import AbacusSource, ConversionResult, Source, abacus, cli, structure
from taesc_prep.output import write_taesc_stru


STRU = """ATOMIC_SPECIES
Li 6.94 Li.upf
H 1.008 H.upf
NUMERICAL_ORBITAL
Li.orb
H.orb
LATTICE_CONSTANT
2.0D0 // Bohr
LATTICE_VECTORS
1 0 0
0.5 1 0
0 0 2
ATOMIC_POSITIONS
Direct
Li
0.0
1
0 0 0 0 0 0
H
0.0
1
-0.25 0.5 0.25 1 1 1 mag 0.0 # trailing properties
"""

# The labels/table layout reproduce ABACUS 3.9 running_scf.log. Geometry is
# synthetic so expected Cartesian coordinates can be checked analytically.
LOG = """READING UNITCELL INFORMATION
Lattice constant (Bohr) = 2
READING ATOM TYPE 1
Atom label = Li
Number of atoms for this type = 1
READING ATOM TYPE 2
Atom label = H
Number of atoms for this type = 1
TOTAL ATOM NUMBER = 2
DIRECT COORDINATES
atom x y z mag
Li 0 0 0 0
H 0.75 0.5 0.25 0
Lattice vectors: (Cartesian coordinate: in unit of a_0)
1 0 0
0.5 1 0
0 0 2
Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)
1 -0.5 0
0 1 0
0 0 0.5
LCAO ALGORITHM ---- ION= 1 ELEC= 1
"""

# symm_rot_out.cpp prints the spatial operation before the orbital matrices.
# The second operation tests a nonzero fractional translation and a matrix
# whose transpose differs. The parser must copy it without transposition.
SYMROT = """Number of irreducible sector: 4
Lmax of AOs: 0
Lmax of ABFs: 0
Format:
The index of the symmetry operation
The rotation matrix of this symmetry operation (3*3)
(The translation vector of this symmetry operation)
Orbital rotation matrix (T) of each angular momentum with size ((2l + 1) * (2l + 1))

0
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
(0.000000 0.000000 0.000000)
(1,0)
1
0.000000 -1.000000 0.000000
1.000000 -1.000000 0.000000
0.000000 0.000000 1.000000
(0.250000 0.500000 0.750000)
(1,0)
"""


def render_structure(geometry):
    with tempfile.TemporaryDirectory() as tmp:
        path = write_taesc_stru(geometry, Path(tmp) / "stru_out")
        return path.read_text()


class TestGeometry(unittest.TestCase):
    def assertRowsAlmostEqual(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for row, ref in zip(actual, expected):
            for value, target in zip(row, ref):
                self.assertAlmostEqual(value, target, places=12)

    def test_skew_lattice_reciprocal_and_wrapping(self):
        geometry = abacus.read_stru(STRU)
        self.assertEqual(geometry["types"], [1, 2])
        self.assertRowsAlmostEqual(geometry["positions"], [[0, 0, 0], [2, 1, 1]])
        b = structure.reciprocal_vectors(geometry["lattice"])
        self.assertRowsAlmostEqual(b, [[math.pi, -math.pi/2, 0],
                                      [0, math.pi, 0], [0, 0, math.pi/2]])
        for i, arow in enumerate(geometry["lattice"]):
            for j, brow in enumerate(b):
                self.assertAlmostEqual(sum(x*y for x, y in zip(arow, brow)),
                                       2*math.pi if i == j else 0)
        out = render_structure(geometry).splitlines()
        self.assertEqual(len(out), 9)  # No legacy k-point block.
        self.assertEqual(out[6], "2")
        self.assertEqual(list(map(float, out[8].split())), [2, 1, 1, 2])

    def test_cartesian_coordinate_units(self):
        for mode, xyz in (("Cartesian", [1, 0.5, 0.5]),
                          ("Cartesian_au", [2, 1, 1]),
                          ("Cartesian_angstrom", [1.058354, 0.529177, 0.529177])):
            with self.subTest(mode=mode):
                text = STRU.replace("\nDirect\n", f"\n{mode}\n").replace(
                    "-0.25 0.5 0.25", " ".join(map(str, xyz)))
                self.assertRowsAlmostEqual(abacus.read_stru(text)["positions"],
                                           [[0, 0, 0], [2, 1, 1]])

    def test_cartesian_center_uses_cartesian_components(self):
        text = STRU.replace("\nDirect\n", "\nCartesian_angstrom_center_xy\n")
        # Full Cartesian center is (1.5, 1, 2), with z suppressed for xy.
        self.assertRowsAlmostEqual([abacus.read_stru(text)["positions"][0]], [[1.5, 1, 0]])

    def test_log_direct_and_old_cartesian_labels(self):
        self.assertEqual(abacus.read_log(LOG), abacus.read_stru(STRU))
        old = LOG.replace("DIRECT COORDINATES", "CARTESIAN COORDINATES ( UNIT = 2 Bohr ).")
        old = old.replace("Li 0 0 0 0", "tauc_Li1 0 0 0 0")
        old = old.replace("H 0.75 0.5 0.25 0", "tauc_H1 1 0.5 0.5 0")
        old = old.replace("Number of atoms", "number of atom").replace("Atom label", "atom label")
        self.assertEqual(abacus.read_log(old), abacus.read_stru(STRU))

    def test_reject_incomplete_or_ambiguous_data(self):
        cases = [(abacus.read_stru, STRU.replace("-0.25 0.5 0.25", "nan 0 0")),
                 (abacus.read_stru, STRU.replace("0 0 2", "0 0 0")),
                 (abacus.read_stru, STRU.replace("\nH\n0.0", "\nHe\n0.0")),
                 (abacus.read_stru, STRU.rsplit("\n", 2)[0]),
                 (abacus.read_log, LOG + "\nLCAO ALGORITHM ---- ION= 2 ELEC= 1"),
                 (abacus.read_log, LOG.replace("TOTAL ATOM NUMBER = 2", "TOTAL ATOM NUMBER = 3")),
                 (abacus.read_log, LOG + "\nDIRECT COORDINATES\n")]
        for reader, text in cases:
            with self.subTest(reader=reader.__name__, text=text[-60:]):
                with self.assertRaises(ValueError):
                    reader(text)

    def test_symmetry_convention_and_translation(self):
        ops = abacus.read_symrot(SYMROT)
        self.assertEqual(len(ops), 2)
        self.assertEqual(ops[1], ([0, -1, 0, 1, -1, 0, 0, 0, 1], [0.25, 0.5, 0.75]))
        geometry = abacus.read_stru(STRU)
        geometry["operations"] = ops
        out = render_structure(geometry).splitlines()
        self.assertEqual(out[9], "2 row")
        self.assertEqual(list(map(float, out[11].split())),
                         [0, -1, 0, 1, -1, 0, 0, 0, 1, 0.25, 0.5, 0.75])
        for bad in (SYMROT.rsplit("\n", 2)[0], SYMROT.replace("0.250000", "nan"),
                    SYMROT.replace("1.000000 0.000000 0.000000", "1.100000 0.000000 0.000000")):
            with self.assertRaises(ValueError):
                abacus.read_symrot(bad)

    def test_symmetry_skips_all_orbital_matrix_rows(self):
        text = SYMROT.replace("Lmax of ABFs: 0", "Lmax of ABFs: 1")
        text = text.replace("(1,0)\n", "(1,0)\n(1,0)(0,0)(0,0)\n"
                            "(0,0)(1,0)(0,0)\n(0,0)(0,0)(1,0)\n")
        self.assertEqual(abacus.read_symrot(text), abacus.read_symrot(SYMROT))


class TestSource(unittest.TestCase):
    def test_source_requires_structure_extraction(self):
        with self.assertRaises(TypeError):
            Source(".")

        class IncompleteSource(Source):
            pass

        with self.assertRaises(TypeError):
            IncompleteSource(".")

    def test_abacus_convert_writes_structure_and_returns_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").write_text(STRU)
            (directory / "symrot_R.txt").write_text(SYMROT)
            source = AbacusSource(tmp)
            self.assertIsInstance(source, Source)
            stdout, stderr = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                result = source.convert()
            self.assertEqual(result.output, directory)
            self.assertEqual((result.n_atoms, result.n_symops, result.notes), (2, 2, ()))
            rows = (result.output / "stru_out").read_text().splitlines()
            self.assertEqual(list(map(float, rows[8].split())), [2, 1, 1, 2])
            self.assertEqual(rows[9], "2 row")
            self.assertEqual(stdout.getvalue(), "")
            self.assertEqual(stderr.getvalue(), "")
            custom = source.convert(str(directory / "custom-output"), include_symmetry=False)
            self.assertEqual(custom.n_symops, 0)
            self.assertEqual(len((custom.output / "stru_out").read_text().splitlines()), 9)

    def test_shared_conversion_with_another_source(self):
        class FixedSource(Source):
            def __init__(self, directory):
                super().__init__(directory)
                self.process_stru = True

            def _extract_stru(self, *, include_symmetry=True):
                return {"lattice": [[2, 0, 0], [0, 3, 0], [0, 0, 4]],
                        "positions": [[1, 1.5, 2]], "types": [1]}

        with tempfile.TemporaryDirectory() as tmp:
            source = FixedSource(tmp)
            result = source.convert()
            rows = (result.output / "stru_out").read_text().splitlines()
            self.assertEqual((result.n_atoms, result.n_symops, result.notes), (1, 0, ()))
            self.assertEqual(list(map(float, rows[3].split())), [math.pi, 0, 0])
            self.assertEqual(list(map(float, rows[7].split())), [1, 1.5, 2, 1])
            converted = source.convert(Path(tmp) / "other")
            self.assertEqual((converted.output / "stru_out").read_text(), (result.output / "stru_out").read_text())

    def test_disabled_structure_conversion_does_not_read_or_write(self):
        class DisabledSource(Source):
            def _extract_stru(self, *, include_symmetry=True):
                raise AssertionError("Disabled extraction must not run")

        with tempfile.TemporaryDirectory() as tmp:
            source = DisabledSource(tmp)
            result = source.convert()
            self.assertIsInstance(result, ConversionResult)
            self.assertEqual((result.output, result.n_atoms, result.n_symops, result.notes),
                             (None, 0, 0, ()))
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_structure_helper_returns_dict_and_convert_returns_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").write_text(STRU)
            source = AbacusSource(tmp)
            extracted = source.convert_stru()
            self.assertIsInstance(extracted, dict)
            self.assertEqual(extracted["types"], [1, 2])
            self.assertTrue((directory / "stru_out").is_file())
            result = source.convert()
            self.assertIsInstance(result, ConversionResult)
            self.assertEqual(result.n_atoms, 2)
            source.process_stru = False
            disabled = source.convert(directory / "unused")
            self.assertEqual((disabled.output, disabled.n_atoms, disabled.n_symops, disabled.notes),
                             (None, 0, 0, ()))
            self.assertFalse((directory / "unused").exists())

    def test_extraction_returns_dict_with_optional_keys(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").write_text(STRU)
            source = AbacusSource(directory)
            extracted = source._extract_stru()
            self.assertIsInstance(extracted, dict)
            self.assertEqual(extracted["types"], [1, 2])
            self.assertEqual(extracted["positions"], [[0, 0, 0], [2, 1, 1]])
            self.assertNotIn("operations", extracted)
            self.assertNotIn("notes", extracted)
            (directory / "symrot_R.txt").write_text(SYMROT)
            extracted = source._extract_stru()
            self.assertEqual(len(extracted["operations"]), 2)
            self.assertNotIn("operations", source._extract_stru(include_symmetry=False))

    def test_missing_required_keys_preserves_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "stru_out"
            for key in ("lattice", "positions", "types"):
                with self.subTest(key=key):
                    geometry = abacus.read_stru(STRU)
                    del geometry[key]
                    output.write_text("existing output\n")
                    with self.assertRaisesRegex(ValueError, "Missing required structure keys: " + key):
                        write_taesc_stru(geometry, output)
                    self.assertEqual(output.read_text(), "existing output\n")

    def test_conversion_does_not_overwrite_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            geometry_path = directory / "STRU"
            symmetry_path = directory / "symrot_R.txt"
            geometry_path.write_text(STRU)
            symmetry_path.write_text(SYMROT)
            source = AbacusSource(directory)
            for path, contents in ((geometry_path, STRU), (symmetry_path, SYMROT)):
                with self.subTest(path=path):
                    output_dir = directory / (path.name + "-output")
                    output_dir.mkdir()
                    (output_dir / "stru_out").symlink_to(path)
                    with self.assertRaisesRegex(ValueError, "Input and output paths must differ"):
                        source.convert(output_dir)
                    self.assertEqual(path.read_text(), contents)

    def test_log_convert_returns_precision_note(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "running_scf.log").write_text(LOG)
            result = AbacusSource(directory).convert()
            self.assertEqual(result.n_atoms, 2)
            self.assertEqual(len(result.notes), 1)
            self.assertIn("precision", result.notes[0])


class TestCommandLine(unittest.TestCase):
    def test_directory_discovery_custom_names_and_symmetry(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "INPUT").write_text("INPUT_PARAMETERS\nstru_file custom.stru\nsuffix test\n")
            (directory / "custom.stru").write_text(STRU)
            outdir = directory / "OUT.test"
            outdir.mkdir()
            (outdir / "symrot_R.txt").write_text(SYMROT)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main([tmp]), 0)
            self.assertIn("2 row\n", (directory / "stru_out").read_text())
            self.assertEqual((directory / "custom.stru").read_text(), STRU)
            with contextlib.redirect_stdout(io.StringIO()):
                cli.main([tmp, "--no-symmetry"])
            self.assertEqual(len((directory / "stru_out").read_text().splitlines()), 9)

    def test_log_only_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "running_scf.log").write_text(LOG)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main([tmp]), 0)
            self.assertEqual((directory / "stru_out").read_text(),
                             render_structure(abacus.read_stru(STRU)))

    def test_configured_log_selection_and_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "INPUT").write_text("INPUT_PARAMETERS\nsuffix test\ncalculation nscf\n")
            outdir = directory / "OUT.test"
            outdir.mkdir()
            (outdir / "running_scf.log").write_text("unused older log")
            (outdir / "running_nscf.log").write_text(LOG)
            output = directory / "new" / "dataset"
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main([tmp, "--output", str(output)]), 0)
            self.assertTrue(output.is_dir())
            self.assertFalse((directory / "stru_out").exists())
            self.assertEqual((output / "stru_out").read_text(),
                             render_structure(abacus.read_stru(STRU)))

    def test_output_argument_rejects_an_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").write_text(STRU)
            output = directory / "existing-file"
            output.write_text("keep this content\n")
            errors = io.StringIO()
            with contextlib.redirect_stderr(errors), self.assertRaises(SystemExit) as error:
                cli.main([tmp, "--output", str(output)])
            self.assertEqual(error.exception.code, 2)
            self.assertIn("Output must be a directory", errors.getvalue())
            self.assertEqual(output.read_text(), "keep this content\n")

    def test_failure_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").write_text(STRU.replace("0 0 2", "0 0 0"))
            output = directory / "stru_out"
            output.write_text("existing output\n")
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                cli.main([tmp])
            self.assertEqual(error.exception.code, 2)
            self.assertEqual(output.read_text(), "existing output\n")

    def test_positional_argument_must_be_a_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "STRU"
            source.write_text(STRU)
            errors = io.StringIO()
            with contextlib.redirect_stderr(errors), self.assertRaises(SystemExit):
                cli.main([str(source)])
            self.assertIn("Expected a calculation directory", errors.getvalue())


if __name__ == "__main__":
    unittest.main()
