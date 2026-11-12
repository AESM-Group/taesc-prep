"""Checks for reference comparison and incomplete-dataset handling."""

import json
from pathlib import Path
import shutil
import tempfile
import unittest

from case_runner import build_suite, compare_stru, compare_text, load_cases


class TestCaseRunner(unittest.TestCase):
    def test_comparison_allows_whitespace_and_input_precision(self):
        with tempfile.TemporaryDirectory() as tmp:
            actual, reference = Path(tmp)/"actual", Path(tmp)/"reference"
            actual.write_text("1  2.00000001\n3.0D-4 row\n")
            reference.write_text("1\n2.0 0.0003   row")
            compare_text(actual, reference)
            with self.assertRaises(AssertionError):
                compare_text(actual, reference, rtol=1e-12, atol=1e-12)

    def test_comparison_rejects_wrong_values_and_nonfinite_data(self):
        pairs = [("1.0000001", "1.0"), ("1.01", "1.0"), ("2 row", "1 row"), ("1", "1 2"),
                 ("col", "row"), ("nan", "nan"), ("inf", "inf")]
        with tempfile.TemporaryDirectory() as tmp:
            actual, reference = Path(tmp)/"actual", Path(tmp)/"reference"
            for got, expected in pairs:
                with self.subTest(got=got, expected=expected):
                    actual.write_text(got)
                    reference.write_text(expected)
                    with self.assertRaises(AssertionError):
                        compare_text(actual, reference)
            actual.write_text("2")
            reference.write_text("1")
            with self.assertRaises(AssertionError):
                compare_text(actual, reference, atol=10)

    def test_declared_case_missing_data_errors(self):
        case = {"source": "abacus", "n_atoms": 2, "n_symops": 48}
        for missing in ("source", "reference"):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)
                (path / "cases.json").write_text(json.dumps({"001_incomplete": case}))
                if missing == "reference":
                    (path / "source" / "001_incomplete").mkdir(parents=True)
                result = unittest.TestResult()
                build_suite(path).run(result)
                self.assertEqual(result.testsRun, 1)
                self.assertEqual(len(result.errors), 1)
                self.assertIn(f"FileNotFoundError: Missing {missing}", result.errors[0][1])
                self.assertFalse(result.failures or result.skipped)

    def test_real_case_requires_complete_matching_reference(self):
        data = Path(__file__).parent / "data"
        name = "001_abacus_LiH_sym"
        case = load_cases(data / "cases.json")[name]
        for change in ("missing_output", "extra_output", "content", "empty_reference"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)
                (path / "source").symlink_to((data / "source").resolve(), target_is_directory=True)
                reference = path / "ref" / name
                shutil.copytree(data / "ref" / name, reference)
                (path / "cases.json").write_text(json.dumps({name: case}))
                if change == "missing_output":
                    shutil.copyfile(reference / "stru_out.txt", reference / "extra_reference")
                elif change == "extra_output":
                    (reference / "bz_sample.txt").unlink()
                elif change == "content":
                    bz = reference / "bz_sample.txt"
                    bz.write_text(bz.read_text() + "\nextra token\n")
                else:
                    for file in reference.iterdir():
                        file.unlink()
                result = unittest.TestResult()
                build_suite(path).run(result)
                self.assertEqual(result.testsRun, 1)
                self.assertFalse(result.skipped)
                if change == "empty_reference":
                    self.assertEqual(len(result.errors), 1)
                    self.assertIn("Missing reference files", result.errors[0][1])
                else:
                    self.assertEqual(len(result.failures), 1)
                    self.assertFalse(result.errors)
                    expected = {"missing_output": "missing=['extra_reference']",
                                "extra_output": "unexpected=['bz_sample.txt']",
                                "content": "token count"}[change]
                    self.assertIn(expected, result.failures[0][1])

    def test_invalid_manifest_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"cases.json"
            for manifest in ({}, {"bad": {"source": "abacus", "n_atoms": 2}},
                             {"bad": {"source": "abacus", "n_atoms": 2,
                                      "n_symops": 1, "atol": -1}},
                             {"bad": {"source": "abacus", "n_atoms": 2,
                                      "n_symops": 1, "skip": "Not ready"}}):
                path.write_text(json.dumps(manifest))
                with self.assertRaises(ValueError):
                    load_cases(path)


class TestStructureComparison(unittest.TestCase):
    """Exercise comparison using the real case's independently supplied reference."""

    def setUp(self):
        reference = Path(__file__).parent / "data/ref/001_abacus_LiH_sym/stru_out.txt"
        self.rows = reference.read_text().splitlines()
        self.symmetry_start = 7 + int(self.rows[6])
        self.convention = self.rows[self.symmetry_start].split()[1]

    def compare_rows(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            actual, reference = Path(tmp) / "actual", Path(tmp) / "reference"
            actual.write_text("\n".join(rows) + "\n")
            reference.write_text("\n".join(self.rows) + "\n")
            compare_stru(actual, reference)

    def test_symmetry_conventions_and_order(self):
        original = list(self.rows)
        for convention in ("row", "col"):
            rows = list(original)
            count = rows[self.symmetry_start].split()[0]
            rows[self.symmetry_start] = f"{count} {convention}"
            if convention != self.convention:
                for index in range(self.symmetry_start + 1, len(rows)):
                    tokens = rows[index].split()
                    tokens[:9] = [tokens[3 * j + i] for i in range(3) for j in range(3)]
                    rows[index] = " ".join(tokens)
            rows[self.symmetry_start + 1:] = reversed(rows[self.symmetry_start + 1:])
            for reference_rows in (original, rows):
                with self.subTest(actual=convention, reference=reference_rows is original):
                    self.rows = reference_rows
                    self.compare_rows(rows if reference_rows is original else original)
        self.rows = original

    def test_changed_or_duplicate_operation_fails(self):
        for change in ("rotation", "translation", "duplicate", "nonfinite"):
            with self.subTest(change=change):
                rows = list(self.rows)
                index = self.symmetry_start + 1
                tokens = rows[index].split()
                if change == "rotation":
                    tokens[0] = str(int(tokens[0]) + 1)
                elif change == "translation":
                    tokens[9] = str(float(tokens[9]) + 0.1)
                elif change == "nonfinite":
                    tokens[9] = "nan"
                else:
                    tokens = rows[index + 1].split()
                rows[index] = " ".join(tokens)
                with self.assertRaises(AssertionError):
                    self.compare_rows(rows)

    def test_missing_or_malformed_symmetry_fails(self):
        for change in ("absent", "truncated", "convention"):
            with self.subTest(change=change):
                rows = list(self.rows)
                if change == "absent":
                    rows = rows[:self.symmetry_start]
                elif change == "truncated":
                    rows.pop()
                else:
                    rows[self.symmetry_start] = "24 invalid"
                with self.assertRaises(AssertionError):
                    self.compare_rows(rows)

    def test_geometry_difference_still_fails(self):
        rows = list(self.rows)
        tokens = rows[7].split()
        tokens[0] = str(float(tokens[0]) + 0.1)
        rows[7] = " ".join(tokens)
        with self.assertRaises(AssertionError):
            self.compare_rows(rows)


if __name__ == "__main__":
    unittest.main()
