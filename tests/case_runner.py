"""Run conversions against independently prepared real-calculation datasets."""

import json
import math
from pathlib import Path
import re
import tempfile
import unittest

from taesc_prep import AbacusSource, detect_source


SOURCE_CLASSES = {"abacus": AbacusSource}
DEFAULT_TOLERANCE = 1e-8


def load_cases(path):
    """Read the case manifest; every declared case must run."""
    cases = json.loads(Path(path).read_text())
    if not isinstance(cases, dict) or not cases:
        raise ValueError("cases.json must contain a nonempty object of named cases")
    for name, case in cases.items():
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or not isinstance(case, dict):
            raise ValueError(f"Invalid case entry: {name!r}")
        if not isinstance(case.get("source"), str) or not case["source"]:
            raise ValueError(f"{name}: source must be a nonempty string")
        for key in ("n_atoms", "n_symops"):
            if type(case.get(key)) is not int or case[key] < (1 if key == "n_atoms" else 0):
                raise ValueError(f"{name}: invalid {key}")
        if "skip" in case:
            raise ValueError(f"{name}: skip is not supported; every declared case must run")
        files = case.get("files", ["stru_out"])
        if not isinstance(files, list) or not files or len(files) != len(set(files)):
            raise ValueError(f"{name}: files must be a nonempty list without duplicates")
        for filename in files:
            if (not isinstance(filename, str) or not filename or filename == "."
                    or Path(filename).is_absolute() or ".." in Path(filename).parts):
                raise ValueError(f"{name}: invalid output filename {filename!r}")
        for key in ("rtol", "atol"):
            value = case.get(key, DEFAULT_TOLERANCE)
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name}: {key} must be finite and nonnegative")
    return cases


def compare_text(actual, reference, *, rtol=DEFAULT_TOLERANCE, atol=DEFAULT_TOLERANCE):
    """Compare text tokens, allowing only the requested numerical differences.

    Whitespace is ignored. Integer pairs and text are compared exactly.
    Other numerical pairs must be finite and satisfy math.isclose.
    """
    actual, reference = Path(actual), Path(reference)
    actual_tokens = actual.read_text().split()
    reference_tokens = reference.read_text().split()
    _compare_tokens(actual_tokens, reference_tokens, reference.name, rtol=rtol, atol=atol)


def _compare_tokens(actual_tokens, reference_tokens, label, *, rtol, atol):
    if len(actual_tokens) != len(reference_tokens):
        raise AssertionError(f"{label}: token count: got {len(actual_tokens)}, "
                             f"expected {len(reference_tokens)}")
    for index, (got, expected) in enumerate(zip(actual_tokens, reference_tokens), 1):
        context = f"{label}, token {index}: got {got!r}, expected {expected!r}"
        if re.fullmatch(r"[+-]?\d+", got) and re.fullmatch(r"[+-]?\d+", expected):
            if int(got) != int(expected):
                raise AssertionError(context)
            continue
        try:
            a = float(got.replace("D", "E").replace("d", "e"))
            b = float(expected.replace("D", "E").replace("d", "e"))
        except ValueError:
            if got != expected:
                raise AssertionError(context)
            continue
        if not math.isfinite(a) or not math.isfinite(b):
            raise AssertionError(context + "; non-finite number")
        if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
            raise AssertionError(context + f"; rtol={rtol}, atol={atol}")



def _read_stru(path):
    """Read geometry tokens and symmetry operations normalized to row convention."""
    path = Path(path)
    tokens = path.read_text().split()
    try:
        n_atoms = int(tokens[18])
        geometry_end = 19 + 4 * n_atoms
        if n_atoms < 1 or len(tokens) < geometry_end:
            raise ValueError("incomplete geometry")
        geometry = tokens[:geometry_end]
        tail = tokens[geometry_end:]
        if not tail:
            return geometry, []
        n_operations = int(tail[0])
        convention = tail[1]
        if convention not in ("row", "col"):
            raise ValueError(f"unknown symmetry convention {convention!r}")
        if n_operations < 0 or len(tail) != 2 + 12 * n_operations:
            raise ValueError("symmetry operation count does not match block length")
        operations = []
        for start in range(2, len(tail), 12):
            rotation = tuple(int(x) for x in tail[start:start + 9])
            if convention == "col":
                rotation = tuple(rotation[3 * j + i] for i in range(3) for j in range(3))
            translation = tuple(float(x.replace("D", "E").replace("d", "e"))
                                for x in tail[start + 9:start + 12])
            if not all(math.isfinite(x) for x in translation):
                raise ValueError("non-finite symmetry translation")
            operations.append((rotation, translation))
        return geometry, operations
    except (IndexError, ValueError) as error:
        raise AssertionError(f"{path}: invalid stru_out: {error}") from error


def _compare_symmetry(actual, reference, label, *, rtol, atol):
    """Match normalized operations one-to-one, retaining duplicate counts."""
    if len(actual) != len(reference):
        raise AssertionError(f"{label}: symmetry count: got {len(actual)}, "
                             f"expected {len(reference)}")
    candidates = [[j for j, (rotation, translation) in enumerate(reference)
                   if got_rotation == rotation
                   and all(math.isclose(a, b, rel_tol=rtol, abs_tol=atol)
                           for a, b in zip(got_translation, translation))]
                  for got_rotation, got_translation in actual]
    # Reassign earlier matches if translations have overlapping tolerance ranges.
    matched = {}

    def match(index, visited):
        for candidate in candidates[index]:
            if candidate in visited:
                continue
            visited.add(candidate)
            if candidate not in matched or match(matched[candidate], visited):
                matched[candidate] = index
                return True
        return False

    for index, operation in enumerate(actual):
        if not match(index, set()):
            raise AssertionError(f"{label}: unmatched symmetry operation {index + 1} "
                                 f"in row convention: {operation}; rtol={rtol}, atol={atol}")


def compare_stru(actual, reference, *, rtol=DEFAULT_TOLERANCE, atol=DEFAULT_TOLERANCE):
    """Compare geometry in order and spatial symmetry independently of convention/order."""
    actual_geometry, actual_operations = _read_stru(actual)
    reference_geometry, reference_operations = _read_stru(reference)
    label = Path(reference).name
    _compare_tokens(actual_geometry, reference_geometry, label, rtol=rtol, atol=atol)
    _compare_symmetry(actual_operations, reference_operations, label, rtol=rtol, atol=atol)


class ConversionCase(unittest.TestCase):
    """One manifest case, with source and reference data left unchanged."""

    def __init__(self, data_dir, name, case):
        super().__init__("runTest")
        self.data_dir = Path(data_dir)
        self.name = name
        self.case = case

    def id(self):
        return f"test_cases.{self.name}"

    def __str__(self):
        return self.id()

    def shortDescription(self):
        return None

    def runTest(self):
        case = self.case
        source_dir = self.data_dir / "source" / self.name
        reference_dir = self.data_dir / "ref" / self.name
        files = case.get("files", ["stru_out"])
        if not source_dir.is_dir():
            raise FileNotFoundError(f"Missing source directory: {source_dir}")
        for filename in files:
            if not (reference_dir / filename).is_file():
                raise FileNotFoundError(f"Missing reference file: {reference_dir / filename}")
        self.assertIn(case["source"], SOURCE_CLASSES,
                      f"Source implementation not registered for {case['source']!r}")
        source = detect_source(source_dir)
        self.assertIsInstance(source, SOURCE_CLASSES[case["source"]])
        with tempfile.TemporaryDirectory(prefix=f"taesc-{self.name}-") as tmp:
            output_dir = Path(tmp) / "output"
            result = source.convert(output_dir)
            self.assertEqual(result.output, output_dir)
            self.assertEqual(result.n_atoms, case["n_atoms"], "Atom count")
            self.assertEqual(result.n_symops, case["n_symops"], "Symmetry-operation count")
            generated = {p.relative_to(output_dir).as_posix()
                         for p in output_dir.rglob("*") if p.is_file()}
            self.assertEqual(generated, set(files), "Generated output files")
            for filename in files:
                with self.subTest(file=filename):
                    compare = compare_stru if filename == "stru_out" else compare_text
                    compare(output_dir / filename, reference_dir / filename,
                                 rtol=case.get("rtol", DEFAULT_TOLERANCE),
                                 atol=case.get("atol", DEFAULT_TOLERANCE))


def build_suite(data_dir):
    """Build one unittest entry per case; manifest errors are test errors."""
    data_dir = Path(data_dir)
    cases = load_cases(data_dir / "cases.json")
    return unittest.TestSuite(ConversionCase(data_dir, name, case)
                              for name, case in sorted(cases.items()))
