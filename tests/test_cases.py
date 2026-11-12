"""Regression tests using real datasets."""

from fnmatch import fnmatchcase
from pathlib import Path
import unittest

import pytest

from case_runner import build_suite


DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.parametrize("case", list(build_suite(DATA_DIR)), ids=lambda case: case.name)
def test_conversion(case):
    case.runTest()


def load_tests(loader, tests, pattern):
    suite = build_suite(DATA_DIR)
    if loader.testNamePatterns:
        suite = unittest.TestSuite(test for test in suite
                                   if any(fnmatchcase(test.id(), pattern)
                                          for pattern in loader.testNamePatterns))
    return suite


if __name__ == "__main__":
    unittest.main()
