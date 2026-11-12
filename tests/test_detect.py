"""Calculation-source detection and CLI dispatch tests."""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from taesc_prep import AbacusSource, Source, detect_source
from taesc_prep import cli


class TestDetectSource(unittest.TestCase):
    def test_abacus_input_and_log_directories(self):
        for filename in ("STRU", "running_scf.log", "OUT.ABACUS/running_nscf.log"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                marker = directory / filename
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.touch()
                source = detect_source(tmp)
                self.assertIsInstance(source, AbacusSource)
                self.assertEqual(source.directory, directory)
                self.assertTrue(source.process_stru)
                self.assertFalse((directory / "stru_out").exists())

    def test_configured_structure_and_output_directory(self):
        for configured in ("custom.stru", "OUT.custom/running_nscf.log"):
            with self.subTest(configured=configured), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                (directory / "INPUT").write_text(
                    "INPUT_PARAMETERS\nstru_file custom.stru\nsuffix custom\ncalculation nscf\n")
                path = directory / configured
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
                self.assertIsInstance(detect_source(directory), AbacusSource)

    def test_ambiguous_logs_keep_discovery_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "running_scf.log").touch()
            (directory / "running_nscf.log").touch()
            with self.assertRaisesRegex(ValueError, "Multiple running logs found"):
                detect_source(directory)

    def test_aims_is_recognized_but_not_supported(self):
        for filename in ("geometry.in", "control.in", "aims.out"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                (directory / filename).touch()
                with self.assertRaisesRegex(ValueError, "FHI-aims directory detected"):
                    detect_source(directory)

    def test_mixed_sources_are_ambiguous(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "STRU").touch()
            (directory / "aims.out").touch()
            with self.assertRaisesRegex(ValueError, "Both ABACUS and FHI-aims"):
                detect_source(directory)

    def test_unknown_directory_and_invalid_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            # Directory names alone must not identify input files.
            (directory / "STRU").mkdir()
            with self.assertRaisesRegex(ValueError, "Cannot detect"):
                detect_source(directory)
            file_path = directory / "file"
            file_path.touch()
            for path in (file_path, directory / "missing"):
                with self.subTest(path=path):
                    with self.assertRaisesRegex(ValueError, "Expected a calculation directory"):
                        detect_source(path)

    def test_main_converts_the_detected_source(self):
        class OtherSource(Source):
            def __init__(self, directory):
                super().__init__(directory)
                self.process_stru = True

            def _extract_stru(self, *, include_symmetry=True):
                return {"lattice": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                        "positions": [[0, 0, 0]], "types": [1]}

        with tempfile.TemporaryDirectory() as tmp:
            source = OtherSource(tmp)
            with patch("taesc_prep.cli.detect_source", return_value=source) as detect:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(cli.main([tmp]), 0)
                detect.assert_called_once_with(Path(tmp))
            rows = (Path(tmp) / "stru_out").read_text().splitlines()
            self.assertEqual(rows[6], "1")
            self.assertEqual(list(map(float, rows[7].split())), [0, 0, 0, 1])


if __name__ == "__main__":
    unittest.main()
