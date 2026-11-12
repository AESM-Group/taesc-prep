# Real-calculation regression tests

Each case converts an existing calculation directory and compares the generated
files with independently prepared references. The tests do not run ABACUS or
FHI-aims, and do not generate or update reference files.

## Layout

```text
tests/
  case_runner.py
  test_cases.py
  test_case_runner.py
  data/
    cases.json
    source/<case-name>/
    ref/<case-name>/
```

`source/<case-name>` contains the real calculation inputs and outputs needed
by the converter. Keep the producing software version and calculation settings
with these files. `ref/<case-name>` contains only the expected TAESC output files,
using the names `stru_out.txt` and `bz_sample.txt`.

The runner calls `detect_source` and `convert`, checks the detected source,
atom and symmetry-operation counts, and compares every file in the reference directory.
Each case has its own temporary output directory, which is removed after the
test. Source and reference data are not modified.

## Case manifest

Every case listed in `cases.json` is run:

```json
{
  "001_abacus_LiH_sym": {
    "source": "abacus",
    "n_atoms": 2,
    "n_symops": 24,
    "n_kpoints": 10,
    "n_full_kpoints": 125
  }
}
```

The counts are expectations supplied by the case author and checked against
the conversion result.

| Key | Required | Meaning |
| --- | --- | --- |
| `source` | Yes | Expected source implementation; currently `abacus`. |
| `n_atoms` | Yes | Expected atom count. |
| `n_symops` | Yes | Expected number of exported spatial symmetry operations. |
| `n_kpoints` | No | Expected SCF k-point count. |
| `n_full_kpoints` | No | Expected full-grid k-point count. |
| `rtol` | No | Relative tolerance for floating-point tokens; defaults to `1e-8`. |
| `atol` | No | Absolute tolerance for floating-point tokens; defaults to `1e-8`. |

Missing source or reference directories, and empty reference directories,
raise `FileNotFoundError`.
Malformed manifests also cause errors. Skipping declared cases is not supported;
an unfinished case makes the test run fail until its data are supplied.

The runner discovers expected files directly from `ref/<case-name>`. Generated
and reference files must have exactly the same relative filenames. A missing
output, an output without a reference, or a content mismatch fails the case.
`cases.json` does not list output filenames.

## Comparisons

The current comparator handles whitespace-separated text files. It ignores
whitespace and checks the token count. Text tokens and pairs of integer tokens
must match exactly. Other numeric tokens must be finite and satisfy
`abs(actual - reference) <= max(atol, rtol * max(abs(actual), abs(reference)))`.
Fortran `D` exponents are accepted.

For `stru_out.txt`, lattice vectors and atoms are compared in their original order.
Symmetry operations are normalized to the `row` convention: a `col` rotation
matrix is transposed, while its translation components are unchanged. Operations
are then matched one-to-one without regard to order. Rotation matrices must match
exactly; translations use the numerical tolerances above. Counts must agree, and
duplicate operations cannot replace missing operations. Unknown conventions and
incomplete blocks cause errors. Other output files retain token-order comparison.

The default `rtol=1e-8` and `atol=1e-8` allow for differences in the precision
of calculation inputs and reference data. The number of digits printed by the
converter does not determine the precision of its inputs. For case `001`, lattice
components differ by up to `4.62e-9` Bohr (about `1.35e-9` relative); `1e-12`
rejects this difference. These tolerances do not allow different atomic positions
or symmetry operations. Any case-specific tolerance must have a reason based
on source precision; do not loosen it to hide a conversion error.
The failure message gives the file, token index, actual value, and reference
value.

## Run

Install the test dependencies and run from the repository root:

```bash
python -m pip install -e ".[test]"
python -m pytest tests/test_cases.py -v
```

Select one case by name:

```bash
python -m pytest tests/test_cases.py -k 001_abacus_LiH_sym -v
```

Check the framework's comparisons, manifest validation, and missing-data errors:

```bash
python -m pytest tests/test_case_runner.py -v
```

Run all tests with `python -m pytest -v`. The existing
`python -m unittest discover -s tests -v` entry point remains available. All conversion tests use the real-data cases; framework unit tests
check comparison and error handling separately.

To add a case, prepare `source/<case-name>` and independently verified output
files in `ref/<case-name>`, then add its metadata to `cases.json`. Register new
source classes in `case_runner.SOURCE_CLASSES` when their implementations become
available. No placeholder FHI-aims conversion is provided.

Case `001` compares both `stru_out.txt` and `bz_sample.txt`. The BZ reference
comes from the supplied ABACUS-native output.
