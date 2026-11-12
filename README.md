# taesc-prep

`taesc-prep` prepares input datasets for TAESC and LibRPA from first-principles
calculations.

## Installation

Python 3.9 or later is required. There are no runtime dependencies.

Create a virtual environment inside the repository and install an editable copy:

```bash
cd ~/codes/taesc-prep
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The environment is stored in `~/codes/taesc-prep/.venv`. The editable
installation picks up source changes automatically.

After activation, run:

```bash
taesc-prep /path/to/calculation
```

Without activating the environment, use the full executable path:

```bash
~/codes/taesc-prep/.venv/bin/taesc-prep /path/to/calculation
```

For a non-editable installation, use `python -m pip install .` instead of
`python -m pip install -e .`.

## Usage

Run the Python preprocessor with the calculation directory as its only
positional argument:

```bash
taesc-prep /path/to/calculation
```

The program detects the source from the directory contents and prepares the
input data.

Use `--no-symmetry` to omit spatial symmetry operations:

```bash
taesc-prep /path/to/calculation --no-symmetry
```

The same command is available as `python3 -m taesc_prep`.

## Supported software and inputs

`X` marks the input file that is handled by `taesc-prep`.
Other essential dataset inputs are generated either directly by the software,
or by another package, e.g. pyatb.

|             | ABACUS | FHI-aims |
| ----------- | ------ | -------- |
| Structure   | X      |          |
| BZ sampling | X      |          |

## Development

Install the test dependencies and run the tests:

```bash
python -m pip install -e ".[test]"
python -m pytest -v
```

See [Real-calculation tests](tests/README.md) for the case manifest and reference-data workflow.

See [Implementation details](docs/implement.md) for source parsing, conversion
interfaces, and output conventions.
