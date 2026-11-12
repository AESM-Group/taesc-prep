# Implementation details

## ABACUS

The program writes `/path/to/calculation/stru_out`. It selects geometry in this
order:

1. `OUT.<suffix>/STRU_FINAL`.
2. The input structure specified by `stru_file`, or `STRU` by default.
3. A single-geometry running log (such as `running_scf.log`) in `OUT.<suffix>`,
   then in the supplied directory.

The `suffix`, `calculation`, and `stru_file` settings come from `INPUT` when
available. The suffix defaults to `ABACUS`. Structure files take priority over
logs to retain their numerical precision. Logs are used only when neither
structure file exists; ambiguous logs cause an error.

The structure reader supports explicit `LATTICE_VECTORS` and the `Direct`,
`Cartesian`, `Cartesian_au`, `Cartesian_angstrom`, and
`Cartesian_angstrom_center_*` coordinate modes. It preserves species and atom
order. For both `STRU_FINAL` and the input structure, lattice vectors in Bohr
are `LATTICE_VECTORS` multiplied by `LATTICE_CONSTANT`. Before converting
fractional positions to Cartesian coordinates, components outside `[0, 1)`
are translated by integer lattice periods into that interval (for example,
`-0.25` becomes `0.75`, and `1.0` becomes `0.0`). Components already within
the interval are unchanged.

Only Angstrom coordinate modes use `BOHR_TO_A = 0.5291770` to convert atomic
coordinates to Bohr. For centered coordinates, the reader selects the requested
Cartesian components of the cell center, as ABACUS does.

For relaxation or dynamics, the selected structure must correspond to the
exported electronic data. Structures specified through `latname` need an
explicit-lattice snapshot or the log route.

The program includes spatial symmetry operations when available in the
ABACUS output. These operations use the
fractional `row` convention. Rotations and translations retain the precision
of the ABACUS export. A symmetry count in `running_*.log` is insufficient
to reconstruct the operations. The program does not generate operations
or spin-symmetry metadata. All input files must belong to the same geometry.
No k-point block is written.

An existing `stru_out` is replaced after parsing succeeds. The program uses
only the Python standard library. The current reader supports ABACUS;
FHI-aims directory support is reserved for a later extension.

The ABACUS conventions follow its
[coordinate conversion source](https://github.com/deepmodeling/abacus-develop/blob/develop/source/source_cell/read_atoms_helper.cpp)
and [symmetry export source](https://github.com/deepmodeling/abacus-develop/blob/develop/source/source_lcao/module_ri/module_exx_symmetry/symm_rot_out.cpp).



## Output conventions

- Lattice vectors and Cartesian positions are in Bohr.
- Reciprocal lattice vectors are in inverse Bohr and include `2*pi`:
  `A @ B.T = 2*pi*I`, with lattice vectors stored as rows.
- Atom types start at 1 and follow the ABACUS species order.
- The optional spatial symmetry block uses fractional row-vector operations:
  `x_new = x @ R + t`.

See the [LibRPA dataset format](https://minyez.github.io/LibRPA/develop/dataset_format.html#stru-out)
for the file specification.

## Python interface

`detect_source(directory)` identifies the calculation and returns a `Source`
instance. The CLI uses this function before conversion. ABACUS detection
reuses `discover_abacus`, including the `stru_file`, `suffix`, and `calculation`
settings from `INPUT`. A directory without discoverable ABACUS input produces
an error.

`Source` implements the common conversion workflow. `AbacusSource` supplies
ABACUS data through `_extract_stru`:

```python
from taesc_prep import Source, detect_source

source: Source = detect_source("/path/to/calculation")
result = source.convert()
print(result.output)
```

`convert(output=None, *, include_symmetry=True)` runs the enabled conversions
and returns a `ConversionResult` with `output`, `n_atoms`, `n_symops`, and `notes`.
`output` is a destination directory and defaults to the input directory.
The returned `result.output` is that directory; the structure file is
`result.output / "stru_out"`. Set `include_symmetry=False` to omit spatial
symmetry operations. If no conversion is enabled, the result has `output=None`,
zero counts, and no notes; no files are read or written.

`Source` initializes `process_stru=False`. A concrete source enables structure
processing in its initializer and implements `_extract_stru`:

```python
class MySource(Source):
    def __init__(self, directory):
        super().__init__(directory)
        self.process_stru = True

    def _extract_stru(self, *, include_symmetry=True) -> dict:
        # Read the calculation files and return the dictionary described below.
        ...
```

Concrete sources inherit `convert` and `convert_stru` without overriding them.
When `process_stru` is True, `convert` calls the shared `convert_stru` method.
That method calls `_extract_stru`, passes its dictionary to `write_taesc_stru`
in `output.py`, and returns the extracted dictionary. Only `convert` constructs
and returns `ConversionResult`.

`_extract_stru` returns a dictionary with the following keys:

| Key | Required | Value |
| --- | --- | --- |
| `lattice` | Yes | Three row vectors, each with three components, in Bohr. |
| `positions` | Yes | One Cartesian coordinate row per atom, in Bohr. |
| `types` | Yes | One integer type per atom, numbered from 1. |
| `operations` | No | `(rotation, translation)` pairs in the fractional row convention. Each rotation is a flat list of nine integers in row-major order; each translation has three components. |
| `notes` | No | Conversion notes as a sequence of strings. |
| `source_paths` | No | Input paths as strings or `Path` objects, protected from overwriting. |

Missing optional keys default to empty sequences. The writer computes reciprocal
vectors from `lattice`. It requires the three geometry keys and uses `operations`
when present. The shared conversion method handles `notes` and `source_paths`.
