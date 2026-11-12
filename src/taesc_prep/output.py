"""Write TAESC input files from extracted calculation data."""

import math
from pathlib import Path

from .structure import reciprocal_vectors


def write_taesc_stru(structure, output):
    """Write a TAESC stru_out file from a dict and return its path.

    ``output`` is the file path selected by the conversion method.
    Required keys are ``lattice``, ``positions``, and ``types``. The optional
    ``operations`` key supplies spatial symmetry. Other metadata is ignored.

    Format all data before replacing an existing file. Lattice vectors and
    positions are in Bohr; reciprocal vectors include 2*pi. Spatial symmetry
    operations use the fractional row convention.
    """
    missing = {"lattice", "positions", "types"} - structure.keys()
    if missing:
        raise ValueError("Missing required structure keys: " + ", ".join(sorted(missing)))
    lattice = structure["lattice"]
    positions = structure["positions"]
    types = structure["types"]
    operations = structure.get("operations", ())
    reciprocal = reciprocal_vectors(lattice)
    rows = [" ".join(f"{x:.16e}" for x in row)
            for row in lattice + reciprocal]
    rows.append(str(len(types)))
    for xyz, it in zip(positions, types):
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError("Non-finite Cartesian position")
        rows.append(" ".join(f"{x:.16e}" for x in xyz) + f" {it}")
    if operations:
        rows.append(f"{len(operations)} row")
        for rotation, translation in operations:
            rows.append(" ".join(str(x) for x in rotation) + " "
                        + " ".join(f"{x:.16e}" for x in translation))
    path = Path(output)
    path.write_text("\n".join(rows) + "\n")
    return path

