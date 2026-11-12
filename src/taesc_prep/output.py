"""Write TAESC input files from extracted calculation data."""

from collections import Counter
import math
from pathlib import Path

from .structure import multiply_row, reciprocal_vectors


def write_taesc_stru(structure, output):
    """Write a TAESC structure file from a dict and return its path.

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


def write_taesc_bz_sampling(sampling, output):
    """Write normalized BZ sampling and its optional explicit full-grid map.

    Required keys: mesh, direct, weights, full_direct, full_map, and lattice.
    Fractional k vectors are converted with the reciprocal lattice including
    2*pi. Mapping indices and the two representative columns are one-based.
    """
    reciprocal = reciprocal_vectors(sampling["lattice"])
    direct, weights = sampling["direct"], sampling["weights"]
    full, mapping = sampling["full_direct"], sampling["full_map"]
    mesh = sampling["mesh"]
    n_reduced = len(direct)
    if (len(mesh) != 3 or any(type(n) is not int or n <= 0 for n in mesh)
            or math.prod(mesh) != len(full) or len(mapping) != len(full)
            or not 0 < n_reduced <= len(full) or len(weights) != n_reduced):
        raise ValueError("Invalid BZ sampling dimensions")
    if any(type(i) is not int or not 1 <= i <= n_reduced for i in mapping):
        raise ValueError("Invalid BZ sampling representative index")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Invalid BZ sampling weight")
    total = sum(weights)
    if not math.isfinite(total):
        raise ValueError("Non-finite total BZ sampling weight")
    normalized = [w / total for w in weights]
    multiplicity = Counter(mapping)
    if any(abs(w - multiplicity[i] / len(full)) > 1e-6
           for i, w in enumerate(normalized, 1)):
        raise ValueError("BZ sampling weight disagrees with full-grid multiplicity")

    def cartesian(point):
        if len(point) != 3 or not all(math.isfinite(x) for x in point):
            raise ValueError("Invalid BZ sampling coordinate")
        return multiply_row(point, reciprocal)

    rows = [" ".join(map(str, mesh)), f"{n_reduced} {n_reduced}"]
    for index, (point, weight) in enumerate(zip(direct, normalized), 1):
        values = [weight, *point, *cartesian(point)]
        rows.append(f"{index} " + " ".join(f"{x:.16e}" for x in values) + f" {index} {index}")
    if len(full) > n_reduced:
        rows.append(f"full_kmap {len(full)}")
        for index, (point, representative) in enumerate(zip(full, mapping), 1):
            rows.append(f"{index} " + " ".join(f"{x:.16e}" for x in cartesian(point))
                        + f" {representative}")
    path = Path(output)
    path.write_text("\n".join(rows) + "\n")
    return path
