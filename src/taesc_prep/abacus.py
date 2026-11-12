"""Read ABACUS structures, running logs, and spatial symmetry operations."""

import math
import re

from .source import Source
from .structure import fractional_coordinates, multiply_row, reciprocal_vectors


# ABACUS source_base/constants.h: BOHR_TO_A.
BOHR_TO_ANGSTROM = 0.5291770
SECTIONS = {
    "ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "NUMERICAL_DESCRIPTOR",
    "PAW_FILES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "LATTICE_PARAMETERS",
    "ATOMIC_POSITIONS",
}


class AbacusSource(Source):
    """Convert the automatically discovered files in an ABACUS directory."""

    def __init__(self, directory):
        super().__init__(directory)
        self.process_stru = True

    def _extract_stru(self, *, include_symmetry=True) -> dict:
        """Read ABACUS geometry and optional spatial symmetry operations."""
        discovered = discover_abacus(self.directory)
        if discovered is None:
            raise ValueError("No ABACUS STRU or running log found; check the directory and stru_file in INPUT")
        source, fmt, outdir = discovered
        text = source.read_text()
        structure = read_stru(text) if fmt == "stru" else read_log(text)
        symfile = None
        if include_symmetry:
            candidates = list(dict.fromkeys([self.directory / "symrot_R.txt",
                                            outdir / "symrot_R.txt"]))
            found = [p for p in candidates if p.is_file()]
            if len(found) > 1:
                raise ValueError("Multiple spatial symmetry exports found; use a directory containing one calculation")
            symfile = found[0] if found else None
        if symfile:
            structure["operations"] = read_symrot(symfile.read_text())
        structure["source_paths"] = (source, symfile) if symfile else (source,)
        if fmt == "log":
            structure["notes"] = ("geometry precision is limited by the running log; prefer STRU.",)
        return structure


def numbers(line, count=3):
    tokens = line.split()
    if len(tokens) < count:
        raise ValueError(f"Expected {count} numbers: {line!r}")
    values = [float(x.replace("D", "E").replace("d", "e"))
              for x in tokens[:count]]
    if not all(math.isfinite(x) for x in values):
        raise ValueError(f"Non-finite number: {line!r}")
    return values



def read_stru(text):
    """Read explicit lattice vectors and preserve ABACUS species/atom order."""
    sections = {}
    current = None
    for raw in text.splitlines():
        line = re.split(r"#|//", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        if line in SECTIONS:
            if line in sections:
                raise ValueError(f"Duplicate STRU section: {line}")
            current = sections[line] = []
        elif current is not None:
            current.append(line)
        else:
            raise ValueError(f"Expected a STRU section: {line!r}")

    for key in ("ATOMIC_SPECIES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS"):
        if not sections.get(key):
            raise ValueError(f"Missing {key}; use an explicit-lattice STRU or a running log")
    if len(sections["LATTICE_CONSTANT"]) != 1 or len(sections["LATTICE_VECTORS"]) != 3:
        raise ValueError("Expected one lattice constant and three lattice vectors")
    scale = numbers(sections["LATTICE_CONSTANT"][0], 1)[0]
    if scale <= 0:
        raise ValueError("LATTICE_CONSTANT must be positive")
    lattice = [[scale*x for x in numbers(row)] for row in sections["LATTICE_VECTORS"]]
    reciprocal = reciprocal_vectors(lattice)
    species = [row.split()[0] for row in sections["ATOMIC_SPECIES"]]
    if len(set(species)) != len(species):
        raise ValueError("Duplicate ATOMIC_SPECIES label")
    rows = iter(sections["ATOMIC_POSITIONS"])
    mode = next(rows).lower()
    factors = {"cartesian": scale, "cartesian_au": 1.0,
               "cartesian_angstrom": 1 / BOHR_TO_ANGSTROM}
    shift = [0.0]*3
    if mode.startswith("cartesian_angstrom_center_"):
        center = mode[len("cartesian_angstrom_center_"):]
        if center not in ("xy", "xz", "yz", "xyz"):
            raise ValueError(f"Unsupported coordinate mode: {mode}")
        # ABACUS selects Cartesian components of the full cell center.
        shift = [sum(row[i] for row in lattice)/2 if axis in center else 0.0
                 for i, axis in enumerate("xyz")]
        mode = "cartesian_angstrom"
    if mode != "direct" and mode not in factors:
        raise ValueError(f"Unsupported coordinate mode: {mode}")
    positions, types = [], []
    try:
        for it, label in enumerate(species, 1):
            if next(rows) != label:
                raise ValueError(f"Expected ATOMIC_POSITIONS species {label}")
            numbers(next(rows), 1)  # Species magnetization.
            count = int(next(rows))
            if count <= 0:
                raise ValueError(f"Atom count for {label} must be positive")
            for _ in range(count):
                xyz = numbers(next(rows))
                if mode == "direct":
                    frac = xyz
                else:
                    xyz = [x*factors[mode] + s for x, s in zip(xyz, shift)]
                    frac = fractional_coordinates(xyz, reciprocal)
                # ABACUS check_dtau moves atoms into the reference unit cell.
                xyz = multiply_row([x % 1.0 for x in frac], lattice)
                positions.append(xyz)
                types.append(it)
    except StopIteration as exc:
        raise ValueError("Incomplete ATOMIC_POSITIONS section") from exc
    if next(rows, None) is not None:
        raise ValueError("Unexpected data after ATOMIC_POSITIONS")
    return {"lattice": lattice, "positions": positions, "types": types}



def read_log(text):
    """Read the setup geometry from a log with one atomic coordinate table.

    ABACUS calls its lattice scaling factor a_0 in the lattice-table header.
    Thus these vectors, and Cartesian atom coordinates, need that factor.
    """
    lines = text.splitlines()

    def unique_match(pattern):
        matches = list(re.finditer(pattern, text, re.MULTILINE | re.IGNORECASE))
        if len(matches) != 1:
            raise ValueError(f"Expected one log entry matching {pattern!r}; found {len(matches)}")
        return matches[0]

    # A relaxation log can contain later positions without another setup table.
    ionic_steps = [int(x) for x in re.findall(r"\bION\s*=\s*(\d+)", text)]
    if any(x > 1 for x in ionic_steps) or re.search(r"\b(?:STEP OF RELAXATION|MD STEP)\b", text):
        raise ValueError("Use a snapshot STRU for relaxation or molecular-dynamics output")
    scale = numbers(unique_match(r"^\s*lattice constant \(Bohr\)\s*=\s*(\S+)").group(1), 1)[0]
    if scale <= 0:
        raise ValueError("Lattice constant must be positive")
    nat = int(unique_match(r"^\s*TOTAL ATOM NUMBER\s*=\s*(\d+)").group(1))
    blocks = re.split(r"(?im)^\s*READING ATOM TYPE\s+(\d+)\s*$", text)
    types, labels = [], []
    for index in range(1, len(blocks), 2):
        it = int(blocks[index])
        body = blocks[index + 1]
        label = re.search(r"(?i)atom label\s*=\s*(\S+)", body)
        count = re.search(r"(?i)number of atoms? for this type\s*=\s*(\d+)", body)
        if not label or not count or it != len(labels) + 1:
            raise ValueError("Incomplete or out-of-order atom-type information in log")
        labels.append(label.group(1))
        types.extend([it] * int(count.group(1)))
    if nat <= 0 or len(types) != nat:
        raise ValueError("Atom-type counts disagree with TOTAL ATOM NUMBER")
    lattice_headers = [i for i, line in enumerate(lines)
                       if line.strip() == "Lattice vectors: (Cartesian coordinate: in unit of a_0)"]
    coord_headers = [i for i, line in enumerate(lines)
                     if re.match(r"^\s*(DIRECT|CARTESIAN) COORDINATES(?:\s|$)", line)]
    if len(lattice_headers) != 1 or len(coord_headers) != 1:
        raise ValueError("Expected one lattice table and one atomic coordinate table; use a snapshot STRU")
    ilat = lattice_headers[0]
    if len(lines[ilat+1:ilat+4]) != 3:
        raise ValueError("Incomplete lattice table")
    lattice = [[scale*x for x in numbers(row)] for row in lines[ilat+1:ilat+4]]
    ic = coord_headers[0]
    direct = lines[ic].strip().startswith("DIRECT")
    rows = lines[ic+2:ic+2+nat]  # Skip column headings.
    if len(rows) != nat:
        raise ValueError("Incomplete atomic coordinate table")
    positions = []
    for row, it in zip(rows, types):
        fields = row.split()
        label = labels[it - 1]
        if not fields or (fields[0] != label and
                          not re.fullmatch(r"tau[cd]_" + re.escape(label) + r"\d+", fields[0])):
            raise ValueError(f"Expected coordinates for {label}: {row!r}")
        xyz = numbers(" ".join(fields[1:]))
        positions.append(multiply_row(xyz, lattice) if direct else [scale*x for x in xyz])
    return {"lattice": lattice, "positions": positions, "types": types}



def read_symrot(text):
    """Extract fractional row-convention operations from the ABACUS spatial symmetry export."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    lmax = []
    for label in ("AOs", "ABFs"):
        matches = re.findall(r"Lmax of " + label + r":\s*(\d+)", text)
        if len(matches) != 1:
            raise ValueError(f"Expected one Lmax of {label} in the spatial symmetry export")
        lmax.append(int(matches[0]))
    start = next((i for i, line in enumerate(lines) if line == "0"), None)
    if start is None:
        raise ValueError("No symmetry operations in the spatial symmetry export")
    # Each angular momentum contributes 2*l+1 orbital-matrix rows.
    orbital_rows = (max(lmax) + 1)**2
    operations = []
    pos = start
    while pos < len(lines):
        if lines[pos] != str(len(operations)) or pos + 5 + orbital_rows > len(lines):
            raise ValueError("Incomplete or out-of-order symmetry operation")
        rotation = [x for row in lines[pos+1:pos+4] for x in numbers(row)]
        if any(abs(x-round(x)) > 1e-8 for x in rotation):
            raise ValueError("Noninteger fractional symmetry rotation")
        rotation = [int(round(x)) for x in rotation]
        translation = lines[pos+4]
        if not (translation.startswith("(") and translation.endswith(")")):
            raise ValueError("Expected a parenthesized symmetry translation")
        translation = numbers(translation[1:-1])
        operations.append((rotation, translation))
        pos += 5 + orbital_rows
    identity = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    if not any(r == identity and all(abs(x-round(x)) < 1e-8 for x in t)
               for r, t in operations):
        raise ValueError("Symmetry operations do not contain the identity")
    return operations



def discover_abacus(directory):
    """Return (source path, format, output directory), or None if absent.

    Raise ValueError for an invalid directory or ambiguous running logs.
    """
    if not directory.is_dir():
        raise ValueError(f"Expected a calculation directory: {directory}")
    settings = {}
    input_file = directory / "INPUT"
    if input_file.is_file():
        for raw in input_file.read_text().splitlines():
            fields = re.split(r"#|//", raw, maxsplit=1)[0].split()
            if len(fields) >= 2:
                settings[fields[0].lower()] = fields[1]
    stru = directory / settings.get("stru_file", "STRU")
    outdir = directory / ("OUT." + settings.get("suffix", "ABACUS"))
    if stru.is_file():
        return stru, "stru", outdir
    calculation = settings.get("calculation")
    logs = ([outdir / f"running_{calculation}.log"] if calculation
            else sorted(outdir.glob("running_*.log")))
    logs = [path for path in logs if path.is_file()]
    if not logs:
        logs = sorted(directory.glob("running_*.log"))
    if len(logs) == 1:
        return logs[0], "log", logs[0].parent
    if len(logs) > 1:
        raise ValueError("Multiple running logs found; specify calculation in INPUT or use a directory containing one log")
    return None

