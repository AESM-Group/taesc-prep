"""Common interface for electronic-structure calculation sources."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .output import write_taesc_stru


@dataclass
class ConversionResult:
    """Destination directory, atom/symmetry counts, and conversion notes."""

    output: Optional[Path] = None
    n_atoms: int = 0
    n_symops: int = 0
    notes: tuple[str, ...] = ()


class Source(ABC):
    """A calculation directory that supplies input data for TAESC."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self.process_stru = False

    def convert(self, output=None, *, include_symmetry=True) -> ConversionResult:
        """Run enabled conversions and return their result.

        Sources enable structure conversion by setting ``process_stru`` to
        True during initialization. ``output`` selects a destination directory
        and defaults to the input directory. If no conversion is enabled,
        return an empty result without reading or writing files.
        """
        result = ConversionResult()
        if self.process_stru:
            output_dir = Path(output) if output is not None else self.directory
            structure = self.convert_stru(output_dir, include_symmetry=include_symmetry)
            result.output = output_dir
            result.n_atoms = len(structure["types"])
            result.n_symops = len(structure.get("operations", ()))
            result.notes = tuple(structure.get("notes", ()))
        return result

    def convert_stru(self, output=None, *, include_symmetry=True) -> dict:
        """Extract and write stru_out, then return the extracted dictionary.

        Write ``stru_out`` inside ``output``, which defaults to the input
        directory. Create the destination directory if necessary. Relative
        directories are resolved from the working directory. Set ``include_symmetry=False``
        to omit spatial symmetry operations. ``convert`` builds the result.
        """
        structure = self._extract_stru(include_symmetry=include_symmetry)
        output_dir = Path(output) if output is not None else self.directory
        if output_dir.exists() and not output_dir.is_dir():
            raise ValueError(f"Output must be a directory: {output_dir}")
        output_path = output_dir / "stru_out"
        if any(Path(path).resolve() == output_path.resolve()
               for path in structure.get("source_paths", ())):
            raise ValueError("Input and output paths must differ")
        output_dir.mkdir(parents=True, exist_ok=True)
        write_taesc_stru(structure, output_path)
        return structure

    @abstractmethod
    def _extract_stru(self, *, include_symmetry=True) -> dict:
        """Return a dict with required geometry and optional metadata.

        Required keys: ``lattice`` (3x3 row vectors in Bohr), ``positions``
        (Nx3 Cartesian coordinates in Bohr), and ``types`` (N atom types,
        numbered from 1).

        Optional keys default to empty sequences: ``operations`` contains
        (rotation, translation) pairs in the fractional row convention, with
        each rotation flattened in row-major order; ``notes`` contains strings;
        ``source_paths`` contains input paths protected from overwriting.

        Concrete sources read their calculation files here without writing
        output. Raise OSError or ValueError when input cannot be read or converted.
        """
        raise NotImplementedError
