"""Common interface for electronic-structure calculation sources."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .output import write_taesc_bz_sampling, write_taesc_stru


@dataclass
class ConversionResult:
    """Destination directory, structure and k-point counts, and conversion notes."""

    output: Optional[Path] = None
    n_atoms: int = 0
    n_symops: int = 0
    notes: tuple[str, ...] = ()
    n_kpoints: int = 0
    n_full_kpoints: int = 0


class Source(ABC):
    """A calculation directory that supplies input data for TAESC."""

    process_stru = False
    stru_filename = None

    process_bz_sampling = False
    bz_sampling_filename = None

    def __init__(self, directory):
        self.directory = Path(directory)

    def convert(self, output=None, *, include_symmetry=True) -> ConversionResult:
        """Run enabled conversions and return their result.

        Sources enable conversion with ``process_stru`` and
        ``process_bz_sampling`` during initialization. ``output`` selects a destination directory
        and defaults to the input directory. If no conversion is enabled,
        return an empty result without reading or writing files.
        """
        result = ConversionResult()
        structure = None
        output_dir = Path(output) if output is not None else self.directory
        if self.process_stru:
            structure = self.convert_stru(output_dir, include_symmetry=include_symmetry)
            result.output = output_dir
            result.n_atoms = len(structure["types"])
            result.n_symops = len(structure.get("operations", ()))
            result.notes = tuple(structure.get("notes", ()))
        if self.process_bz_sampling:
            sampling = self.convert_bz_sampling(output_dir, structure=structure)
            result.output = output_dir
            result.n_kpoints = len(sampling["direct"])
            result.n_full_kpoints = len(sampling["full_direct"])
            result.notes += tuple(sampling.get("notes", ()))
        return result

    def convert_stru(self, output=None, *, include_symmetry=True) -> dict:
        """Extract and write the structure file, then return the extracted dictionary.

        Write ``stru_filename`` inside ``output``, which defaults to the input
        directory. Create the destination directory if necessary. Relative
        directories are resolved from the working directory. Set ``include_symmetry=False``
        to omit spatial symmetry operations. ``convert`` builds the result.
        """
        structure = self._extract_stru(include_symmetry=include_symmetry)
        output_path = self._output_path(output, self.stru_filename, structure)
        write_taesc_stru(structure, output_path)
        return structure

    def convert_bz_sampling(self, output=None, *, structure=None) -> dict:
        """Extract BZ sampling, add the structure's lattice, and write its file."""
        if structure is None:
            structure = self._extract_stru(include_symmetry=False)
        sampling = self._extract_bz_sampling()
        sampling["lattice"] = structure["lattice"]
        sampling["source_paths"] = (*sampling.get("source_paths", ()),
                                    *structure.get("source_paths", ()))
        sampling["notes"] = tuple(sampling.get("notes", ()))
        if not self.process_stru:
            sampling["notes"] += tuple(structure.get("notes", ()))
        output_path = self._output_path(output, self.bz_sampling_filename, sampling)
        write_taesc_bz_sampling(sampling, output_path)
        return sampling

    def _output_path(self, output, filename, data):
        output_dir = Path(output) if output is not None else self.directory
        if output_dir.exists() and not output_dir.is_dir():
            raise ValueError(f"Output must be a directory: {output_dir}")
        output_path = output_dir / filename
        if any(Path(path).resolve() == output_path.resolve()
               for path in data.get("source_paths", ())):
            raise ValueError("Input and output paths must differ")
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_path

    @abstractmethod
    def _extract_bz_sampling(self) -> dict:
        """Return mesh, direct, weights, full_direct, and one-based full_map.

        Optional notes and source_paths have the same meaning as in _extract_stru.
        Only sources that set process_bz_sampling need to implement this method.
        """
        raise NotImplementedError

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
