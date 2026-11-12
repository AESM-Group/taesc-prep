"""Detect the calculation source from files in a directory."""

from pathlib import Path

from .source import Source
from .abacus import AbacusSource, discover_abacus


def detect_source(directory) -> Source:
    """Return a source instance without converting the calculation data.

    Recognize ABACUS input directories and standalone running-log directories.
    Reject unknown or mixed sources. FHI-aims markers are recognized, but its
    source implementation is not available yet.
    """
    directory = Path(directory)
    abacus = discover_abacus(directory) is not None
    if abacus:
        return AbacusSource(directory)
    raise ValueError(f"Cannot detect a supported calculation source in {directory}")
