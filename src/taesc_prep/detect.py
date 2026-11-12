"""Detect the calculation source from files in a directory."""

from pathlib import Path

from .source import Source
from .abacus import AbacusSource, discover_abacus
from .aims import AimsSource, discover_aims


def detect_source(directory) -> Source:
    """Identify the calculation source without choosing conversion inputs."""
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"Expected a calculation directory: {directory}")
    if discover_abacus(directory):
        return AbacusSource(directory)
    if discover_aims(directory):
        return AimsSource(directory)
    raise ValueError(
        f"Cannot detect a supported calculation source in {directory}")
