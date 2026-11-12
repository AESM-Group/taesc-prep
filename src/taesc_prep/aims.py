"""Read FHI-aims calculation directory."""

from .source import Source

__all__ = ["AimsSource", "discover_aims"]


class AimsSource(Source):

    def __init__(self, directory):
        super().__init__(directory)


def discover_aims(directory) -> bool:
    """Identify a FHI-aims directory from file markers without selecting inputs."""
    return False
