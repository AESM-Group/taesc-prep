"""Lattice and coordinate conversion helpers."""

import math


def multiply_row(vector, matrix):
    return [sum(vector[i] * matrix[i][j] for i in range(3)) for j in range(3)]



def fractional_coordinates(xyz, reciprocal):
    return [sum(x*y for x, y in zip(xyz, row)) / (2*math.pi)
            for row in reciprocal]



def reciprocal_vectors(lattice):
    """Return row vectors B satisfying A B^T = 2*pi I."""
    def cross(a, b):
        return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]]

    a, b, c = lattice
    cofactors = [cross(b, c), cross(c, a), cross(a, b)]
    volume = sum(x*y for x, y in zip(a, cofactors[0]))
    if not math.isfinite(volume) or volume == 0:
        raise ValueError("Lattice vectors must define a finite, nonzero volume")
    result = [[2 * math.pi * x / volume for x in row] for row in cofactors]
    if not all(math.isfinite(x) for row in result for x in row):
        raise ValueError("Non-finite reciprocal lattice")
    return result
