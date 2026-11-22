"""Spinor (SU(2)) representations of spatial symmetry operations.

ABACUS's nspin=4 (SOC) producer exports a ``spin_symmetry`` block next to the
spatial operations. Each row carries the antiunitary flag and the eight real
components of a 2x2 complex matrix, stored as
``(Re u00, Im u00, Re u01, Im u01, Re u10, Im u10, Re u11, Im u11)``.

The matrices follow ``U sigma_i U^dagger = sum_j W_ij sigma_j`` with
``W = R_proper^T``, which is the convention LibRPA's spinor reader expects.
This module mirrors ABACUS's ``ModuleSymmetry::SpinRotation::so3_to_su2``.
"""

import math

__all__ = ["rotation_to_su2", "proper_part"]

EPS = 1.0e-10


def proper_part(rotation):
    """Return the proper (det=+1) part of an orthogonal Cartesian rotation."""
    det = (rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
           - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
           + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]))
    if det < 0.0:
        return [[-x for x in row] for row in rotation]
    return [list(row) for row in rotation]


def _su2_from_axis_angle(nx, ny, nz, theta):
    c = math.cos(0.5 * theta)
    s = math.sin(0.5 * theta)
    return (complex(c, -s * nz), complex(-s * ny, -s * nx),
            complex(s * ny, -s * nx), complex(c, s * nz))


def rotation_to_su2(cartesian, eps=EPS):
    """Convert a Cartesian rotation to its SU(2) representative.

    ``cartesian`` is a 3x3 row-major matrix in Cartesian coordinates.
    Return ``(u00, u01, u10, u11)`` as complex numbers.
    """
    r = proper_part(cartesian)
    trace = r[0][0] + r[1][1] + r[2][2]
    cos_theta = max(-1.0, min(1.0, 0.5 * (trace - 1.0)))
    theta = math.acos(cos_theta)
    sin_theta = math.sin(theta)

    if theta < eps:
        return (complex(1.0, 0.0), complex(0.0, 0.0),
                complex(0.0, 0.0), complex(1.0, 0.0))

    if abs(theta - math.pi) < eps or abs(sin_theta) < eps:
        nn = [0.5 * (r[i][i] + 1.0) for i in range(3)]
        nn = [x if x > 0.0 else 0.0 for x in nn]
        imax = nn.index(max(nn))
        n = [0.0, 0.0, 0.0]
        n[imax] = math.sqrt(nn[imax])
        for j in range(3):
            if j == imax:
                continue
            off = 0.25 * (r[imax][j] + r[j][imax])
            n[j] = off / n[imax] if n[imax] > 0.0 else 0.0
        norm = math.sqrt(sum(x * x for x in n))
        if norm > 0.0:
            n = [x / norm for x in n]
        return _su2_from_axis_angle(n[0], n[1], n[2], math.pi)

    inv = 1.0 / (2.0 * sin_theta)
    nx = (r[1][2] - r[2][1]) * inv
    ny = (r[2][0] - r[0][2]) * inv
    nz = (r[0][1] - r[1][0]) * inv
    return _su2_from_axis_angle(nx, ny, nz, theta)


def flatten_su2(u):
    """Flatten (u00, u01, u10, u11) to the eight real values ABACUS writes."""
    out = []
    for value in u:
        out.extend((value.real, value.imag))
    return out


def cartesian_rotation(lattice, fractional_rotation):
    """``latvec^-1 * R_frac * latvec`` for row-vector fractional operations."""
    inverse = _inverse3(lattice)
    tmp = _matmul(fractional_rotation, lattice)
    return _matmul(inverse, tmp)


def _matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def _inverse3(m):
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    if not math.isfinite(det) or abs(det) < 1.0e-30:
        raise ValueError("Lattice vectors are singular")
    cof = [[(m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3]
             - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]) / det
            for j in range(3)] for i in range(3)]
    return cof
