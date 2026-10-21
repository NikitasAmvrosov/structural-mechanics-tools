"""
Kern area, principal axes and neutral axis of a polygonal cross-section.

Coordinate system (as in German TM courses):
    y to the right, z downwards, x along the beam axis.
    Iy  =  integral(z^2 dA)
    Iz  =  integral(y^2 dA)
    Iyz = -integral(y*z dA)        (minus convention)

Bending normal stress for centroidal principal axes (Iyz = 0):
    sigma(y, z) = N/A + My/Iy * z - Mz/Iz * y
The neutral axis is the line where sigma = 0.

Usage:
    python kern_area.py                  # runs the example below and shows the plot
    python kern_area.py --save out.png   # saves the plot instead

Author: Nikitas Amvrosov
"""

import argparse
from math import atan2, cos, degrees, sin, sqrt

import matplotlib.pyplot as plt

EPS = 1e-9


# ------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------

def close_polygon(points):
    pts = list(points)
    if pts and pts[0] != pts[-1]:
        pts.append(pts[0])
    return pts


def unique_points(points, tol=1e-7):
    out = []
    for p in points:
        if all(abs(p[0] - q[0]) > tol or abs(p[1] - q[1]) > tol for q in out):
            out.append(p)
    return out


def sort_by_angle(points):
    if len(points) <= 1:
        return list(points)
    cy = sum(p[0] for p in points) / len(points)
    cz = sum(p[1] for p in points) / len(points)
    return sorted(points, key=lambda p: atan2(p[1] - cz, p[0] - cy))


def section_properties(points):
    """Area, centroid and centroidal second moments of a simple polygon.

    Uses the shoelace formula, so the polygon may be convex or not,
    and its vertices may be given clockwise or counter-clockwise.
    Returns (A, yS, zS, Iy, Iz, Iyz) with Iyz in the minus convention.
    """
    if len(points) < 3:
        raise ValueError("A polygon needs at least 3 points.")

    q = close_polygon(points)
    a2 = cy = cz = iy0 = iz0 = iyz0 = 0.0
    for (y1, z1), (y2, z2) in zip(q[:-1], q[1:]):
        c = y1 * z2 - y2 * z1
        a2 += c
        cy += (y1 + y2) * c
        cz += (z1 + z2) * c
        iy0 += (z1 * z1 + z1 * z2 + z2 * z2) * c
        iz0 += (y1 * y1 + y1 * y2 + y2 * y2) * c
        iyz0 += (2 * y1 * z1 + y1 * z2 + y2 * z1 + 2 * y2 * z2) * c

    if a2 < 0:  # clockwise input -> flip all signs
        a2, cy, cz, iy0, iz0, iyz0 = -a2, -cy, -cz, -iy0, -iz0, -iyz0

    area = 0.5 * a2
    ys = cy / (6.0 * area)
    zs = cz / (6.0 * area)

    # second moments about the origin, then shifted to the centroid (Steiner)
    iy = iy0 / 12.0 - area * zs * zs
    iz = iz0 / 12.0 - area * ys * ys
    int_yz = iyz0 / 24.0 - area * ys * zs
    iyz = -int_yz  # minus convention

    return area, ys, zs, iy, iz, iyz


def principal_axes(iy, iz, iyz):
    """Principal moments I1 >= I2 and angle phi (rad) of the I1 axis from y."""
    root = sqrt((iy - iz) ** 2 + 4.0 * iyz * iyz)
    i1 = 0.5 * (iy + iz) + 0.5 * root
    i2 = 0.5 * (iy + iz) - 0.5 * root
    if root < EPS * max(abs(iy), abs(iz), 1.0):
        phi = 0.0  # Iy = Iz and Iyz = 0: every axis is a principal axis
    else:
        phi = 0.5 * atan2(2.0 * iyz, iy - iz)
    return i1, i2, phi


def convex_hull(points):
    """Andrew's monotone chain. Returns hull vertices counter-clockwise."""
    pts = sorted(unique_points(points))
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


# ------------------------------------------------------------
# Kern area
# ------------------------------------------------------------

def kern_area(points_centroidal, area, iy, iz):
    """Kern (core) area of a section, in centroidal principal axes.

    Every edge of the convex hull is treated as a possible neutral axis.
    If it cuts the y-axis at y0 and the z-axis at z0, the matching load
    point on the kern boundary is
        yK = -iz^2 / y0,   zK = -iy^2 / z0,   with iy^2 = Iy/A, iz^2 = Iz/A.
    Only valid if y and z are principal axes (Iyz = 0).
    """
    hull = convex_hull(points_centroidal)
    if len(hull) < 3:
        return []

    iy2, iz2 = iy / area, iz / area
    kern = []
    h = close_polygon(hull)
    for (y1, z1), (y2, z2) in zip(h[:-1], h[1:]):
        # line a*y + b*z + c = 0 through both hull points
        a, b, c = z1 - z2, y2 - y1, y1 * z2 - y2 * z1
        if abs(c) < EPS:  # edge passes through the centroid -> no kern point
            continue
        yk = 0.0 if abs(a) < EPS else iz2 / c * a   # = -iz2 / y0 with y0 = -c/a
        zk = 0.0 if abs(b) < EPS else iy2 / c * b   # = -iy2 / z0 with z0 = -c/b
        kern.append((yk, zk))

    return sort_by_angle(unique_points(kern))


# ------------------------------------------------------------
# Neutral axis
# ------------------------------------------------------------

def neutral_axis(area, iy, iz, n, my, mz):
    """Coefficients (a, b, c) of the neutral axis a*y + b*z + c = 0.

    From sigma = N/A + My/Iy*z - Mz/Iz*y = 0. Only valid for Iyz = 0.
    """
    a = -mz / iz if abs(iz) > EPS else 0.0
    b = my / iy if abs(iy) > EPS else 0.0
    return a, b, n / area


def stress(y, z, area, iy, iz, n, my, mz):
    return n / area + my / iy * z - mz / iz * y


# ------------------------------------------------------------
# Main analysis + plot
# ------------------------------------------------------------

def analyse(points, n=0.0, my=0.0, mz=0.0, section=None, plot=True, save=None):
    """Run the full analysis.

    points  : section outline (y, z), any orientation, any origin
    n, my, mz: normal force and bending moments (consistent units, e.g. N, Nmm)
    section : optional (A, Iy, Iz, Iyz) to override the polygon values,
              e.g. for thin-walled profiles where the outline is only used
              for the kern area hull
    """
    a_poly, ys, zs, iy_poly, iz_poly, iyz_poly = section_properties(points)
    pc = [(y - ys, z - zs) for y, z in points]  # centroidal coordinates

    if section is None:
        area, iy, iz, iyz = a_poly, iy_poly, iz_poly, iyz_poly
    else:
        area, iy, iz, iyz = section

    i1, i2, phi = principal_axes(iy, iz, iyz)

    print("=== Section properties (centroidal axes) ===")
    print(f"A   = {area:.6g}")
    print(f"yS  = {ys:.6g}   zS = {zs:.6g}   (in input coordinates)")
    print(f"Iy  = {iy:.6g}")
    print(f"Iz  = {iz:.6g}")
    print(f"Iyz = {iyz:.6g}")
    print(f"I1  = {i1:.6g}   I2 = {i2:.6g}   phi* = {degrees(phi):.4g} deg")

    principal = abs(iyz) <= 1e-9 * max(abs(iy), abs(iz), 1.0)
    if not principal:
        print("\nWarning: Iyz != 0, so y and z are not principal axes.")
        print("Kern area and neutral axis below assume Iyz = 0 and are not valid.")

    kern = kern_area(pc, area, iy, iz)
    print("\n=== Kern area (centroidal) ===")
    for k, (y, z) in enumerate(kern, 1):
        print(f"K{k} = ({y:.4g}, {z:.4g})")

    na = neutral_axis(area, iy, iz, n, my, mz)
    a, b, c = na
    print("\n=== Neutral axis:  a*y + b*z + c = 0 ===")
    print(f"a = {a:.6g}   b = {b:.6g}   c = {c:.6g}")
    if abs(a) > EPS:
        print(f"crosses y-axis at y0 = {-c / a:.4g}")
    if abs(b) > EPS:
        print(f"crosses z-axis at z0 = {-c / b:.4g}")

    sig = [stress(y, z, area, iy, iz, n, my, mz) for y, z in pc]
    print(f"\nsigma_max = {max(sig):.4g}   sigma_min = {min(sig):.4g}")

    if plot or save:
        plot_section(pc, kern, na, phi, save=save)


def plot_section(pc, kern, na, phi, save=None):
    fig, ax = plt.subplots(figsize=(7, 6))

    sec = close_polygon(pc)
    ax.fill([p[0] for p in sec], [p[1] for p in sec], color="#cfd8e3", zorder=1)
    ax.plot([p[0] for p in sec], [p[1] for p in sec], color="#2b4a6f", lw=1.5,
            label="Cross-section", zorder=2)

    if len(kern) >= 3:
        k = close_polygon(kern)
        ax.fill([p[0] for p in k], [p[1] for p in k], color="#b08cd9", alpha=0.45, zorder=3)
        ax.plot([p[0] for p in k], [p[1] for p in k], color="#6a3d9a", lw=1.5,
                label="Kern area", zorder=4)

    ys = [p[0] for p in pc]
    zs = [p[1] for p in pc]
    span = max(max(ys) - min(ys), max(zs) - min(zs))
    pad = 0.25 * span
    ymin, ymax = min(ys) - pad, max(ys) + pad
    zmin, zmax = min(zs) - pad, max(zs) + pad

    # principal axes through the centroid
    length = 0.6 * span
    c, s = cos(phi), sin(phi)
    ax.plot([-length * c, length * c], [-length * s, length * s], "--", color="#c0392b",
            lw=1, label="Principal axis 1")
    ax.plot([length * s, -length * s], [-length * c, length * c], "--", color="#e67e22",
            lw=1, label="Principal axis 2")

    # neutral axis clipped to the plot window
    a, b, c0 = na
    if abs(b) > EPS:
        yy = [ymin, ymax]
        zz = [-(a * y + c0) / b for y in yy]
        ax.plot(yy, zz, color="#27ae60", lw=2, label="Neutral axis")
    elif abs(a) > EPS:
        ax.axvline(-c0 / a, color="#27ae60", lw=2, label="Neutral axis")

    ax.plot(0, 0, "k+", ms=10, mew=1.5, label="Centroid S")
    ax.set_xlim(ymin, ymax)
    ax.set_ylim(zmax, zmin)  # z points down
    ax.set_aspect("equal")
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax.set_title("Kern area and neutral axis (centroidal axes, z down)")
    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150)
        print(f"\nPlot saved to {save}")
    else:
        plt.show()


# ------------------------------------------------------------
# Example: T-section, units mm / N / Nmm
# ------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Kern area and neutral axis of a polygonal section")
    parser.add_argument("--save", help="save the plot to this file instead of showing it")
    args = parser.parse_args()

    t_section = [
        (-40, 0), (40, 0), (40, 12),   # flange
        (6, 12), (6, 70), (-6, 70),    # web
        (-6, 12), (-40, 12),
    ]

    analyse(
        t_section,
        n=-60e3,      # 60 kN compression
        my=1.2e6,     # 1.2 kNm
        mz=0.4e6,     # 0.4 kNm
        save=args.save,
    )
