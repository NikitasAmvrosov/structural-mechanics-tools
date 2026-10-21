from math import *
import ti_plotlib as plt

# ============================================================
# TI-Nspire CX II version of kern_area.py (same maths,
# matplotlib replaced by ti_plotlib). Edit the inputs below.
#
# z-down convention: ti_plotlib windows must have ymin < ymax,
# so there is no invert_yaxis(). Instead every point (y,z) is
# plotted at (y, -z). All plotting goes through _line/_pt which
# apply that flip once, so the math below stays in real (y,z).
#
# ROBUSTNESS: all numeric results are printed BEFORE plotting.
# The plot is wrapped in try/except -- if a ti_plotlib call
# behaves differently on your OS, you still get every number.
# ============================================================

# ---- consistent units: geometry mm, N in N, moments N*mm ----
# Points in order around the boundary (CW or CCW). (y, z).

# USE_CENTROID_SHIFT: True if the points are NOT given in centroidal axes.
# USE_MANUAL_SECTION: True to use A/Iy/Iz/Iyz below instead of the polygon
#   values (e.g. thin-walled profiles where the outline is only the hull).
USE_CENTROID_SHIFT = True
USE_MANUAL_SECTION = False

A_manual   = 0.0
Iy_manual  = 0.0
Iz_manual  = 0.0
Iyz_manual = 0.0

# example: T-section in mm (same as the PC version)
pts = [
    (-40, 0), (40, 0), (40, 12),
    (6, 12), (6, 70), (-6, 70),
    (-6, 12), (-40, 12),
]

N  = -60000     # N
My = 1200000    # Nmm
Mz = 400000     # Nmm

PLOT_SECTION        = True
PLOT_NEUTRAL_AXIS   = True
PLOT_KERNEL         = True
PLOT_PRINCIPAL_AXES = True

# screen is wider than tall; stretch x so units look ~square.
# tweak if the section looks squashed/stretched on the handheld.
ASPECT = 1.5

# ============================================================
# BASIC HELPERS
# ============================================================

EPS = 1e-9

def almost_zero(x):
    return abs(x) < EPS

def copy_points(P):
    return [(p[0], p[1]) for p in P]

def close_polygon(P):
    if len(P) == 0:
        return []
    if P[0] == P[-1]:
        return copy_points(P)
    Q = copy_points(P)
    Q.append(Q[0])
    return Q

def unique_points(P, tol=1e-7):
    out = []
    for p in P:
        ok = True
        for q in out:
            if abs(p[0]-q[0]) < tol and abs(p[1]-q[1]) < tol:
                ok = False
                break
        if ok:
            out.append(p)
    return out

def centroid_of_points(P):
    n = len(P)
    if n == 0:
        return (0.0, 0.0)
    sy = 0.0
    sz = 0.0
    for y, z in P:
        sy += y
        sz += z
    return (sy/n, sz/n)

def sort_by_angle(P):
    if len(P) <= 1:
        return copy_points(P)
    cy, cz = centroid_of_points(P)
    tmp = []
    for y, z in P:
        a = atan2(z - cz, y - cy)
        tmp.append((a, (y, z)))
    tmp.sort(key=lambda x: x[0])
    return [t[1] for t in tmp]

# ============================================================
# GEOMETRY: POLYGON AREA, CENTROID, INERTIA
# ============================================================

def polygon_area_centroid_inertia(P):
    if len(P) < 3:
        raise ValueError("Need at least 3 points.")

    Q = close_polygon(P)

    A2 = 0.0
    Cy_num = 0.0
    Cz_num = 0.0
    Iy0 = 0.0
    Iz0 = 0.0
    Iyz0 = 0.0

    for i in range(len(Q)-1):
        y1, z1 = Q[i]
        y2, z2 = Q[i+1]
        c = y1*z2 - y2*z1

        A2 += c
        Cy_num += (y1 + y2) * c
        Cz_num += (z1 + z2) * c

        Iy0  += (z1*z1 + z1*z2 + z2*z2) * c
        Iz0  += (y1*y1 + y1*y2 + y2*y2) * c
        Iyz0 += (2*y1*z1 + y1*z2 + y2*z1 + 2*y2*z2) * c

    A = 0.5 * A2

    if A < 0:
        Pr = copy_points(P)
        Pr.reverse()
        return polygon_area_centroid_inertia(Pr)

    yS = Cy_num / (6.0 * A)
    zS = Cz_num / (6.0 * A)

    Iy0  /= 12.0
    Iz0  /= 12.0
    Iyz0 /= 24.0

    IyS  = Iy0  - A * zS*zS
    IzS  = Iz0  - A * yS*yS
    IyzS = -(Iyz0 - A * yS*zS)   # minus convention: Iyz = -int(y*z dA)

    return A, yS, zS, IyS, IzS, IyzS

def shift_points(P, y0, z0):
    return [(y-y0, z-z0) for y, z in P]

# ============================================================
# PRINCIPAL AXES
# ============================================================

def principal_axes(Iy, Iz, Iyz):
    d = sqrt((Iy - Iz)*(Iy - Iz) + 4.0*Iyz*Iyz)
    I1 = 0.5*(Iy + Iz) + 0.5*d
    I2 = 0.5*(Iy + Iz) - 0.5*d
    if almost_zero(d):
        phi = 0.0
    else:
        phi = 0.5 * atan2(2.0*Iyz, Iy - Iz)
    return I1, I2, phi

# ============================================================
# CONVEX HULL (for kernel support lines)
# ============================================================

def cross(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def convex_hull(P):
    ptsL = unique_points(P)
    if len(ptsL) <= 2:
        return ptsL

    ptsL.sort()
    lower = []
    for p in ptsL:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(ptsL):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[:-1] + upper[:-1]

# ============================================================
# KERNEL AREA
# ============================================================

def line_from_points(p1, p2):
    y1, z1 = p1
    y2, z2 = p2
    a = z1 - z2
    b = y2 - y1
    c = y1*z2 - y2*z1
    return a, b, c

def kernel_point_from_support_line(a, b, c, iy2, iz2):
    if almost_zero(a):
        y0 = None
    else:
        y0 = -c / a
    if almost_zero(b):
        z0 = None
    else:
        z0 = -c / b

    if y0 is None:
        yK = 0.0
    elif almost_zero(y0):
        return None
    else:
        yK = -iz2 / y0

    if z0 is None:
        zK = 0.0
    elif almost_zero(z0):
        return None
    else:
        zK = -iy2 / z0

    return (yK, zK)

def compute_kernel_polygon(Pc, A, Iy, Iz):
    hull = convex_hull(Pc)
    if len(hull) < 3:
        return []

    iy2 = Iy / A
    iz2 = Iz / A

    K = []
    H = close_polygon(hull)
    for i in range(len(H)-1):
        a, b, c = line_from_points(H[i], H[i+1])
        kp = kernel_point_from_support_line(a, b, c, iy2, iz2)
        if kp is not None:
            K.append(kp)

    K = unique_points(K)
    K = sort_by_angle(K)
    return K

# ============================================================
# NEUTRAL AXIS  (Formelsammlung 3.2, gilt nur fuer Iyz = 0)
# ============================================================

def neutral_axis_coeffs(A, Iy, Iz, N, My, Mz):
    aa = -Mz / Iz if not almost_zero(Iz) else 0.0
    bb =  My / Iy if not almost_zero(Iy) else 0.0
    cc =  N / A
    return aa, bb, cc

def line_box_intersections(a, b, c, ymin, ymax, zmin, zmax):
    cand = []
    if not almost_zero(b):
        z = -(a*ymin + c)/b
        if zmin - 1e-8 <= z <= zmax + 1e-8:
            cand.append((ymin, z))
        z = -(a*ymax + c)/b
        if zmin - 1e-8 <= z <= zmax + 1e-8:
            cand.append((ymax, z))
    if not almost_zero(a):
        y = -(b*zmin + c)/a
        if ymin - 1e-8 <= y <= ymax + 1e-8:
            cand.append((y, zmin))
        y = -(b*zmax + c)/a
        if ymin - 1e-8 <= y <= ymax + 1e-8:
            cand.append((y, zmax))
    cand = unique_points(cand)
    if len(cand) >= 2:
        return cand[0], cand[1]
    return None

# ============================================================
# PLOTTING (ti_plotlib) -- all (y,z) flipped to (y,-z) here
# ============================================================

def _line(y1, z1, y2, z2):
    plt.line(y1, -z1, y2, -z2, "")

def _pt(y, z, mark):
    plt.plot(y, -z, mark)

def extents(P):
    ys = [p[0] for p in P]
    zs = [p[1] for p in P]
    return min(ys), max(ys), min(zs), max(zs)

def plot_all(Pc, Kc, na_coeffs, phi):
    all_pts = copy_points(Pc) + list(Kc)
    ymin, ymax, zmin, zmax = extents(all_pts)
    dy = max(ymax - ymin, 1.0)
    dz = max(zmax - zmin, 1.0)
    span = max(dy, dz)
    half = 0.5 * span + 0.25 * span
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)

    xlo = cy - ASPECT * half
    xhi = cy + ASPECT * half
    # vertical axis represents -z
    vlo = -cz - half
    vhi = -cz + half

    plt.cls()
    plt.window(xlo, xhi, vlo, vhi)

    # axes-mode string + grid are the OS-fragile calls; guard them
    try:
        plt.axes("on")
    except:
        pass

    # origin cross
    plt.color(0, 0, 0)
    _line(xlo, 0, xhi, 0)
    _line(0, -vhi, 0, -vlo)

    if PLOT_SECTION:
        plt.color(0, 0, 220)
        Q = close_polygon(Pc)
        for i in range(len(Q)-1):
            _line(Q[i][0], Q[i][1], Q[i+1][0], Q[i+1][1])
        for p in Pc:
            _pt(p[0], p[1], "o")

    if PLOT_PRINCIPAL_AXES:
        L = 0.9 * half
        c1, s1 = cos(phi), sin(phi)
        plt.color(220, 0, 0)
        _line(-L*c1, -L*s1, L*c1, L*s1)   # I1
        plt.color(240, 140, 0)
        _line(L*s1, -L*c1, -L*s1, L*c1)   # I2 (perpendicular)

    if PLOT_NEUTRAL_AXIS:
        a, b, c_ = na_coeffs
        seg = line_box_intersections(a, b, c_, xlo, xhi, -vhi, -vlo)
        if seg is not None:
            plt.color(0, 170, 0)
            _line(seg[0][0], seg[0][1], seg[1][0], seg[1][1])

    if PLOT_KERNEL and len(Kc) >= 3:
        plt.color(150, 0, 180)
        K = close_polygon(Kc)
        for i in range(len(K)-1):
            _line(K[i][0], K[i][1], K[i+1][0], K[i+1][1])
        for p in Kc:
            _pt(p[0], p[1], "x")

    try:
        plt.title("blau=QS gruen=nF lila=Kern rot=HA")
    except:
        pass

    plt.show_plot()

# ============================================================
# OUTPUT
# ============================================================

def fmt(x):
    return str(round(x, 6))

def print_results(A, yS, zS, Iy, Iz, Iyz, I1, I2, phi, na):
    print("=== Querschnittswerte ===")
    print("A   =", fmt(A))
    print("yS  =", fmt(yS))
    print("zS  =", fmt(zS))
    print("Iy  =", fmt(Iy))
    print("Iz  =", fmt(Iz))
    print("Iyz =", fmt(Iyz))

    if not almost_zero(Iyz):
        print("!! Iyz != 0: neutrale Faser via 3.2 UNGUELTIG,")
        print("   3.3 (beliebige Achsen) noetig.")

    print("=== Hauptachsen ===")
    print("I1 =", fmt(I1))
    print("I2 =", fmt(I2))
    print("phi* rad =", fmt(phi))
    print("phi* deg =", fmt(degrees(phi)))

    print("=== Biegenormalspannung ===")
    a, b, c = na
    print("a*y + b*z + c = 0")
    print("a =", fmt(a))
    print("b =", fmt(b))
    print("c =", fmt(c))
    if not almost_zero(b):
        print("z(y) =", fmt(-a/b), "* y +", fmt(-c/b))
    elif not almost_zero(a):
        print("vertikal: y =", fmt(-c/a))
    else:
        print("keine eindeutige neutrale Faser")

    print("--- Schnittpunkte ---")
    if not almost_zero(a):
        print("y-Achse (z=0): y0 =", fmt(-c/a))
    else:
        print("y-Achse: kein Schnitt")
    if not almost_zero(b):
        print("z-Achse (y=0): z0 =", fmt(-c/b))
    else:
        print("z-Achse: kein Schnitt")

def print_kernel(Kc, yS, zS):
    print("=== Kernflaeche (Schwerpunktachsen) ===")
    if len(Kc) < 3:
        print("keine gueltige Kernflaeche.")
        return
    for i in range(len(Kc)):
        y, z = Kc[i]
        print("K", i+1, "= (", fmt(y), ",", fmt(z), ")")
    print("=== Kernflaeche (global) ===")
    for i in range(len(Kc)):
        y, z = Kc[i]
        print("Kg", i+1, "= (", fmt(y+yS), ",", fmt(z+zS), ")")

# ============================================================
# MAIN
# ============================================================

def run():
    A_calc, yS, zS, Iy_calc, Iz_calc, Iyz_calc = \
        polygon_area_centroid_inertia(pts)

    if USE_CENTROID_SHIFT:
        Pc = shift_points(pts, yS, zS)
    else:
        Pc = copy_points(pts)
        yS = 0.0
        zS = 0.0

    if USE_MANUAL_SECTION:
        A, Iy, Iz, Iyz = A_manual, Iy_manual, Iz_manual, Iyz_manual
    else:
        A, Iy, Iz, Iyz = A_calc, Iy_calc, Iz_calc, Iyz_calc

    I1, I2, phi = principal_axes(Iy, Iz, Iyz)
    na = neutral_axis_coeffs(A, Iy, Iz, N, My, Mz)
    Kc = compute_kernel_polygon(Pc, A, Iy, Iz)

    print_results(A, yS, zS, Iy, Iz, Iyz, I1, I2, phi, na)
    print_kernel(Kc, yS, zS)

    if PLOT_SECTION:
        try:
            plot_all(Pc, Kc, na, phi)
        except Exception as e:
            print("Plot uebersprungen:", e)

run()

