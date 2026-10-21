from math import pi, sqrt, atan2, sin, cos

# =========================================================
#   COMPOSITE SECTION CALCULATOR
#   Area, centroid, second moments of area (Iy, Iz, Iyz) and
#   principal axes of a cross-section built from standard shapes.
#   Author: Nikitas Amvrosov
# ---------------------------------------------------------
#   CONTROLS (at every prompt):
#     xx -> re-enter the current part / prompt
#     q  -> quit (the calculation is saved on the calculator)
#
#   REVIEW LOOP:  -2 cycles  short -> detail -> STEPS
#     STEPS = full worked solution line by line
#     (formula -> numbers -> result), to check hand calculations.
#
#   CONVENTION (z points down, minus sign for Iyz):
#     Iyz = -integral(y*z dA),  Steiner: Iyz_st = -A*dy*dz
#
#   THIN-WALLED: 'duenn stab' neglects t^3, 'rechteck' keeps it.
#   PARAMETER (a): set at the start, then type '2*a', 'a/2', ...
#   The prompts are in German.
# =========================================================

SYM = {}
UNIT = ""

# --- only used for the step-by-step output, not for the results ---
TXT = {}      # id(teil-tupel) -> liste von textzeilen
LASTTXT = []  # von der zuletzt eingegebenen form gefuellt

try:
    import ti_system as _tsys
    HAVE_STORE = True
except:
    HAVE_STORE = False


class Zurueck(Exception):
    pass


class Ende(Exception):
    pass


def F(x):
    # zahl kurz formatiert (5 signifikante stellen)
    try:
        if x == 0:
            return "0"
        return "{:.5g}".format(x)
    except:
        return str(round(x, 4))


def P(x):
    # negative zahl in klammern (fuer produkte)
    if x < 0:
        return "(" + F(x) + ")"
    return F(x)


def T(s):
    # zeile fuer handrechnung merken
    LASTTXT.append(s)


# -----------------------------------------
#   SPEICHERN / LADEN (defensiv; fehler = kein persist)
# -----------------------------------------

def save_state(parts, modus, neg):
    if not HAVE_STORE:
        return False
    try:
        flat = []
        for p in parts:
            for v in p:
                flat.append(float(v))
        if len(flat) == 0:
            return False
        _tsys.store_list("ftmflat", flat)
        _tsys.store_list("ftmmeta", [float(modus), float(neg), float(len(parts))])
        return True
    except:
        return False


def load_state():
    if not HAVE_STORE:
        return None
    try:
        flat = _tsys.recall_list("ftmflat")
        meta = _tsys.recall_list("ftmmeta")
    except:
        return None
    try:
        modus = int(meta[0])
        neg = int(meta[1])
        n = int(meta[2])
        if len(flat) < 6 * n or n < 1:
            return None
        parts = []
        for k in range(n):
            b = 6 * k
            parts.append((flat[b], flat[b+1], flat[b+2],
                          flat[b+3], flat[b+4], flat[b+5]))
        return parts, modus, neg
    except:
        return None


# -----------------------------------------
#   EINGABE-HELFER
# -----------------------------------------

def num(prompt):
    while True:
        s = input(prompt).strip()
        if s == "q":
            raise Ende()
        if s == "xx":
            raise Zurueck()
        try:
            return float(s)
        except ValueError:
            pass
        try:
            return float(eval(s, dict(SYM)))
        except:
            print("ungueltig  (xx=neu, q=beenden)")


def intval(prompt):
    while True:
        s = input(prompt).strip()
        if s == "q":
            raise Ende()
        if s == "xx":
            raise Zurueck()
        try:
            return int(s)
        except ValueError:
            print("ganze zahl  (xx=neu, q=beenden)")


def ask_int(prompt):
    while True:
        try:
            return intval(prompt)
        except Zurueck:
            print("xx hier ohne wirkung  (q=beenden)")


def lab(exp):
    if UNIT == "":
        return ""
    return "* " + UNIT + "^" + str(exp)


# -----------------------------------------
#   KERNFUNKTIONEN
# -----------------------------------------

def vorzeichen_teil(p, neg_flag):
    if neg_flag == 1:
        return (-p[0], p[1], p[2], -p[3], -p[4], -p[5])
    return p


def rotate_ftm(iyy, izz, iyz, phi_deg):
    phi = phi_deg * pi / 180
    iyy_r = 0.5 * (iyy + izz) + 0.5 * (iyy - izz) * cos(2 * phi) + iyz * sin(2 * phi)
    izz_r = 0.5 * (iyy + izz) - 0.5 * (iyy - izz) * cos(2 * phi) - iyz * sin(2 * phi)
    iyz_r = -0.5 * (iyy - izz) * sin(2 * phi) + iyz * cos(2 * phi)
    return iyy_r, izz_r, iyz_r


def steineranteil(a, y, z, ys, zs):
    dy = y - ys
    dz = z - zs
    return a * dz**2, a * dy**2, -a * dy * dz


def frag_drehung(iyy, izz, iyz):
    if intval("geneigt? 1=ja 0=nein: ") == 1:
        phi = num("phi_deg=")
        r = rotate_ftm(iyy, izz, iyz, phi)
        T("gedreht (FS 1.3), phi=" + F(phi) + " deg:")
        T(" Ieta=(Iy+Iz)/2+(Iy-Iz)/2*cos2phi+Iyz*sin2phi")
        T("  =(" + F(iyy) + "+" + F(izz) + ")/2+(" + F(iyy) + "-" + F(izz)
          + ")/2*cos(" + F(2 * phi) + ")+" + F(iyz) + "*sin(" + F(2 * phi) + ")")
        T("  =" + F(r[0]))
        T(" Izeta=(Iy+Iz)/2-(Iy-Iz)/2*cos2phi-Iyz*sin2phi =" + F(r[1]))
        T(" Ietazeta=-(Iy-Iz)/2*sin2phi+Iyz*cos2phi =" + F(r[2]))
        return r
    return iyy, izz, iyz


# --- formen: liefern (A, y_s, z_s, Iyy, Izz, Iyz) ---

def rechteck():
    print("rechteck (mit t^3-term)")
    b = num("b=")
    h = num("h=")
    y = num("y_i=")
    z = num("z_i=")
    iyy = b * h**3 / 12
    izz = h * b**3 / 12
    T("rechteck b=" + F(b) + " h=" + F(h))
    T(" A=b*h=" + F(b) + "*" + F(h) + "=" + F(b * h))
    T(" Iy_eig=b*h^3/12=" + F(b) + "*" + F(h) + "^3/12=" + F(iyy))
    T(" Iz_eig=h*b^3/12=" + F(h) + "*" + F(b) + "^3/12=" + F(izz))
    T(" Iyz_eig=0 (symmetrie)")
    iyy, izz, iyz = frag_drehung(iyy, izz, 0.0)
    return b * h, y, z, iyy, izz, iyz


def vollkreis():
    print("vollkreis")
    r = num("r=")
    y = num("y_i=")
    z = num("z_i=")
    T("vollkreis r=" + F(r))
    T(" A=pi*r^2=pi*" + F(r) + "^2=" + F(pi * r**2))
    T(" Iy_eig=Iz_eig=pi*r^4/4=pi*" + F(r) + "^4/4=" + F(pi * r**4 / 4))
    T(" Iyz_eig=0")
    return pi * r**2, y, z, pi * r**4 / 4, pi * r**4 / 4, 0.0


def kreisring():
    print("kreisring")
    ra = num("ra=")
    ri = num("ri=")
    y = num("y_i=")
    z = num("z_i=")
    i = pi * (ra**4 - ri**4) / 4
    T("kreisring ra=" + F(ra) + " ri=" + F(ri))
    T(" A=pi*(ra^2-ri^2)=pi*(" + F(ra) + "^2-" + F(ri) + "^2)="
      + F(pi * (ra**2 - ri**2)))
    T(" Iy_eig=Iz_eig=pi*(ra^4-ri^4)/4=pi*(" + F(ra) + "^4-" + F(ri)
      + "^4)/4=" + F(i))
    T(" Iyz_eig=0")
    return pi * (ra**2 - ri**2), y, z, i, i, 0.0


def duenner_kreisring():
    print("duenner kreisring (t << Rm)")
    rm = num("Rm=")
    t = num("t=")
    y = num("y_i=")
    z = num("z_i=")
    T("duenner kreisring Rm=" + F(rm) + " t=" + F(t))
    T(" A=2*pi*Rm*t=2*pi*" + F(rm) + "*" + F(t) + "=" + F(2 * pi * rm * t))
    T(" Iy_eig=Iz_eig=pi*Rm^3*t=pi*" + F(rm) + "^3*" + F(t) + "="
      + F(pi * rm**3 * t))
    T(" Iyz_eig=0")
    return 2 * pi * rm * t, y, z, pi * rm**3 * t, pi * rm**3 * t, 0.0


def rechtw_dreieck():
    print("rechtw dreieck (rechter winkel unten links, schenkel +y und -z)")
    b = num("b=")
    h = num("h=")
    y = num("y_i=")
    z = num("z_i=")
    iyy = b * h**3 / 36
    izz = h * b**3 / 36
    iyz = -(b**2 * h**2) / 72
    T("rechtw dreieck b=" + F(b) + " h=" + F(h))
    T(" A=b*h/2=" + F(b) + "*" + F(h) + "/2=" + F(b * h / 2))
    T(" Iy_eig=b*h^3/36=" + F(b) + "*" + F(h) + "^3/36=" + F(iyy))
    T(" Iz_eig=h*b^3/36=" + F(h) + "*" + F(b) + "^3/36=" + F(izz))
    T(" Iyz_eig=-b^2*h^2/72=-" + F(b) + "^2*" + F(h) + "^2/72=" + F(iyz))
    iyy, izz, iyz = frag_drehung(iyy, izz, iyz)
    return b * h / 2, y, z, iyy, izz, iyz


def dreieck_allgemein():
    print("dreieck allgemein")
    b = num("b=")
    h = num("h=")
    a_shift = num("a_shift=")
    y = num("y_i=")
    z = num("z_i=")
    iyy = b * h**3 / 36
    izz = h * (b**3 - b**2 * a_shift + b * a_shift**2) / 36
    iyz = -b * h**2 * (b - 2 * a_shift) / 72
    T("dreieck b=" + F(b) + " h=" + F(h) + " a=" + F(a_shift))
    T(" A=b*h/2=" + F(b) + "*" + F(h) + "/2=" + F(b * h / 2))
    T(" Iy_eig=b*h^3/36=" + F(b) + "*" + F(h) + "^3/36=" + F(iyy))
    T(" Iz_eig=b*h/36*(b^2-b*a+a^2)")
    T("  =" + F(b) + "*" + F(h) + "/36*(" + F(b) + "^2-" + F(b) + "*"
      + F(a_shift) + "+" + F(a_shift) + "^2)=" + F(izz))
    T(" Iyz_eig=-b*h^2/72*(b-2a)")
    T("  =-" + F(b) + "*" + F(h) + "^2/72*(" + F(b) + "-2*" + F(a_shift)
      + ")=" + F(iyz))
    iyy, izz, iyz = frag_drehung(iyy, izz, iyz)
    return b * h / 2, y, z, iyy, izz, iyz


def halbkreis():
    print("halbkreis (gerade kante VERTIKAL, y = symmetrieachse)")
    r = num("r=")
    y = num("y_i=")
    z = num("z_i=")
    iyy = pi * r**4 / 8
    izz = r**4 * (9 * pi**2 - 64) / (72 * pi)
    T("halbkreis r=" + F(r) + " (kante vertikal)")
    T(" A=pi*r^2/2=pi*" + F(r) + "^2/2=" + F(pi * r**2 / 2))
    T(" Iy_eig=pi*r^4/8=pi*" + F(r) + "^4/8=" + F(iyy))
    T(" Iz_eig=r^4/(72pi)*(9pi^2-64)=" + F(r) + "^4/(72pi)*(9pi^2-64)="
      + F(izz))
    T(" Iyz_eig=0 (symmetrie)")
    iyy, izz, iyz = frag_drehung(iyy, izz, 0.0)
    return pi * r**2 / 2, y, z, iyy, izz, iyz


def viertelkreis():
    print("viertelkreis (ecke = punkt mit rechtem winkel)")
    r = num("r=")
    yc = num("y_ecke=")
    zc = num("z_ecke=")
    print("richtung ab ecke (+1/-1):")
    sy = num("sy=")
    sz = num("sz=")
    e = 4 * r / (3 * pi)
    y = yc + sy * e
    z = zc + sz * e
    i = pi * r**4 / 16 - 4 * r**4 / (9 * pi)
    iyz = -sy * sz * r**4 * (0.125 - 4 / (9 * pi))
    print("  -> schwerpunkt y=", round(y, 4), " z=", round(z, 4))
    T("viertelkreis r=" + F(r))
    T(" A=pi*r^2/4=pi*" + F(r) + "^2/4=" + F(pi * r**2 / 4))
    T(" e=4r/(3pi)=4*" + F(r) + "/(3pi)=" + F(e))
    T(" y_i=" + F(yc) + "+(" + F(sy) + ")*" + F(e) + "=" + F(y))
    T(" z_i=" + F(zc) + "+(" + F(sz) + ")*" + F(e) + "=" + F(z))
    T(" (eckachsen: I=pi*r^4/16, Iyz=-r^4/8*sy*sz)")
    T(" Iy_eig=Iz_eig=pi*r^4/16-A*e^2")
    T("  =pi*" + F(r) + "^4/16-4*" + F(r) + "^4/(9pi)=" + F(i))
    T(" Iyz_eig=-sy*sz*r^4*(1/8-4/(9pi))=" + F(iyz))
    return pi * r**2 / 4, y, z, i, i, iyz


def duenn_stab():
    print("duenn stab (ohne t^3-term)  1=horizontal 2=vertikal")
    typ = intval("typ=")
    l = num("l=")
    t = num("t=")
    y = num("y_i=")
    z = num("z_i=")
    if typ == 1:
        iyy, izz = 0.0, t * l**3 / 12
        T("duenn stab HORIZONTAL l=" + F(l) + " t=" + F(t))
        T(" A=l*t=" + F(l) + "*" + F(t) + "=" + F(l * t))
        T(" Iy_eig=0 (t^3 vernachl.)")
        T(" Iz_eig=t*l^3/12=" + F(t) + "*" + F(l) + "^3/12=" + F(izz))
    elif typ == 2:
        iyy, izz = t * l**3 / 12, 0.0
        T("duenn stab VERTIKAL l=" + F(l) + " t=" + F(t))
        T(" A=l*t=" + F(l) + "*" + F(t) + "=" + F(l * t))
        T(" Iy_eig=t*l^3/12=" + F(t) + "*" + F(l) + "^3/12=" + F(iyy))
        T(" Iz_eig=0 (t^3 vernachl.)")
    else:
        print("ungueltiger typ")
        raise Zurueck()
    T(" Iyz_eig=0")
    iyy, izz, iyz = frag_drehung(iyy, izz, 0.0)
    return l * t, y, z, iyy, izz, iyz


def eingabe_teil():
    global LASTTXT
    print()
    print("1 rechteck          5 halbkreis")
    print("2 vollkreis         6 duenn stab")
    print("3 kreisring         7 dreieck allg")
    print("4 rechtw dreieck    8 viertelkreis")
    print("                    9 duenner kreisring")
    w = intval("wahl=")
    formen = {1: rechteck, 2: vollkreis, 3: kreisring, 4: rechtw_dreieck,
              5: halbkreis, 6: duenn_stab, 7: dreieck_allgemein,
              8: viertelkreis, 9: duenner_kreisring}
    if w in formen:
        LASTTXT = []
        p = formen[w]()
        try:
            TXT[id(p)] = LASTTXT
        except:
            pass
        return p
    print("ungueltige wahl")
    raise Zurueck()


def neuer_teil(neg_global):
    while True:
        try:
            p = eingabe_teil()
            if neg_global == 1:
                p0 = p
                p = vorzeichen_teil(p, intval("negativ? 1=ja 0=nein: "))
                if p is not p0:
                    try:
                        TXT[id(p)] = TXT.get(id(p0), []) + [
                            " NEGATIV (aussparung): alle werte mit minus"]
                    except:
                        pass
            return p
        except Zurueck:
            print("-- teil neu --")


# -----------------------------------------
#   BERECHNUNG (rein)
# -----------------------------------------

def rechne(parts, modus, detail):
    a_ges = 0.0
    for p in parts:
        a_ges += p[0]
    if a_ges == 0:
        print("fehler: a_ges = 0")
        return None

    if modus == 1:
        sy = sz = 0.0
        for p in parts:
            sy += p[0] * p[1]
            sz += p[0] * p[2]
        ys, zs = sy / a_ges, sz / a_ges
        print()
        print("schwerpunkt gesamt")
        print("A_ges =", a_ges, lab(2))
        print("ys =", ys, lab(1))
        print("zs =", zs, lab(1))
    else:
        ys = zs = 0.0
        print()
        print("bezug: gegebene achsen (ursprung y=0 z=0)")
        print("A_ges =", a_ges, lab(2))

    iy = iz = iyz = 0.0
    for k in range(len(parts)):
        p = parts[k]
        iy_st, iz_st, iyz_st = steineranteil(p[0], p[1], p[2], ys, zs)
        iy_k = p[3] + iy_st
        iz_k = p[4] + iz_st
        iyz_k = p[5] + iyz_st
        iy += iy_k
        iz += iz_k
        iyz += iyz_k
        if detail == 1:
            print()
            print("teil", k + 1)
            print("  A_i  =", p[0])
            print("  y_i  =", p[1], "  z_i =", p[2])
            print("  Iy_i =", p[3], " Iz_i =", p[4], " Iyz_i =", p[5])
            print("  Iy_st=", iy_st, " Iz_st=", iz_st, " Iyz_st=", iyz_st)
            print("  -> Iy=", iy_k, " Iz=", iz_k, " Iyz=", iyz_k)
        else:
            print("teil", k + 1, ": Iy=", round(iy_k, 3),
                  " Iz=", round(iz_k, 3), " Iyz=", round(iyz_k, 3))

    print()
    print("ergebnisse")
    print("Iy  =", iy, lab(4))
    print("Iz  =", iz, lab(4))
    print("Iyz =", iyz, lab(4))
    print("Ip  =", iy + iz, lab(4))

    skala = max(abs(iy), abs(iz))
    eps = 1e-9 * skala if skala != 0 else 1e-12
    if abs(iyz) <= eps:
        print("Iyz ~ 0 -> y,z sind bereits hauptachsen")

    if modus == 1:
        wurzel = sqrt((iy - iz)**2 + 4 * iyz**2)
        i1 = 0.5 * (iy + iz) + 0.5 * wurzel
        i2 = 0.5 * (iy + iz) - 0.5 * wurzel
        if abs(iy - iz) <= eps and abs(iyz) <= eps:
            print("I1 =", iy, lab(4))
            print("I2 =", iz, lab(4))
            print("phi_deg = beliebig (I in jeder richtung gleich)")
        else:
            phi = 0.5 * atan2(2 * iyz, iy - iz)
            print("I1 =", i1, lab(4))
            print("I2 =", i2, lab(4))
            print("phi_deg =", phi * 180 / pi, " (richtung von I1/Imax)")
    else:
        print("(I1/I2/phi nur im schwerpunkt -> modus 1)")

    return iy, iz, iyz


# -----------------------------------------
#   STEP-BY-STEP OUTPUT (recomputes everything independently,
#   so it must match 'ergebnisse' -> built-in cross-check)
# -----------------------------------------

def summe_txt(werte):
    s = ""
    for k in range(len(werte)):
        v = werte[k]
        if k == 0:
            s = F(v)
        elif v < 0:
            s = s + "-" + F(-v)
        else:
            s = s + "+" + F(v)
    return s


def hand(parts, modus):
    n = len(parts)
    print()
    print("========== HANDRECHNUNG ==========")
    print("konvention: Iyz=-int(y*z dA), Iyz_st=-A*dy*dz")
    print("u = einheit wie eingegeben", lab(1))
    a_ges = 0.0
    aw = []
    for p in parts:
        a_ges += p[0]
        aw.append(p[0])
    if a_ges == 0:
        print("fehler: a_ges = 0")
        return

    # --- 1) teile einzeln (formel -> zahlen) ---
    print()
    print("--- 1) TEILFLAECHEN ---")
    for k in range(n):
        p = parts[k]
        print()
        print("teil", k + 1, ":")
        zeilen = TXT.get(id(p), None)
        if zeilen:
            for z in zeilen:
                print(z)
        else:
            print(" (formel-text n.v., nur werte)")
            print(" A=" + F(p[0]))
            print(" Iy_eig=" + F(p[3]) + " Iz_eig=" + F(p[4])
                  + " Iyz_eig=" + F(p[5]))
        print(" lage: y_i=" + F(p[1]) + "  z_i=" + F(p[2]))

    # --- 2) schwerpunkt ---
    print()
    print("--- 2) SCHWERPUNKT ---")
    print("A_ges=sum A_i=" + summe_txt(aw) + "=" + F(a_ges))
    if modus == 1:
        ay = []
        az = []
        ty = ""
        tz = ""
        for k in range(n):
            p = parts[k]
            ay.append(p[0] * p[1])
            az.append(p[0] * p[2])
            if k > 0:
                ty += "+"
                tz += "+"
            ty += P(p[0]) + "*" + P(p[1])
            tz += P(p[0]) + "*" + P(p[2])
        ys = 0.0
        zs = 0.0
        for v in ay:
            ys += v
        for v in az:
            zs += v
        ys = ys / a_ges
        zs = zs / a_ges
        print("ys=sum(A_i*y_i)/A_ges")
        print(" =(" + ty + ")/" + F(a_ges))
        print(" =(" + summe_txt(ay) + ")/" + F(a_ges) + "=" + F(ys))
        print("zs=sum(A_i*z_i)/A_ges")
        print(" =(" + tz + ")/" + F(a_ges))
        print(" =(" + summe_txt(az) + ")/" + F(a_ges) + "=" + F(zs))
    else:
        ys = 0.0
        zs = 0.0
        print("gegebene achsen -> bezug ys=0, zs=0")

    # --- 3) steiner pro teil ---
    print()
    print("--- 3) STEINER (FS 1.2) ---")
    liy = []
    liz = []
    liyz = []
    for k in range(n):
        p = parts[k]
        dy = p[1] - ys
        dz = p[2] - zs
        iy_k = p[3] + p[0] * dz**2
        iz_k = p[4] + p[0] * dy**2
        iyz_k = p[5] - p[0] * dy * dz
        liy.append(iy_k)
        liz.append(iz_k)
        liyz.append(iyz_k)
        print()
        print("teil", k + 1, ":")
        print(" dy=y_i-ys=" + F(p[1]) + "-(" + F(ys) + ")=" + F(dy))
        print(" dz=z_i-zs=" + F(p[2]) + "-(" + F(zs) + ")=" + F(dz))
        print(" Iy=Iy_eig+A*dz^2")
        print("  =" + F(p[3]) + "+(" + F(p[0]) + ")*(" + F(dz) + ")^2")
        print("  =" + F(p[3]) + "+(" + F(p[0] * dz**2) + ")=" + F(iy_k))
        print(" Iz=Iz_eig+A*dy^2")
        print("  =" + F(p[4]) + "+(" + F(p[0]) + ")*(" + F(dy) + ")^2")
        print("  =" + F(p[4]) + "+(" + F(p[0] * dy**2) + ")=" + F(iz_k))
        print(" Iyz=Iyz_eig-A*dy*dz")
        print("  =" + F(p[5]) + "-(" + F(p[0]) + ")*(" + F(dy) + ")*("
              + F(dz) + ")")
        print("  =" + F(p[5]) + "+(" + F(-p[0] * dy * dz) + ")=" + F(iyz_k))

    # --- 4) summe ---
    iy = 0.0
    iz = 0.0
    iyz = 0.0
    for v in liy:
        iy += v
    for v in liz:
        iz += v
    for v in liyz:
        iyz += v
    print()
    print("--- 4) GESAMT ---")
    print("Iy=" + summe_txt(liy) + "=" + F(iy), lab(4))
    print("Iz=" + summe_txt(liz) + "=" + F(iz), lab(4))
    print("Iyz=" + summe_txt(liyz) + "=" + F(iyz), lab(4))
    print("Ip=Iy+Iz=" + F(iy) + "+" + F(iz) + "=" + F(iy + iz), lab(4))

    # --- 5) hauptachsen ---
    if modus == 1:
        print()
        print("--- 5) HAUPTACHSEN (FS 1.5) ---")
        w = sqrt((iy - iz)**2 + 4 * iyz**2)
        print("I1,2=(Iy+Iz)/2 +- 1/2*sqrt((Iy-Iz)^2+4Iyz^2)")
        print(" =(" + F(iy) + "+" + F(iz) + ")/2 +- 1/2*sqrt(("
              + F(iy) + "-" + F(iz) + ")^2+4*(" + F(iyz) + ")^2)")
        print(" =" + F(0.5 * (iy + iz)) + " +- " + F(0.5 * w))
        print("I1=" + F(0.5 * (iy + iz) + 0.5 * w)
              + "  I2=" + F(0.5 * (iy + iz) - 0.5 * w), lab(4))
        if abs(iy - iz) < 1e-12 * max(1.0, abs(iy)) and abs(iyz) < 1e-12 * max(1.0, abs(iy)):
            print("phi beliebig (Iy=Iz, Iyz=0)")
        else:
            phi = 0.5 * atan2(2 * iyz, iy - iz) * 180 / pi
            print("tan(2phi*)=2Iyz/(Iy-Iz)=2*(" + F(iyz) + ")/("
                  + F(iy) + "-" + F(iz) + ")")
            print(" -> phi*=" + F(phi) + " deg (achse von I1)")
    print("==================================")
    print("(muss mit 'ergebnisse' oben uebereinstimmen)")


def entfernen(parts, modus, detail):
    # Gesamt-FTM MIT einem entfernten teilstueck.
    # aendert 'parts' NICHT -> volltotal (aufgabe a) bleibt danach erhalten.
    print()
    print("gesamt OHNE ein teilstueck")
    print(" 1 = ganzes vorhandenes teil (nummer) entfernen")
    print(" 2 = teilstueck abziehen (z.b. ausschnitt A aus einem teil)")
    w = ask_int("wahl: ")
    if w == 1:
        k = ask_int("welche teil-nr entfernen: ")
        if not (1 <= k <= len(parts)):
            print("gibt es nicht")
            return
        rest = []
        for i in range(len(parts)):
            if i != k - 1:
                rest.append(parts[i])
        if len(rest) == 0:
            print("nichts mehr uebrig")
            return
        print()
        print("=== GESAMT OHNE TEIL", k, "===")
        rechne(rest, modus, detail)
    elif w == 2:
        print("teilstueck (A) eingeben, das ABGEZOGEN wird:")
        try:
            s = eingabe_teil()
        except Zurueck:
            print("abgebrochen")
            return
        neg = (-s[0], s[1], s[2], -s[3], -s[4], -s[5])
        print("  abzug: A=", round(s[0], 4),
              " y=", round(s[1], 3), " z=", round(s[2], 3))
        print()
        print("=== GESAMT OHNE DIESES TEILSTUECK ===")
        # in modus 2 zieht der negative teil seinen anteil um die
        # gegebenen achsen ab; in modus 1 wird der restschwerpunkt neu bestimmt.
        rechne(parts + [neg], modus, detail)
    else:
        print("ungueltig")


def zeige_parts(parts):
    print()
    print("teile:")
    for k in range(len(parts)):
        p = parts[k]
        print(" ", k + 1, "A=", round(p[0], 3),
              " y=", round(p[1], 3), " z=", round(p[2], 3))


# -----------------------------------------
#   KORREKTUR-SCHLEIFE (von frisch UND von reload benutzt)
# -----------------------------------------

def review(parts, modus, neg_global):
    detail = 0
    while True:
        save_state(parts, modus, neg_global)   # immer neueste version sichern
        rechne(parts, modus, 1 if detail == 1 else 0)
        if detail == 2:
            try:
                hand(parts, modus)
            except Exception as e:
                print("handrechnung-fehler:", e)
                print("(ergebnisse oben sind trotzdem gueltig)")
        zeige_parts(parts)
        print()
        print("modus anzeige:", ["kurz", "detail", "HAND"][detail])
        print("0=fertig  -1=teil anhaengen  -2=kurz/detail/HAND")
        print("-3=gesamt ohne teil  (q=beenden)")
        k = ask_int("korrigieren? teil-nr: ")
        if k == 0:
            break
        elif k == -2:
            detail = (detail + 1) % 3
        elif k == -3:
            entfernen(parts, modus, 1 if detail == 1 else 0)
        elif k == -1:
            parts.append(neuer_teil(neg_global))
        elif 1 <= k <= len(parts):
            print("teil", k, "neu eingeben")
            parts[k - 1] = neuer_teil(neg_global)
        else:
            print("gibt es nicht")
    ok = save_state(parts, modus, neg_global)
    if ok:
        print("fertig. gespeichert -> beim naechsten start ladbar")
    else:
        print("fertig. (speichern nicht verfuegbar / fehlgeschlagen)")


def ftm():
    global UNIT, SYM
    SYM = {}
    UNIT = ""
    print("steuerung:  xx = teil/prompt neu   q = beenden")

    s = input("parameter name (leer=numerisch): ").strip()
    if s == "q":
        raise Ende()
    if s != "" and s != "xx":
        SYM[s] = 1.0
        UNIT = s
        print("rechne mit", s, "= 1  -> ergebnisse * ", s)
        print("ACHTUNG: dann ALLE masse in", s, ", nicht mischen")

    modus = ask_int("modus 1=schwerpunkt 2=gegebene achsen: ")
    if modus != 1 and modus != 2:
        print("ungueltiger modus")
        return
    neg_global = ask_int("negative flaechen? 1=ja 0=nein: ")

    parts = []
    i = 0
    while True:
        i += 1
        print()
        print("teil", i)
        p = neuer_teil(neg_global)
        parts.append(p)
        print("  -> A=", round(p[0], 3),
              " y=", round(p[1], 3), " z=", round(p[2], 3))
        c = ask_int("noch ein teil? 1=ja 0=nein: ")
        if c == 0:
            break

    review(parts, modus, neg_global)


def steineranteil_interaktiv():
    print("steineranteil einzeln  (xx/q moeglich)")
    a = num("A_i=")
    y = num("y_i=")
    z = num("z_i=")
    ys = num("y_s_ges=")
    zs = num("z_s_ges=")
    iy_st, iz_st, iyz_st = steineranteil(a, y, z, ys, zs)
    print("Iy_st  =", iy_st)
    print("Iz_st  =", iz_st)
    print("Iyz_st =", iyz_st)


def main():
    try:
        st = load_state()
        if st is not None:
            parts, modus, neg = st
            print("letzte rechnung gefunden:", len(parts), "teile, modus", modus)
            if ask_int("laden? 1=ja 0=neue rechnung: ") == 1:
                print("(reload numerisch; parameter-label ggf. weg)")
                review(parts, modus, neg)
                return
        print("1 ftm (vollstaendig)   2 nur steineranteil   (q=beenden)")
        wahl = ask_int("wahl=")
        if wahl == 1:
            ftm()
        elif wahl == 2:
            steineranteil_interaktiv()
        else:
            print("ungueltige wahl")
    except Ende:
        print("== programm beendet (q) ==")


if __name__ == "__main__":
    main()
