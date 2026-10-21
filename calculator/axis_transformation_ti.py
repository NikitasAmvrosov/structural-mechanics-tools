"""
Axis transformation for second moments of area and bending moments.

  1 = rotate Iy, Iz, Iyz by an angle phi
  2 = principal moments I1, I2 and principal angle phi*
  3 = rotate the bending moments My, Mz by phi

Convention: z points down, Iyz = -integral(y*z dA).
TI-Nspire CX II version (runs on PC too).
Author: Nikitas Amvrosov
"""

from math import pi, sin, cos, atan2, sqrt

def deg2rad(phi):
    return phi * pi / 180.0

def ftm_transform():
    print("")
    print("--- FTM Transformation ---")

    Iyy = float(input("Iyy = "))
    Izz = float(input("Izz = "))
    Iyz = float(input("Iyz = "))
    phi = float(input("phi_deg = "))

    p = deg2rad(phi)

    I_etaeta = 0.5*(Iyy + Izz) + 0.5*(Iyy - Izz)*cos(2*p) + Iyz*sin(2*p)
    I_zetazeta = 0.5*(Iyy + Izz) - 0.5*(Iyy - Izz)*cos(2*p) - Iyz*sin(2*p)
    I_etazeta = -0.5*(Iyy - Izz)*sin(2*p) + Iyz*cos(2*p)

    print("")
    print("Ergebnis:")
    print("I_etaeta =", I_etaeta)
    print("I_zetazeta =", I_zetazeta)
    print("I_etazeta =", I_etazeta)

def haupttraegheit():
    print("")
    print("--- Haupttraegheitsmomente ---")

    Iyy = float(input("Iyy = "))
    Izz = float(input("Izz = "))
    Iyz = float(input("Iyz = "))

    root = sqrt((Iyy - Izz)**2 + 4*Iyz**2)

    I1 = 0.5*(Iyy + Izz) + 0.5*root
    I2 = 0.5*(Iyy + Izz) - 0.5*root
    phi = 0.5 * atan2(2*Iyz, Iyy - Izz) * 180.0 / pi

    print("")
    print("Ergebnis:")
    print("I1 =", I1)
    print("I2 =", I2)
    print("phi_star_deg =", phi)

def moment_transform():
    print("")
    print("--- Momententransformation ---")

    My = float(input("My = "))
    Mz = float(input("Mz = "))
    phi = float(input("phi_deg = "))

    p = deg2rad(phi)

    M_eta = My*cos(p) + Mz*sin(p)
    M_zeta = -My*sin(p) + Mz*cos(p)

    print("")
    print("Ergebnis:")
    print("M_eta =", M_eta)
    print("M_zeta =", M_zeta)

def trafo():
    while True:
        print("=================================")
        print("TRAFO TOOL")
        print("=================================")
        print("1 = FTM Transformation")
        print("2 = Haupttraegheitsmomente")
        print("3 = Momententransformation")
        print("0 = Exit")

        choice = int(input("Was willst du machen? "))

        if choice == 1:
            ftm_transform()
        elif choice == 2:
            haupttraegheit()
        elif choice == 3:
            moment_transform()
        elif choice == 0:
            print("Bye")
            break
        else:
            print("Ungueltig")


trafo()
