# --- Import libraries ---
from pathlib import Path
import csv
from math import sin, cos, atan, radians, atan, tan, pi
import numpy as np

# --- Define file paths ---
root_dir = Path(__file__).resolve().parent.parent
harris_data = root_dir / "PHYS3116_Group_Code Files" / "HarrisParts.csv"
krause_data = root_dir / "PHYS3116_Group_Code Files" / "Krause21.csv"
vandenberg_data = root_dir / "PHYS3116_Group_Code Files" / "vandenberg_table2.csv"

#Information
def get_csv_data(file_path, skip_lines: int = 0):
    """
    Read a CSV file and return a list of rows (tuples).
    Row 0 is the header if your file has one.
    Use skip_lines>0 to skip header or extra lines at the top.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"CSV not found: {p}")

    csv_datapoints = []
    with open(p, "r", encoding="utf-8", newline="") as file:
        csv_reader = csv.reader(file)
        for _ in range(skip_lines):
            next(csv_reader, None)
        for row in csv_reader:
            csv_datapoints.append(tuple(row))
    return csv_datapoints

# References the data (Harris)
root_dir = Path(__file__).resolve().parent.parent
csv_data = root_dir / "PHYS3116_Group_Code Files" / "HarrisParts.csv"
csv_data_points = get_csv_data(csv_data, skip_lines=0)

# Krause21
krause_csv = root_dir / "PHYS3116_Group_Code Files" / "Krause21.csv"
krause_data_points = get_csv_data(krause_csv, skip_lines=0)

# VandenBerg (note the capital B in your file)
vandenberg_csv = root_dir / "PHYS3116_Group_Code Files" / "vandenBerg_table2.csv"
vandenberg_data_points = get_csv_data(vandenberg_csv, skip_lines=0)

# Column Data Values (Harris CSV)
"""
0  - ID
1  - Name
2  - (RA) Right Ascension
3  - (DEC) Declination
4  - (L) Galactic longitude
5  - (B) Galactic Latitude
6  - (R_Sun) Distance from Sun
7  - (R_gc) Distance from Galactic Centre
8  - (X) X-Galactic distance from Sun
9  - (Y) Y-Galactic distance from Sun
10 - (Z) Z-Galactic distance from Sun
11 - [Fe/H] Metallicity
12 - (wt) Weight of mean metallicity
13 - E(B-V) Foreground reddening
14 - (V_HB) V magnitude level of the Horizontal branch
15 - (m-M)V Apparent visual distance modulus
16 - (V_t) Integrated V magnitude of the cluster
17 - (M_V,t) Absolute visual magnitude (cluster luminosity), M_V,t = V_t - (m-M)V
18 - (U-B) Integrated color indices
19 - (B-V) Integrated color indices
20 - (V-R) Integrated color indices
21 - (V-I) Integrated color indices
22 - (spt) Spectral type of the integrated cluster light
23 - (ellip) Projected eccentricity of isophotes
24 - (v_r) Heliocentric radial velocity (km/s)
25 - (+/-) Observational (internal) uncertainty in radial velocity
26 - (v_LSR) Radial velocity relative to Solar neighborhood LSR
27 - (sig_v) Central velocity dispersion sig_v (km/s)
28 - (+/-) Observational (internal) uncertainty in velocity dispersion
29 - (c) King-model central concentration
30 - (r_c) Core radius in arcmin
31 - (r_h) Half-light radius in arcmin
32 - (mu_V) Central surface brightness, V magnitudes per square arcsecond
33 - (rho_0) Central luminosity density, log10(L☉ pc⁻3)
34 - (lg(tc)) Core relaxation time t(r_c), in log10(years)
35 - (lg(th)) Median relaxation time t(r_h), in log10(years)
"""

#   0  -> if the cell is blank
#   1  -> if the value is between -1.5 and -0.9 (inclusive)
#  -1  -> otherwise
def function(n, c):
    values = []
    for i in range(1, n + 1):                 # treat row 0 as header
        cell = csv_data_points[i][c]
        if len(cell) == 0:
            values.append(0)
        elif -1.5 <= float(cell) <= -0.9:
            values.append(1)
        else:
            values.append(-1)
    return values

def Name_list(n, c=0):
    values = []
    for i in range(1, n + 1):   # row 0 is header
        values.append(csv_data_points[i][c])
    return values

def Galatic_distance(n, c=7):
    values = []
    for i in range(1, n + 1):                 # row 0 is header
        cell = csv_data_points[i][c]
        if len(cell) == 0:
            values.append(0)
        elif float(cell) > 10:
            values.append(2)
        else:
            values.append(-2)
    return values

from math import atan, pi  # make sure this import is present at the top

def Gal_inclination(n, c1=7, c2=10):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        if len(csv_data_points[i][c1]) == 0 or len(csv_data_points[i][c2]) == 0:
            values.append(0)
        elif abs(atan(float(csv_data_points[i][c2]) / float(csv_data_points[i][c1])) * 180 / pi) > 60:
            values.append(1)
        elif abs(atan(float(csv_data_points[i][c2]) / float(csv_data_points[i][c1])) * 180 / pi) > 30:
            values.append(0)
        else:
            values.append(-1)
    return values
#Fe_H ratio
def Fe_H_ratio_function(n, c=11):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        if len(csv_data_points[i][c]) == 0:
            values.append(0)
        elif float(csv_data_points[i][c]) < -1.5 or -0.9 < float(csv_data_points[i][c]):
            values.append(2)
        else:
            values.append(-2)
    return values

def Projected_eccentricity(n, c=23):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        if len(csv_data_points[i][c]) == 0:
            values.append(0)
        elif float(csv_data_points[i][c]) < 0.1:
            values.append(0)
        elif float(csv_data_points[i][c]) < 0.2:
            values.append(0.5)
        elif float(csv_data_points[i][c]) < 0.7:
            values.append(1)
        else:
            values.append(1.5)
    return values
#Core radius

def Gal_Orbit2(n, l=4, s=6, g=7, x=8, y=9, v=24):
    """
    Direction of Galactic Orbit 
    Harris columns:
      l=L (deg), s=R_Sun, g=R_gc, x=X, y=Y, v=v_r (heliocentric km/s)
    Returns a score per row (excludes header).
    """
    values = []
    for i in range(1, n + 1):  # row 0 is header
        Ls = csv_data_points[i][l]
        Rs = csv_data_points[i][s]
        Rg = csv_data_points[i][g]
        Xs = csv_data_points[i][x]
        Ys = csv_data_points[i][y]
        Vr = csv_data_points[i][v]

        # blank checks 
        if (len(Ls) == 0 or len(Rs) == 0 or len(Rg) == 0 or
            len(Xs) == 0 or len(Ys) == 0 or len(Vr) == 0):
            values.append(0)
            continue

        # parse numerics
        L  = float(Ls)     # degrees
        RS = float(Rs)
        RG = float(Rg)
        X  = float(Xs)
        Y  = float(Ys)
        VR = float(Vr)

        # choose phi by the sign of Y (0 < phi < 360)
        phi = L if Y > 0 else (360.0 - L)

        # region branches
        if X >= 8:
            # L1 or R1
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))
            v_helio_corrected = VR * cos(psi)
            v_orbit = 44.0 * RG if RG < 5.0 else 220.0

            if (Y > 0 and VR > 0) or (Y < 0 and VR < 0):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) > v_orbit):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) < v_orbit):
                values.append(3)
            elif (Y < 0 and abs(v_helio_corrected) < v_orbit):
                values.append(-3)
            elif (Y < 0 and abs(v_helio_corrected) > v_orbit):
                values.append(3)
            else:
                values.append(0)

        elif 0 < X < 8:
            # L2 or R2
            psi = (8.0 / RS) * sin(radians(phi)) - (np.pi / 2.0)
            v_helio_corrected = VR * cos(psi)
            v_orbit = 44.0 * RG if RG < 5.0 else 220.0

            if (Y > 0 and VR > 0) or (Y < 0 and VR < 0):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) > v_orbit):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) < v_orbit):
                values.append(3)
            elif (Y < 0 and abs(v_helio_corrected) < v_orbit):
                values.append(-3)
            elif (Y < 0 and abs(v_helio_corrected) > v_orbit):
                values.append(3)
            else:
                values.append(0)

        else:
            # L3 or R3 (X <= 0)
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))
            v_helio_corrected = VR * cos(psi)
            v_orbit = 44.0 * RG if RG < 5.0 else 220.0

            if (Y > 0 and VR > 0) or (Y < 0 and VR < 0):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) < v_orbit):
                values.append(-3)
            elif (Y > 0 and abs(v_helio_corrected) > v_orbit):
                values.append(3)
            elif (Y < 0 and abs(v_helio_corrected) < v_orbit):
                values.append(-3)
            elif (Y < 0 and abs(v_helio_corrected) > v_orbit):
                values.append(3)
            else:
                values.append(0)

    return values





