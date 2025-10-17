# --- Import libraries ---
from pathlib import Path
import csv
from math import sin, cos, atan, radians, tan, pi
import numpy as np

# --- helpers/loaders ---
def get_csv_data(file_path, skip_lines: int = 0):
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"CSV not found: {p}")
    rows = []
    with open(p, "r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for _ in range(skip_lines):
            next(r, None)
        for row in r:
            rows.append(tuple(row))
    return rows

def to_float(s):
    try:
        return float(str(s).strip())
    except Exception:
        return None
    
# --- Define file paths ---
root_dir = Path(__file__).resolve().parent.parent
harris_data = root_dir / "PHYS3116_Group_Code Files" / "HarrisParts.csv"
krause_data = root_dir / "PHYS3116_Group_Code Files" / "Krause21.csv"
vandenberg_data = root_dir / "PHYS3116_Group_Code Files" / "vandenberg_table2.csv"

# Harris -> working table
csv_data_points = get_csv_data(harris_data, skip_lines=0)

# Krause21 + VandenBerg
krause_data_points     = get_csv_data(krause_data, skip_lines=0)
vandenberg_data_points = get_csv_data(vandenberg_data, skip_lines=0)

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
    for i in range(1, n + 1):  # row 0 is header
        val = to_float(csv_data_points[i][c])
        if val is None:
            values.append(0)       # blank or junk like "c:" -> neutral
        elif val > 10:
            values.append(2)
        else:
            values.append(-2)
    return values


def Gal_inclination(n, c1=7, c2=10):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        r_gc = to_float(csv_data_points[i][c1])
        zval = to_float(csv_data_points[i][c2])
        if r_gc is None or zval is None or r_gc == 0:
            values.append(0)
            continue
        angle = abs(atan(zval / r_gc) * 180.0 / np.pi)
        if angle > 60:
            values.append(1)
        elif angle > 30:
            values.append(0)
        else:
            values.append(-1)
    return values

def Fe_H_ratio_function(n, c=11):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        feh = to_float(csv_data_points[i][c])
        if feh is None:
            values.append(0)          # blank/junk like 'rho_0' -> neutral
        elif feh < -1.5 or feh > -0.9:
            values.append(2)
        else:
            values.append(-2)
    return values

def Projected_eccentricity(n, c=23):
    values = []
    for i in range(1, n + 1):  # row 0 is header
        val = to_float(csv_data_points[i][c])
        if val is None:
            values.append(0)
        elif val < 0.1:
            values.append(0)
        elif val < 0.2:
            values.append(0.5)
        elif val < 0.7:
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

# --- quick smoke test (optional) ---
n_harris = len(csv_data_points) - 1
print("Harris rows (excl header):", n_harris)

# ===== HARRIS scoring (csv_data_points currently = Harris) =====
harris_features = [
    Galatic_distance(n_harris, c=7),                 # R_gc
    Gal_inclination(n_harris, c1=7, c2=10),         # inclination
    Fe_H_ratio_function(n_harris, c=11),            # [Fe/H]
    Projected_eccentricity(n_harris, c=23),         # ellip
    Gal_Orbit2(n_harris, l=4, s=6, g=7, x=8, y=9, v=24),  # orbit
]
Sum_harris = [sum(x) for x in zip(*harris_features)]
HARRIS_NAME_SCORE = list(zip([csv_data_points[i][0] for i in range(1, n_harris + 1)], Sum_harris))
HARRIS_NAME_SCORE.sort(key=lambda t: t[1], reverse=True)

print("\n=== HARRIS (top 5) ===")
print(HARRIS_NAME_SCORE[:5])
print("=== HARRIS (bottom 5) ===")
print(HARRIS_NAME_SCORE[-5:])

# ===== KRAUSE scoring =====
csv_data_points = krause_data_points
n_krause = len(krause_data_points) - 1

def Age_and_Fe_H_ratio_function(n, c1=6, c2=7):
    vals = []
    for i in range(1, n + 1):
        a = to_float(csv_data_points[i][c1])
        f = to_float(csv_data_points[i][c2])
        if a is None or f is None:
            vals.append(0)
        elif (a > 12.5 and f > -0.9) or (a <= 12.5 and f < -1.0):
            vals.append(1)
        else:
            vals.append(-1)
    return vals

KRAUSE_SCORES = Age_and_Fe_H_ratio_function(n_krause, c1=6, c2=7)
KRAUSE_NAME_SCORE = list(zip([krause_data_points[i][1] for i in range(1, n_krause + 1)], KRAUSE_SCORES))
KRAUSE_NAME_SCORE.sort(key=lambda t: t[1], reverse=True)

print("\n=== KRAUSE (top 5) ===")
print(KRAUSE_NAME_SCORE[:5])
print("=== KRAUSE (bottom 5) ===")
print(KRAUSE_NAME_SCORE[-5:])

# ===== Combine by ID + bucket like PDF =====
combined = {}
for name, s in HARRIS_NAME_SCORE:
    combined[name] = combined.get(name, 0.0) + float(s)
for name, s in KRAUSE_NAME_SCORE:
    combined[name] = combined.get(name, 0.0) + float(s)

FINAL = sorted(combined.items(), key=lambda t: t[1], reverse=True)

certain   = [n for n, s in FINAL if s >= 7.5]
likely    = [n for n, s in FINAL if 3 <= s < 7.5]
possible  = [n for n, s in FINAL if -1 <= s < 3]
unlikely  = [n for n, s in FINAL if -5.5 <= s < -1]
not_gc    = [n for n, s in FINAL if s < -5.5]

print("\n=== BINS ===")
print("Certain :", certain[:10], ("... (+more)" if len(certain) > 10 else ""))
print("Likely  :", likely[:10], ("... (+more)" if len(likely) > 10 else ""))
print("Possible:", possible[:10], ("... (+more)" if len(possible) > 10 else ""))
print("Unlikely:", unlikely[:10], ("... (+more)" if len(unlikely) > 10 else ""))
print("Not_GC  :", not_gc[:10], ("... (+more)" if len(not_gc) > 10 else ""))

# ===== Save CSVs for Excel =====
out_dir = Path(__file__).resolve().parent.parent / "outputs"
out_dir.mkdir(exist_ok=True)

with open(out_dir / "Harris_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["ID", "Harris_Score"]); w.writerows(HARRIS_NAME_SCORE)

with open(out_dir / "Krause_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["ID", "Krause_Score"]); w.writerows(KRAUSE_NAME_SCORE)

with open(out_dir / "Combined_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["ID", "Combined_Score"]); w.writerows(FINAL)

print(f"\n📁 Saved outputs to: {out_dir}")
