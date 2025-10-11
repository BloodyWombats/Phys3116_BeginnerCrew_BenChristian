# --- Import libraries ---
from pathlib import Path
import csv
from math import sin, cos, atan, radians
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

