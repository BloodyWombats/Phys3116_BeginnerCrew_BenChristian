# --- Import libraries ---
from pathlib import Path
import csv
from math import sin, cos, atan, radians
import numpy as np

# ---------------------------
# Helpers / CSV loading
# ---------------------------
def get_csv_data(file_path, skip_lines: int = 0):
    """
    Read a CSV file and return a list of tuples (rows).
    Uses UTF-8 with BOM handling and skips `skip_lines` rows if asked.
    """
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
    """Safe float converter. Returns None for blanks/junk."""
    try:
        return float(str(s).strip())
    except Exception:
        return None
    # ---------------------------
# Paths & data loads (your names)
# ---------------------------
root_dir = Path(__file__).resolve().parent.parent

harris_data     = root_dir / "PHYS3116_Group_Code Files" / "HarrisParts.csv"
krause_data     = root_dir / "PHYS3116_Group_Code Files" / "Krause21.csv"
vandenberg_data = root_dir / "PHYS3116_Group_Code Files" / "vandenberg_table2.csv"

# Harris is the initial working table
csv_data_points = get_csv_data(harris_data, skip_lines=0)
krause_data_points     = get_csv_data(krause_data, skip_lines=0)
vandenberg_data_points = get_csv_data(vandenberg_data, skip_lines=0)

# ---------------------------
# Harris column map (comment)
# 0:ID 1:Name 2:RA 3:DEC 4:L 5:B 6:R_Sun 7:R_gc 8:X 9:Y 10:Z
# 11:[Fe/H] 12:wt 13:E(B-V) 14:V_HB 15:(m-M)V 16:V_t 17:M_V,t
# 18:U-B 19:B-V 20:V-R 21:V-I 22:spt 23:ellip 24:v_r 25:+/- (v_r)
# 26:v_LSR 27:sig_v 28:+/- (sig_v) 29:c 30:r_c 31:r_h 32:mu_V
# 33:rho_0 34:lg(tc) 35:lg(th)
# ---------------------------

# ---------------------------
# Harris feature functions
# ---------------------------
def function(n, c):
    """
    PDF’s simple scorer:
      0  -> blank
      1  -> -1.5 <= value <= -0.9
     -1  -> otherwise
    (made robust with to_float)
    """
    values = []
    for i in range(1, n + 1):
        v = to_float(csv_data_points[i][c]) if c < len(csv_data_points[i]) else None
        if v is None:
            values.append(0)
        elif -1.5 <= v <= -0.9:
            values.append(1)
        else:
            values.append(-1)
    return values

def Name_list(n, c=0):
    values = []
    for i in range(1, n + 1):
        values.append(csv_data_points[i][c] if c < len(csv_data_points[i]) else "")
    return values

def Galatic_distance(n, c=7):
    values = []
    for i in range(1, n + 1):
        val = to_float(csv_data_points[i][c]) if c < len(csv_data_points[i]) else None
        if val is None:
            values.append(0)       # blank/junk -> neutral
        elif val > 10:
            values.append(2)
        else:
            values.append(-2)
    return values

def Gal_inclination(n, c1=7, c2=10):
    values = []
    for i in range(1, n + 1):
        r_gc = to_float(csv_data_points[i][c1]) if c1 < len(csv_data_points[i]) else None
        zval = to_float(csv_data_points[i][c2]) if c2 < len(csv_data_points[i]) else None
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
    for i in range(1, n + 1):
        feh = to_float(csv_data_points[i][c]) if c < len(csv_data_points[i]) else None
        if feh is None:
            values.append(0)          # blank/junk -> neutral
        elif feh < -1.5 or feh > -0.9:
            values.append(2)
        else:
            values.append(-2)
    return values

def Projected_eccentricity(n, c=23):
    values = []
    for i in range(1, n + 1):
        val = to_float(csv_data_points[i][c]) if c < len(csv_data_points[i]) else None
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

def Gal_Orbit2(n, l=4, s=6, g=7, x=8, y=9, v=24):
    """
    Direction of Galactic Orbit (robust to short/messy rows).
    Uses the PDF’s L1/L2/L3 vs R1/R2/R3 logic with phi & psi as specified.
    """
    values = []
    for i in range(1, n + 1):
        row = csv_data_points[i]

        # ensure needed columns exist
        needed = (l, s, g, x, y, v)
        if any(idx >= len(row) for idx in needed):
            values.append(0)
            continue

        Ls, Rs, Rg, Xs, Ys, Vr = row[l], row[s], row[g], row[x], row[y], row[v]

        # blank checks
        if (Ls == "" or Rs == "" or Rg == "" or Xs == "" or Ys == "" or Vr == ""):
            values.append(0)
            continue

        # parse
        try:
            L  = float(Ls)   # deg
            RS = float(Rs)
            RG = float(Rg)
            X  = float(Xs)
            Y  = float(Ys)
            VR = float(Vr)
        except Exception:
            values.append(0)
            continue

        phi = L if Y > 0 else (360.0 - L)

        if X >= 8:
            # L1 / R1
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
            # L2 / R2
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
            # L3 / R3  (X <= 0)
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

# ---------------------------
# Main: scoring, combine, export
# ---------------------------
if __name__ == "__main__":
    # quick smoke
    n_harris = len(csv_data_points) - 1
    print("Harris rows (excl header):", n_harris)

    # ----- Harris scoring -----
    harris_features = [
        Galatic_distance(n_harris, c=7),                 # R_gc
        Gal_inclination(n_harris, c1=7, c2=10),          # inclination
        Fe_H_ratio_function(n_harris, c=11),             # [Fe/H]
        Projected_eccentricity(n_harris, c=23),          # ellip
        Gal_Orbit2(n_harris, l=4, s=6, g=7, x=8, y=9, v=24),  # orbit
    ]
    Sum_harris = [sum(x) for x in zip(*harris_features)]

    HARRIS_NAME_SCORE = list(
        zip([csv_data_points[i][0] for i in range(1, n_harris + 1)], Sum_harris)
    )
    HARRIS_NAME_SCORE.sort(key=lambda t: t[1], reverse=True)

    print("\n=== HARRIS (top 5) ===")
    print(HARRIS_NAME_SCORE[:5])
    print("=== HARRIS (bottom 5) ===")
    print(HARRIS_NAME_SCORE[-5:])

    # ----- Krause scoring -----
    csv_data_points = krause_data_points   # switch working table
    n_krause = len(krause_data_points) - 1

    def Age_and_Fe_H_ratio_function(n, c1=6, c2=7):
        vals = []
        for i in range(1, n + 1):
            a = to_float(csv_data_points[i][c1]) if c1 < len(csv_data_points[i]) else None
            f = to_float(csv_data_points[i][c2]) if c2 < len(csv_data_points[i]) else None
            if a is None or f is None:
                vals.append(0)
            elif (a > 12.5 and f > -0.9) or (a <= 12.5 and f < -1.0):
                vals.append(1)
            else:
                vals.append(-1)
        return vals

    KRAUSE_SCORES = Age_and_Fe_H_ratio_function(n_krause, c1=6, c2=7)
    KRAUSE_NAME_SCORE = list(
        zip([krause_data_points[i][1] for i in range(1, n_krause + 1)], KRAUSE_SCORES)
    )
    KRAUSE_NAME_SCORE.sort(key=lambda t: t[1], reverse=True)

    print("\n=== KRAUSE (top 5) ===")
    print(KRAUSE_NAME_SCORE[:5])
    print("=== KRAUSE (bottom 5) ===")
    print(KRAUSE_NAME_SCORE[-5:])

    # ----- Combine by ID/name -----
    combined = {}
    for name, s in HARRIS_NAME_SCORE:
        combined[name] = combined.get(name, 0.0) + float(s)
    for name, s in KRAUSE_NAME_SCORE:
        combined[name] = combined.get(name, 0.0) + float(s)

    FINAL = sorted(combined.items(), key=lambda t: t[1], reverse=True)

    # ----- Bin labels (like the PDF) -----
    def label_bin(score: float) -> str:
        if score >= 7.5: return "Certain"
        if 3 <= score < 7.5: return "Likely"
        if -1 <= score < 3: return "Possible"
        if -5.5 <= score < -1: return "Unlikely"
        return "Not_GC"

    BINNED = [(name, score, label_bin(score)) for name, score in FINAL]

    print("\n=== BINS (first few) ===")
    print(BINNED[:10])

    # ----- OPTIONAL: Join VandenBerg onto FINAL by cluster name -----
    def find_col_idx(header_row, candidates):
        """Return first matching column index by (case-insensitive) header name."""
        hdr_lc = [h.strip().lower() for h in header_row]
        for want in candidates:
            w = want.strip().lower()
            for i, h in enumerate(hdr_lc):
                if h == w:
                    return i
        return None

    vb_rows = vandenberg_data_points
    vb_header = vb_rows[0] if len(vb_rows) > 0 else []

    name_idx = find_col_idx(vb_header, ["Name", "Cluster", "ID"]) or 0
    age_idx  = find_col_idx(vb_header, ["Age", "Age(Gyr)", "Age_gyr", "Age [Gyr]"])
    feh_idx  = find_col_idx(vb_header, ["[Fe/H]", "FeH", "Fe_H", "Fe/H"])

    vb_dict = {}
    for r in vb_rows[1:]:
        nm = r[name_idx] if name_idx < len(r) else ""
        if not nm:
            continue
        age_val = r[age_idx] if (age_idx is not None and age_idx < len(r)) else ""
        feh_val = r[feh_idx] if (feh_idx is not None and feh_idx < len(r)) else ""
        vb_dict[nm.strip().lower()] = (age_val, feh_val)

    COMBINED_WITH_VB = []
    for name, score in FINAL:
        age_val, feh_val = ("", "")
        tup = vb_dict.get(name.strip().lower())
        if tup is not None:
            age_val, feh_val = tup
        COMBINED_WITH_VB.append((name, score, age_val, feh_val))

    vb_age_header = "VB_Age" if age_idx is not None else "VB_Age(NA)"
    vb_feh_header = "VB_[Fe/H]" if feh_idx is not None else "VB_[Fe/H](NA)"

    # ----- Save CSVs for Excel -----
    out_dir = root_dir / "outputs"
    out_dir.mkdir(exist_ok=True)

    with open(out_dir / "Harris_Scores.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Harris_Score"])
        w.writerows(HARRIS_NAME_SCORE)

    with open(out_dir / "Krause_Scores.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Krause_Score"])
        w.writerows(KRAUSE_NAME_SCORE)

    with open(out_dir / "Combined_Scores.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Combined_Score"])
        w.writerows(FINAL)

    with open(out_dir / "Binned_Combined.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Combined_Score", "Bin"])
        w.writerows(BINNED)

    with open(out_dir / "Combined_with_VandenBerg.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "Combined_Score", vb_age_header, vb_feh_header])
        w.writerows(COMBINED_WITH_VB)

    print(f"\n📁 Saved outputs to: {out_dir}")