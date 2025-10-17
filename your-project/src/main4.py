# --- Import libraries ---
from pathlib import Path
import csv
from math import sin, cos, atan, radians, pi
import numpy as np

# ========== helpers ==========
def get_csv_data(file_path, skip_lines: int = 0):
    """Read a CSV file (UTF-8/UTF-8-SIG) into a list of tuples."""
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
    """Best-effort float parse; returns None on failure/blank."""
    try:
        return float(str(s).strip())
    except Exception:
        return None

def get_cell_safe(row, idx):
    """Return cell if index exists, else ''."""
    return row[idx] if 0 <= idx < len(row) else ""

def get_float_safe(row, idx):
    """Return float value for cell idx, or None if missing/bad."""
    return to_float(get_cell_safe(row, idx))

# ========== paths & loads (your exact layout) ==========
root_dir   = Path(__file__).resolve().parent.parent
data_dir   = root_dir / "PHYS3116_Group_Code Files"
harris_p   = data_dir / "HarrisParts.csv"
krause_p   = data_dir / "Krause21.csv"
vberg_p    = data_dir / "vandenberg_table2.csv"

harris_data_points = get_csv_data(harris_p, skip_lines=0)
krause_data_points = get_csv_data(krause_p, skip_lines=0)
vberg_data_points  = get_csv_data(vberg_p,  skip_lines=0)

# ========== HARRIS column map (by index in your HarrisParts.csv) ==========
# (matches your earlier doc)
# 4 L, 6 R_Sun, 7 R_gc, 8 X, 9 Y, 10 Z, 11 [Fe/H], 23 ellip, 24 v_r
H_L     = 4
H_RSUN  = 6
H_RGC   = 7
H_X     = 8
H_Y     = 9
H_Z     = 10
H_FEH   = 11
H_ELLIP = 23
H_VR    = 24
H_NAME  = 0  # we’ll use column 0 as the ID/Name key

# ========== HARRIS scoring functions ==========
def Name_list(n, c=H_NAME, table=None):
    t = table if table is not None else harris_data_points
    return [get_cell_safe(t[i], c) for i in range(1, n + 1)]

def Galactic_distance(n, c=H_RGC, table=None):
    t = table if table is not None else harris_data_points
    vals = []
    for i in range(1, n + 1):
        v = get_float_safe(t[i], c)
        if v is None:
            vals.append(0)
        elif v > 10:
            vals.append(2)
        else:
            vals.append(-2)
    return vals

def Gal_inclination(n, c1=H_RGC, c2=H_Z, table=None):
    t = table if table is not None else harris_data_points
    vals = []
    for i in range(1, n + 1):
        r_gc = get_float_safe(t[i], c1)
        zval = get_float_safe(t[i], c2)
        if r_gc is None or zval is None or r_gc == 0:
            vals.append(0)
            continue
        angle = abs(atan(zval / r_gc) * 180.0 / np.pi)
        if angle > 60:
            vals.append(1)
        elif angle > 30:
            vals.append(0)
        else:
            vals.append(-1)
    return vals

def Fe_H_ratio_function(n, c=H_FEH, table=None):
    t = table if table is not None else harris_data_points
    vals = []
    for i in range(1, n + 1):
        feh = get_float_safe(t[i], c)
        if feh is None:
            vals.append(0)
        elif feh < -1.5 or feh > -0.9:
            vals.append(2)
        else:
            vals.append(-2)
    return vals

def Projected_eccentricity(n, c=H_ELLIP, table=None):
    t = table if table is not None else harris_data_points
    vals = []
    for i in range(1, n + 1):
        e = get_float_safe(t[i], c)
        if e is None:
            vals.append(0)
        elif e < 0.1:
            vals.append(0)
        elif e < 0.2:
            vals.append(0.5)
        elif e < 0.7:
            vals.append(1)
        else:
            vals.append(1.5)
    return vals

def Gal_Orbit2(n, l=H_L, s=H_RSUN, g=H_RGC, x=H_X, y=H_Y, v=H_VR, table=None):
    """
    Direction-of-rotation proxy from the PDF logic (cleaned for Python).
    Returns a signed score per row: positive ~ co-rotating, negative ~ counter.
    """
    t = table if table is not None else harris_data_points
    out = []
    for i in range(1, n + 1):
        L  = get_float_safe(t[i], l)
        RS = get_float_safe(t[i], s)
        RG = get_float_safe(t[i], g)
        X  = get_float_safe(t[i], x)
        Y  = get_float_safe(t[i], y)
        VR = get_float_safe(t[i], v)

        if None in (L, RS, RG, X, Y, VR) or RS == 0:
            out.append(0)
            continue

        phi = L if Y > 0 else (360.0 - L)

        if X >= 8:
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))
        elif 0 < X < 8:
            psi = (8.0 / RS) * sin(radians(phi)) - (np.pi / 2.0)
        else:
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))

        v_helio_corrected = VR * cos(psi)
        v_orbit = 44.0 * RG if RG < 5.0 else 220.0

        if (Y > 0 and VR > 0) or (Y < 0 and VR < 0):
            out.append(-3)
        elif (Y > 0 and abs(v_helio_corrected) > v_orbit):
            out.append(-3)
        elif (Y > 0 and abs(v_helio_corrected) < v_orbit):
            out.append(3)
        elif (Y < 0 and abs(v_helio_corrected) < v_orbit):
            out.append(-3)
        elif (Y < 0 and abs(v_helio_corrected) > v_orbit):
            out.append(3)
        else:
            out.append(0)
    return out

# ========== HARRIS scoring ==========
n_harris = len(harris_data_points) - 1
HARRIS_IDS = Name_list(n_harris, c=H_NAME, table=harris_data_points)

harris_feats = [
    Galactic_distance(n_harris, table=harris_data_points),         # R_gc
    Gal_inclination(n_harris, table=harris_data_points),           # inclination
    Fe_H_ratio_function(n_harris, table=harris_data_points),       # [Fe/H]
    Projected_eccentricity(n_harris, table=harris_data_points),    # ellip
    Gal_Orbit2(n_harris, table=harris_data_points),                # rotation proxy
]
Sum_harris = [sum(x) for x in zip(*harris_feats)]
HARRIS_NAME_SCORE = list(zip(HARRIS_IDS, Sum_harris))
HARRIS_NAME_SCORE.sort(key=lambda t: t[1], reverse=True)

# rotation flag by ID (for the report)
gal_orb = harris_feats[-1]
rot_flag_by_id = {}
for idx, sc in enumerate(gal_orb):
    cid = HARRIS_IDS[idx]
    if sc >= 2:
        rot_flag_by_id[cid] = 1   # co-rotating-ish
    elif sc <= -2:
        rot_flag_by_id[cid] = -1  # counter-ish
    else:
        rot_flag_by_id[cid] = 0   # unknown/weak

# ========== VANDENBERG 2013: AMR accretion signal ==========
vb_header = [h.strip() for h in vberg_data_points[0]]
vb_idx = {h.lower(): i for i, h in enumerate(vb_header)}
VB_NAME = vb_idx.get("name", 1)
VB_FEH  = vb_idx.get("feh", 2)
VB_AGE  = vb_idx.get("age", 3)
VB_HBT  = vb_idx.get("hbtype", 8)
VB_RG   = vb_idx.get("r_g", 9)

def vberg_rows_numeric(tbl):
    rows = []
    for i in range(1, len(tbl)):
        name = str(get_cell_safe(tbl[i], VB_NAME)).strip()
        feh  = get_float_safe(tbl[i], VB_FEH)
        age  = get_float_safe(tbl[i], VB_AGE)
        hb   = get_float_safe(tbl[i], VB_HBT) if VB_HBT is not None else None
        rg   = get_float_safe(tbl[i], VB_RG)  if VB_RG  is not None else None
        if name:
            rows.append((name, feh, age, hb, rg))
    return rows

vb_rows = vberg_rows_numeric(vberg_data_points)

def amr_outlier_scores(rows):
    """
    Fit Age ~ a + b*[Fe/H] then residual = pred - obs (positive => younger).
    Score:
      +2 if residual >= 0.8 Gyr
      +1 if 0.4 <= residual < 0.8
      +1 if (Fe/H <= -1.3 and HBtype <= 0)    # red HB at low metallicity
      +1 if R_G >= 15 kpc                      # outer halo
    """
    xs, ys = [], []
    for name, feh, age, hb, rg in rows:
        if feh is None or age is None:
            continue
        if -2.4 <= feh <= -0.8:
            xs.append(feh); ys.append(age)
    if len(xs) >= 5:
        A = np.vstack([np.ones(len(xs)), xs]).T
        coef, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
        a, b = coef[0], coef[1]
    else:
        a, b = 12.0, 0.0  # fallback flat line

    out = []
    for name, feh, age, hb, rg in rows:
        score = 0
        residual = None
        if feh is not None and age is not None:
            pred = a + b * feh
            residual = pred - age
            if residual >= 0.8:
                score += 2
            elif residual >= 0.4:
                score += 1
        if feh is not None and hb is not None and feh <= -1.3 and hb <= 0:
            score += 1
        if rg is not None and rg >= 15:
            score += 1
        out.append((name, float(score), residual))
    out.sort(key=lambda t: t[1], reverse=True)
    return out

VBERG_AMR = amr_outlier_scores(vb_rows)

# ========== (Optional) KRAUSE tiny score ==========
# Expect header like: #NGC,Name,FeH,Age,Age_err,Method,Figs,Range,HBtype,R_G,M_V,v_e0,log_sigma_0
kr_header = [h.strip() for h in krause_data_points[0]]
kr_idx = {h.lower(): i for i, h in enumerate(kr_header)}
KR_NAME = kr_idx.get("name", 1)
KR_AGE  = kr_idx.get("age", 3)
KR_FEH  = kr_idx.get("feh", 2)

def Age_and_Fe_H_ratio_function_krause(tbl):
    vals = []
    for i in range(1, len(tbl)):
        name = str(get_cell_safe(tbl[i], KR_NAME)).strip()
        a    = get_float_safe(tbl[i], KR_AGE)
        f    = get_float_safe(tbl[i], KR_FEH)
        if not name:
            continue
        if a is None or f is None:
            score = 0
        elif (a > 12.5 and f > -0.9) or (a <= 12.5 and f < -1.0):
            score = 1
        else:
            score = -1
        vals.append((name, score))
    vals.sort(key=lambda t: t[1], reverse=True)
    return vals

KRAUSE_NAME_SCORE = Age_and_Fe_H_ratio_function_krause(krause_data_points)

# ========== Combine Harris + VandenBerg (+ Krause) ==========
combined_scores = {name: float(s) for name, s in HARRIS_NAME_SCORE}

# add Krause
for name, s in KRAUSE_NAME_SCORE:
    combined_scores[name] = combined_scores.get(name, 0.0) + float(s)

# add VandenBerg AMR (weight 1.0; change if you want)
VB_WEIGHT = 1.0
amr_score = {}
amr_resid = {}
for name, s, r in VBERG_AMR:
    combined_scores[name] = combined_scores.get(name, 0.0) + VB_WEIGHT * float(s)
    amr_score[name] = float(s)
    amr_resid[name] = None if r is None else float(r)

FINAL = sorted(combined_scores.items(), key=lambda t: t[1], reverse=True)

def bin_label(score):
    if score >= 9.0:  return "Certain"
    if score >= 4.0:  return "Likely"
    if score >= 0.0:  return "Possible"
    if score >= -4.0: return "Unlikely"
    return "Not_GC"

final_rows = []
for name, total in FINAL:
    rot_flag = rot_flag_by_id.get(name, 0)
    vb_s     = amr_score.get(name, 0.0)
    vb_res   = amr_resid.get(name, None)
    final_rows.append([
        name,
        round(total, 2),
        bin_label(total),
        rot_flag,
        vb_s,
        (None if vb_res is None else round(vb_res, 2)),
    ])

# ========== Save outputs ==========
out_dir = root_dir / "outputs"
out_dir.mkdir(exist_ok=True)

# Harris only
with open(out_dir / "Harris_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["ID/Name", "Harris_Score"])
    w.writerows(HARRIS_NAME_SCORE)

# Krause only
with open(out_dir / "Krause_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["Name", "Krause_Score"])
    w.writerows(KRAUSE_NAME_SCORE)

# VandenBerg AMR table
with open(out_dir / "VandenBerg_AMR_Scores.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["Name", "AMR_Score(0-4)", "AgeResidual_Gyr(+=>younger)"])
    for name, s, r in VBERG_AMR:
        w.writerow([name, s, (None if r is None else round(r, 2))])

# Combined final
with open(out_dir / "Combined_With_VandenBerg.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["ID/Name", "TotalScore", "Bin", "RotationFlag", "AMR_Score", "AgeResidual_Gyr(+=>younger)"])
    w.writerows(final_rows)

print("\n=== HARRIS (top 5) ===")
print(HARRIS_NAME_SCORE[:5])
print("=== VandenBerg AMR (top 5) ===")
print(VBERG_AMR[:5])
print("=== Combined (top 10) ===")
print(FINAL[:10])
print(f"\n📁 Saved outputs to: {out_dir}")
