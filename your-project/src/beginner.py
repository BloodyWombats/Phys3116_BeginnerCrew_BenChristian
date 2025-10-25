from pathlib import Path
import csv, re
from math import sin, cos, radians
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "PHYS3116_Group_Code Files"
OUT  = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

HARRIS_CSV = DATA / "HarrisParts.csv"
VBERG_CSV  = DATA / "vandenberg_table2.csv"
KRAUSE_CSV = DATA / "Krause21.csv"

def read_csv(path):
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for row in r:
            rows.append(tuple(row))
    return rows

def to_float(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None

def canon_ngc(name):
    if not name:
        return None
    m = re.search(r'\bNGC\s*([0-9]{1,4})\b', str(name).upper())
    if m:
        return f"NGC {int(m.group(1))}"
    return None

def canon_simple(name):
    return re.sub(r"\s+", " ", str(name).strip().upper()) if name else None

def canon_key(name):
    return canon_ngc(name) or canon_simple(name)

harris = read_csv(HARRIS_CSV)
vberg  = read_csv(VBERG_CSV)
krause = read_csv(KRAUSE_CSV)
if not harris: raise SystemExit("HarrisParts.csv is empty/missing")
if not vberg:  raise SystemExit("vandenberg_table2.csv is empty/missing")
if not krause: raise SystemExit("Krause21.csv is empty/missing")

H_L     = 4
H_RSUN  = 6
H_RGC   = 7
H_X     = 8
H_Y     = 9
H_VR    = 24
H_IDCOL = 0

def harris_rotation_flags(hrows):
    flags = {}
    if len(hrows) <= 1:
        return flags
    for i in range(1, len(hrows)):
        row = hrows[i]
        name = (row[H_IDCOL] if len(row) > H_IDCOL else "").strip()
        key  = canon_key(name)
        if not key:
            continue
        def gf(idx):
            return to_float(row[idx]) if idx < len(row) else None
        L  = gf(H_L);  RS = gf(H_RSUN); RG = gf(H_RGC)
        X  = gf(H_X);  Y  = gf(H_Y);    VR = gf(H_VR)
        if None in (L, RS, RG, X, Y, VR) or RS == 0:
            flags[key] = 0
            continue
        phi = L if Y > 0 else (360.0 - L)
        if X >= 8:
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))
        elif 0 < X < 8:
            psi = (8.0 / RS) * sin(radians(phi)) - (np.pi / 2.0)
        else:
            psi = np.pi/2.0 - (8.0 / RS) * sin(radians(phi))
        v_corr = VR * cos(psi)
        v_orb  = 44.0 * RG if RG < 5.0 else 220.0
        if (Y > 0 and abs(v_corr) < v_orb) or (Y < 0 and abs(v_corr) > v_orb):
            flags[key] = +1
        elif (Y > 0 and abs(v_corr) > v_orb) or (Y < 0 and abs(v_corr) < v_orb):
            flags[key] = -1
        else:
            flags[key] = 0
    return flags

rot_flag = harris_rotation_flags(harris)

vb_header = [h.strip() for h in vberg[0]]
vb_idx = {h.lower(): i for i, h in enumerate(vb_header)}
VB_NAME = vb_idx.get("name", 1)
VB_FEH  = vb_idx.get("feh", 2)
VB_AGE  = vb_idx.get("age", 3)

def vberg_amr_scores(vrows):
    points = []
    for i in range(1, len(vrows)):
        r = vrows[i]
        name = (r[VB_NAME] if VB_NAME is not None and VB_NAME < len(r) else "").strip()
        key  = canon_key(name)
        if not key:
            continue
        feh = to_float(r[VB_FEH]) if VB_FEH is not None else None
        age = to_float(r[VB_AGE]) if VB_AGE is not None else None
        if feh is not None and age is not None and -2.5 <= feh <= -0.5:
            points.append((key, feh, age))
    if len(points) >= 5:
        X = np.array([p[1] for p in points])
        Y = np.array([p[2] for p in points])
        A = np.vstack([np.ones(len(X)), X]).T
        coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
        a, b = coef[0], coef[1]
    else:
        a, b = 12.0, 0.0
    scores = {}
    for key, feh, age in points:
        pred = a + b * feh
        resid = pred - age
        if resid >= 0.8:   scr = 2
        elif resid >= 0.4: scr = 1
        else:              scr = 0
        scores[key] = (scr, resid)
    return scores

amr = vberg_amr_scores(vberg)

kr_header = [h.strip() for h in krause[0]]
kr_idx = {h.lower(): i for i, h in enumerate(kr_header)}
KR_NAME = kr_idx.get("name", 1)
KR_FEH  = kr_idx.get("feh", 2)
KR_AGE  = kr_idx.get("age", 3)

def krause_scores(krows):
    out = {}
    for i in range(1, len(krows)):
        r = krows[i]
        name = (r[KR_NAME] if KR_NAME is not None and KR_NAME < len(r) else "").strip()
        key  = canon_key(name)
        if not key:
            continue
        age = to_float(r[KR_AGE]) if KR_AGE is not None else None
        feh = to_float(r[KR_FEH]) if KR_FEH is not None else None
        if age is None or feh is None:
            out[key] = out.get(key, 0) + 0
            continue
        if (age > 12.5 and feh > -0.9) or (age <= 12.5 and feh < -1.0):
            out[key] = out.get(key, 0) + 1
        else:
            out[key] = out.get(key, 0) - 1
    return out

krause_score = krause_scores(krause)

combined = {}
details  = {}

for key, rf in rot_flag.items():
    add = 1 if rf == -1 else 0
    combined[key] = combined.get(key, 0) + add
    details[key]  = [rf, 0, None, 0]

for key, (scr, resid) in amr.items():
    combined[key] = combined.get(key, 0) + scr
    if key not in details:
        details[key] = [0, scr, resid, 0]
    else:
        details[key][1] = scr
        details[key][2] = resid

for key, ks in krause_score.items():
    combined[key] = combined.get(key, 0) + ks
    if key not in details:
        details[key] = [0, 0, None, ks]
    else:
        details[key][3] = ks

def label(score):
    if score >= 4: return "Certain"
    if score >= 2: return "Likely"
    if score >= 1: return "Possible"
    if score >= 0: return "Unlikely"
    return "Not_GC"

final = sorted(combined.items(), key=lambda t: t[1], reverse=True)

out_path = OUT / "Combined_Simple_WithKrause.csv"
with open(out_path, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow([
        "ID/Name", "TotalScore", "Bin",
        "RotationFlag(-1/0/+1)",
        "VandenBerg_AMR_Score(0-2)", "AMR_AgeResidual_Gyr(+=>younger)",
        "Krause_Score(-1/0/+1)"
    ])
    for key, total in final:
        rf, a_s, a_r, k_s = details.get(key, [0,0,None,0])
        w.writerow([key, total, label(total), rf, a_s,
                    (None if a_r is None else round(a_r, 2)), k_s])

print("Saved:", out_path)
print("Top 10:", final[:10])
