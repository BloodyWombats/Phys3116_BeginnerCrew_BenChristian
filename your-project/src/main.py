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

# --- Function to safely read CSV files ---
def get_csv_data(file_path):
    """Reads a CSV and returns a list of rows."""
    if not file_path.exists():
        print(f"⚠️  File not found: {file_path}")
        return []
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        data = [row for row in reader]
    print(f"✅ Loaded {file_path.name} with {len(data)-1} data rows.")
    return data

# --- Load all your tables ---
harris_table = get_csv_data(harris_data)
krause_table = get_csv_data(krause_data)
vandenberg_table = get_csv_data(vandenberg_data)

# --- Confirm everything loaded ---
if not harris_table or not krause_table or not vandenberg_table:
    print("\n❌ One or more files did not load correctly. Check paths or file names.")
else:
    print("\n🎉 All three CSV files loaded successfully!")
    print("Harris rows:", len(harris_table)-1)
    print("Krause rows:", len(krause_table)-1)
    print("Vandenberg rows:", len(vandenberg_table)-1)

