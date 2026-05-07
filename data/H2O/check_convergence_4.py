"""
check_convergence_4.py
-----------------------
Checks all .out files in 4/CCSD_T/aug-cc-pvdz/ for the presence of
SCF, MP2, and CCSD(T) energy lines.

Run from bsse_db/data/H2O/
"""

from pathlib import Path

here     = Path(__file__).resolve().parent        # H2O/
out_dir  = here / "4" / "CCSD_T" / "aug-cc-pvdz"

scf_markers   = ["Total SCF energy", "Total RHF energy"]
mp2_markers   = ["Total MP2 energy"]
ccsdt_markers = ["CCSD(T) total energy", "Total CCSD(T) energy"]

missing_scf   = []
missing_mp2   = []
missing_ccsdt = []

out_files = sorted(out_dir.glob("*.out"))
print(f"Checking {len(out_files)} output files...\n")

for f in out_files:
    text = f.read_text(errors="replace")

    has_scf   = any(m in text for m in scf_markers)
    has_mp2   = any(m in text for m in mp2_markers)
    has_ccsdt = any(m in text for m in ccsdt_markers)

    if not has_scf:
        missing_scf.append(f.name)
    if not has_mp2:
        missing_mp2.append(f.name)
    if not has_ccsdt:
        missing_ccsdt.append(f.name)

# --- Report ---
print(f"Missing SCF   ({len(missing_scf)}):")
for name in missing_scf:
    print(f"  {name}")

print(f"\nMissing MP2   ({len(missing_mp2)}):")
for name in missing_mp2:
    print(f"  {name}")

print(f"\nMissing CCSD(T) ({len(missing_ccsdt)}):")
for name in missing_ccsdt:
    print(f"  {name}")

if not any([missing_scf, missing_mp2, missing_ccsdt]):
    print("All 65 files have all three energy markers. ✅")