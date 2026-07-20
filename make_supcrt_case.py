#!/usr/bin/env python3
"""
make_supcrt_case.py -- build a parallel case mirroring every core study deck but running on the
SUPCRT-HPT thermodynamic database instead of hanford.dat, to test whether diopside stays
supersaturated under an independent, high-P/T (HKF-coefficient) dataset.

Run from $MYSCRATCH/WAG (alongside 01_baseline/, ..., run_study_setonix.sh):  python3 make_supcrt_case.py

Each deck is copied to supcrt_case/<study>/decks/ with only these changes:
  * DATABASE hanford.dat  ->  GEOTHERMAL_HPT + DATABASE supcrt-hpt.dat
      GEOTHERMAL_HPT switches PFLOTRAN v6 to the HKF-coefficient reader (log K computed at cell P,T).
  * RENAME  : species/minerals named differently in SUPCRT-HPT.
  * REMOVE  : aqueous complexes SUPCRT-HPT does not define.
Nothing else (grid, injection, kinetics, surface areas) changes.

NOTE on the database file itself (done manually, once, NOT here): the distributed supcrt-hpt.dat
ships with a tabulated-style header ('temperature points' 8 ...) even though its body is
17-coefficient HPT. PFLOTRAN v6's reader takes that leading integer literally, so the file's first
line must be corrected to  'Number of Parameters' 17  before use. Its body is field-for-field
identical to the v6-validated geothermal-hpt.dat, so this is a header correction, not a data change.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE / "supcrt_case"
NEW_DB = "supcrt-hpt.dat"

RENAME = {
    "SiO2(am)": "Chalcedony",   # SUPCRT-HPT has no amorphous silica; chalcedony is the low-T proxy
    "H3SiO4-":  "HSiO3-",       # SUPCRT-HPT's silica-hydrolysis species
}
REMOVE_SECONDARY = {"CaHCO3+", "MgHCO3+", "Al(OH)2+", "Al(OH)3(aq)", "Al(OH)4-"}  # absent from SUPCRT-HPT

def core_studies():
    out = []
    for p in sorted(HERE.glob("*/decks")):
        name = p.parent.name
        if re.match(r"^\d{2}_", name) and not any(k in name for k in ("boundary", "database")):
            if list(p.glob("*.in")):
                out.append(name)
    return out

def transform(text):
    if "DATABASE hanford.dat" not in text:
        return None
    text = text.replace("DATABASE hanford.dat", f"GEOTHERMAL_HPT\n  DATABASE {NEW_DB}")
    for a, b in RENAME.items():
        text = text.replace(a, b)
    if REMOVE_SECONDARY:
        text = "\n".join(ln for ln in text.split("\n") if ln.strip() not in REMOVE_SECONDARY)
    return text

def main():
    studies = core_studies()
    if not studies:
        print("No core study decks found. Run from $MYSCRATCH/WAG."); return
    total = 0; per = {}
    for study in studies:
        src = HERE / study / "decks"
        dst = CASE / study / "decks"; dst.mkdir(parents=True, exist_ok=True)
        n = 0
        for deck in sorted(src.glob("*.in")):
            t = transform(deck.read_text())
            if t is None:
                print(f"    [skip] {study}/{deck.name}: no 'DATABASE hanford.dat'"); continue
            (dst / deck.name).write_text(t); n += 1
        per[study] = n; total += n
        print(f"  {study:<24}: {n:>3} decks")
    print(f"\n{total} decks -> {CASE}/  (DATABASE {NEW_DB}, GEOTHERMAL_HPT, renames+removals applied)")
    print("\nBefore running, ensure the header-corrected database is in place:")
    print("  cp /software/projects/pawsey1284/ychen6/pflotran-v6/pflotran/database/supcrt-hpt.dat supcrt-hpt.dat")
    print("  sed -i \"1s/.*/'Number of Parameters' 17/\" supcrt-hpt.dat")
    print("\nSmoke-test ONE deck (wait for the job to leave squeue before reading the log):")
    print("  sbatch --job-name=supcrt_test --array=0-0 --ntasks=16 --time=01:30:00 \\")
    print("    --export=ALL,DECKS_DIR=$PWD/supcrt_case/01_baseline/decks,\\")
    print("RUN_ROOT=$PWD/supcrt_case/01_baseline/runs,HANFORD_DB=$PWD/supcrt-hpt.dat run_study_setonix.sh")
    print("\nArray sizes for the full launch:")
    for s, n in per.items():
        print(f"    {s:<24} --array=0-{n-1}")

if __name__ == "__main__":
    main()
