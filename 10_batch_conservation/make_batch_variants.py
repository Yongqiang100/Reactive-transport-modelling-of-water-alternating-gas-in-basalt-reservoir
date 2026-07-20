#!/usr/bin/env python3
"""Build carbon-coupling fix variants from decks/batch_cell.in."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); DK = os.path.join(HERE,"decks")
base = open(os.path.join(DK,"batch_cell.in")).read()
def retitle(s, t): return s.replace("SINGLE-CELL BATCH CARBON-CONSERVATION TEST", t)

v1 = base.replace("  PASSIVE_GAS_SPECIES\n    CO2(g)\n    O2(g)\n  /",
                  "  PASSIVE_GAS_SPECIES\n    O2(g)\n  /")
assert v1 != base, "v1 patch failed"
open(os.path.join(DK,"v1_noco2gas.in"),"w").write(retitle(v1,"VARIANT v1: GENERAL, no reactive CO2(g) passive gas"))

v2 = base
v2 = v2.replace("  PRIMARY_SPECIES\n    H+\n    Ca++\n    Mg++\n    Fe++\n    Na+\n    K+\n    Al+++\n    SiO2(aq)\n    HCO3-\n    SO4--\n    Cl-\n    O2(aq)\n  /",
                "  PRIMARY_SPECIES\n    H+\n    Ca++\n    Mg++\n    Fe++\n    Na+\n    K+\n    Al+++\n    SiO2(aq)\n    CO2(aq)\n    SO4--\n    Cl-\n    O2(aq)\n  /")
v2 = v2.replace("  SECONDARY_SPECIES\n    OH-\n    CO2(aq)\n    CO3--",
                "  SECONDARY_SPECIES\n    OH-\n    HCO3-\n    CO3--")
v2 = v2.replace("    HCO3-     3.0d-3  T", "    CO2(aq)   3.0d-3  T")
assert v2.count("CO2(aq)") >= 2 and "\n    HCO3-\n    CO3--" in v2, "v2 swap failed"
open(os.path.join(DK,"v2_co2aq_primary.in"),"w").write(retitle(v2,"VARIANT v2: GENERAL, CO2(aq) primary carbon species"))

v3 = base.replace("      MODE GENERAL", "      MODE TH")
v3 = v3.replace("  PASSIVE_GAS_SPECIES\n    CO2(g)\n    O2(g)\n  /\n", "")
v3 = v3.replace(
"""FLOW_CONDITION initial
  TYPE
    LIQUID_PRESSURE DIRICHLET
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 6.0d6
  MOLE_FRACTION 1.d-8
  TEMPERATURE 60.d0
/""",
"""FLOW_CONDITION initial
  TYPE
    LIQUID_PRESSURE DIRICHLET
    TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 6.0d6
  TEMPERATURE 60.d0
/""")
assert "MODE TH" in v3 and "MOLE_FRACTION" not in v3 and "CO2(g)" not in v3, "v3 patch failed"
open(os.path.join(DK,"v3_th_nogas.in"),"w").write(retitle(v3,"VARIANT v3: single-phase TH, no gas component (diagnostic)"))
print("wrote v1_noco2gas.in, v2_co2aq_primary.in, v3_th_nogas.in")
