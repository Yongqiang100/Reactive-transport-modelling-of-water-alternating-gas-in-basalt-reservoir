#!/usr/bin/env python3
"""Regenerate the single-cell batch conservation deck from the CORRECTED build_deck.
Run from anywhere; it imports build_deck from ../03_dape/generate_dape_decks.py,
extracts the CHEMISTRY / database / material / brine-constraint blocks VERBATIM, and
wraps them in a 1-cell closed reactor. Writes decks/batch_cell.in."""
import importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "..", "03_dape", "generate_dape_decks.py")
spec = importlib.util.spec_from_file_location("gdd", GEN)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
deck = m.build_deck("dissolved", rate_mult=0.0)     # corrected zero-seed, no injection
L = deck.splitlines()
HEADERS = ["NUMERICAL_METHODS FLOW","NUMERICAL_METHODS TRANSPORT","GRID","MATERIAL_PROPERTY",
           "CHARACTERISTIC_CURVES","FLUID_PROPERTY","CHEMISTRY","OUTPUT","REGION ","FLOW_CONDITION ",
           "TRANSPORT_CONDITION ","CONSTRAINT ","STRATA","INITIAL_CONDITION","BOUNDARY_CONDITION",
           "SOURCE_SINK","TIME","END_SUBSURFACE"]
starts=[]
for i,ln in enumerate(L):
    if ln and not ln[0].isspace():
        for h in HEADERS:
            if ln.startswith(h): starts.append((i, ln.strip())); break
def card(hs):
    for k,(i,h) in enumerate(starts):
        if h.startswith(hs) or h==hs:
            nxt = starts[k+1][0] if k+1 < len(starts) else len(L)
            return "\n".join(L[i:nxt]).rstrip()
nm_flow=card("NUMERICAL_METHODS FLOW"); nm_tran=card("NUMERICAL_METHODS TRANSPORT")
mat=card("MATERIAL_PROPERTY"); cc=card("CHARACTERISTIC_CURVES"); fluid=card("FLUID_PROPERTY")
chem=card("CHEMISTRY"); tc_init=card("TRANSPORT_CONDITION initial_tc"); con_fw=card("CONSTRAINT formation_water")
batch = f"""# SINGLE-CELL BATCH CARBON-CONSERVATION TEST (identical chemistry, 1 m^3 closed cell)
SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE GENERAL
    /
    SUBSURFACE_TRANSPORT transport
      MODE GIRT
    /
  /
END

SUBSURFACE

{nm_flow}

{nm_tran}

GRID
  TYPE STRUCTURED
  NXYZ 1 1 1
  DXYZ
    1.d0
    1.d0
    1.d0
  /
END

{mat}

{cc}

{fluid}

{chem}

REGION all
  COORDINATES
    -1.d20 -1.d20 -1.d20
     1.d20  1.d20  1.d20
  /
END

FLOW_CONDITION initial
  TYPE
    LIQUID_PRESSURE DIRICHLET
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 6.0d6
  MOLE_FRACTION 1.d-8
  TEMPERATURE 60.d0
/

{tc_init}

{con_fw}

STRATA
  REGION all
  MATERIAL basalt
END

INITIAL_CONDITION
  FLOW_CONDITION initial
  TRANSPORT_CONDITION initial_tc
  REGION all
END

OUTPUT
  SNAPSHOT_FILE
    PERIODIC TIME 50. y
    FORMAT HDF5
  /
  MASS_BALANCE_FILE
    PERIODIC TIME 5. y
  /
END

TIME
  FINAL_TIME 100. y
  INITIAL_TIMESTEP_SIZE 1.d-3 y
  MAXIMUM_TIMESTEP_SIZE 1. y
END

END_SUBSURFACE
"""
os.makedirs(os.path.join(HERE,"decks"), exist_ok=True)
open(os.path.join(HERE,"decks","batch_cell.in"),"w").write(batch)
print("wrote", os.path.join(HERE,"decks","batch_cell.in"), "(%d lines)"%len(batch.splitlines()))
