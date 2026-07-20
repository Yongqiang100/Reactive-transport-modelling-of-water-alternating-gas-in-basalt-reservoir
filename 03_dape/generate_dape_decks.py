#!/usr/bin/env python3
"""
generate_dape_decks.py — REBUILT from the working pflotran.in baseline

Creates 18 PFLOTRAN input decks for the Da-Pe disentangling experiment.

Suite A: Vary kinetic rate constants by factor κ (0.1, 0.3, 1.0, 3.0, 10.0)
         at fixed flow rate → varies Da while holding Pe constant.

Suite B: Vary molecular diffusion coefficient by factor D/D₀ (0.1, 1, 10, 100)
         at fixed flow rate and kinetics → holds Da nearly constant while
         varying Pe.

The two endmember configurations are tested:
  - dissolved: pure aqueous injection (MOLE_FRACTION 0.04, no free gas)
  - scco2:    pure supercritical CO2 injection (MOLE_FRACTION 0.99)

Total: 5 × 2 + 4 × 2 = 18 simulations.

Usage:
    python3 generate_dape_decks.py
Output:
    ./decks/*.in
"""

import os
import math
from pathlib import Path

# --------------------------------------------------------------------
# Baseline kinetic parameters (log₁₀ k₂₅ in mol/m²/s)
# Pulled directly from the working pflotran.in template
# --------------------------------------------------------------------
BASE_KINETICS = {
    'forsterite_acid':   {'log_k':  -6.05, 'Ea': 70.5e3, 'has_H': True,  'alpha':  0.50},
    'forsterite_alk':    {'log_k':  -4.65, 'Ea': 66.3e3, 'has_H': True,  'alpha':  0.25},
    'anorthite_neutral': {'log_k':  -9.12, 'Ea': 74.5e3, 'has_H': False, 'alpha':  0.0},
    'anorthite_acid':    {'log_k':  -3.50, 'Ea': 69.4e3, 'has_H': True,  'alpha':  1.411},
    'diopside_neutral':  {'log_k': -11.11, 'Ea': 40.6e3, 'has_H': False, 'alpha':  0.0},
    'diopside_acid':     {'log_k':  -6.36, 'Ea': 96.1e3, 'has_H': True,  'alpha':  0.7},
    'calcite':           {'log_k':  -8.0,  'Ea': 48.0e3},
    'magnesite':         {'log_k': -10.0,  'Ea': 62.0e3},
    'siderite':          {'log_k':  -9.0,  'Ea': 52.0e3},
    'dolomite_ord':      {'log_k': -11.0,  'Ea': 55.0e3},
    'kaolinite':         {'log_k': -13.0,  'Ea': 62.0e3},
    'sio2_am':           {'log_k': -10.0,  'Ea': 50.0e3},
}


def scale_log_k(log_k_base, kappa):
    """k_new = κ × k_base ⇒ log₁₀(k_new) = log₁₀(k_base) + log₁₀(κ)"""
    return log_k_base + math.log10(kappa)


# Simulated duration (total) and injection duration, in years. Injection runs
# for INJECTION_END_YR with each scenario's schedule, then stops; the run
# continues to FINAL_TIME_YR with zero injection (post-injection monitoring,
# Awolayo-style). Edit these two values to change the protocol; the WAG/adaptive
# schedules, the output times, and FINAL_TIME all follow automatically. Set
# INJECTION_END_YR = FINAL_TIME_YR for continuous injection throughout.
FINAL_TIME_YR = 100.0
INJECTION_END_YR = 30.0


def build_deck(scenario, kappa=1.0, diffusion=1.0e-9,
               nxyz="250 1 50",
               dxyz="    100@1.d0 50@3.d0 50@9.d0 50@26.d0\n    1.d0\n    50@2.d0",
               rate_mult=1.0, seed_carb_vf="0.0d0"):
    """
    scenario:  'dissolved' or 'scco2'
    kappa:     kinetic-rate multiplier (1.0 = baseline)
    diffusion: molecular diffusion coefficient (m²/s), baseline 1e-9
    nxyz/dxyz: grid specification (defaults reproduce the baseline mesh:
               2 m cells over 0-100 m near the well, widening to 30 m far-field;
               25 cells of 4 m in z). Override for grid-refinement studies.
    rate_mult: injection-rate multiplier (1.0 = baseline; 0.0 = no-injection
               control). Scales the well RATE only; kinetics/times unchanged.
    seed_carb_vf: initial volume fraction of the secondary carbonates
               (Calcite/Magnesite/Siderite/Dolomite-ord) in the formation_water
               constraint. CORRECTED DEFAULT '0.0d0' (zero-seed): carbonate VF then
               reflects only injection-driven precipitation. The original (now
               superseded) decks used '1.0d-4', which seeded ~20 m^3 of each phase
               out of equilibrium with the brine+silicates and re-equilibrated to
               dolomite (~6-27 t apparent carbonate, injection-independent),
               swamping the real uptake. Surface area is kept nonzero so kinetic
               nucleation still proceeds. Pass '1.0d-4' to reproduce the old runs.
    """
    def _fd(x):
        return f"{x:.6e}".replace("e", "d")
    def _ft(t):
        return f"{t:g}d0"
    # Simulated duration -> FINAL_TIME tag + output snapshot times (base set to
    # 30 yr, then every 25 yr out to FINAL_TIME_YR).
    _ft_tag = f"{FINAL_TIME_YR:g}.d0"
    _times = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
    _t = 50.0
    while _t <= FINAL_TIME_YR + 1e-9:
        _times.append(_t); _t += 25.0
    _out_times = " ".join(f"{v:g}" for v in _times if v <= FINAL_TIME_YR + 1e-9)
    # Apply κ to every mineral rate constant
    kin = {name: dict(props) for name, props in BASE_KINETICS.items()}
    for name in kin:
        kin[name]['log_k'] = scale_log_k(kin[name]['log_k'], kappa)

    # --- WAG schedule builders (manuscript base rates q_l=1e-5, q_g=5e-6) ----
    # Water slug:  q_l = 1e-5, x_CO2 = 0.04 ; Gas slug: q_g = 5e-6, x_CO2 = 0.99.
    # All scaled by rate_mult. SWAG co-injects both phases simultaneously.
    QL0, QG0 = 1.0e-5, 5.0e-6

    def _row_water(t):
        return (f"    {_ft(t)} {_fd(QL0 * rate_mult)} 0.d0 0.d0",
                f"    {_ft(t)} 0.04d0")

    def _row_gas(t):
        return (f"    {_ft(t)} 0.d0 {_fd(QG0 * rate_mult)} 0.d0",
                f"    {_ft(t)} 0.99d0")

    def _alt_schedule(half_cycle, total=INJECTION_END_YR):
        rr, rx, t, water = [], [], 0.0, True
        while t < total - 1e-9:
            r, x = _row_water(t) if water else _row_gas(t)
            rr.append(r); rx.append(x); water = not water; t += half_cycle
        return rr, rx

    def _adaptive_schedule():
        rr, rx = [], []
        def add(p): rr.append(p[0]); rx.append(p[1])
        add(_row_water(0.0))                       # 0-2 yr: dissolved priming
        t, water = 2.0, True                       # 2-10 yr: WAG 1:1 (6-mo)
        while t < 10.0 - 1e-9:
            add(_row_water(t) if water else _row_gas(t)); water = not water; t += 0.5
        t = 10.0                                   # 10-20 yr: WAG 1:2 (4mo/8mo)
        while t < 20.0 - 1e-9:
            add(_row_water(t)); add(_row_gas(t + 4.0 / 12.0)); t += 1.0
        t = 20.0                                   # 20 yr -> inj. end: WAG 1:3 (3mo/9mo)
        while t < INJECTION_END_YR - 1e-9:
            add(_row_water(t)); add(_row_gas(t + 3.0 / 12.0)); t += 1.0
        return rr, rx

    # Monitoring tail: injection drops to zero at INJECTION_END_YR and stays off
    # to FINAL_TIME_YR (rates held piecewise-constant between LIST entries).
    def _with_monitor(rr, rx):
        return (list(rr) + [f"    {_ft(INJECTION_END_YR)} 0.d0 0.d0 0.d0"],
                list(rx) + [f"    {_ft(INJECTION_END_YR)} 0.04d0"])

    def _list_block(rr, rx):
        rate = ("RATE LIST\n    TIME_UNITS yr\n    DATA_UNITS m^3/s m^3/s W\n"
                + "\n".join(rr) + "\n  /")
        xco2 = "MOLE_FRACTION LIST\n    TIME_UNITS yr\n" + "\n".join(rx) + "\n  /"
        return f"{rate}\n  {xco2}"

    # Fixed-rate scenarios as a pulse: constant (ql,qg) at xco2 until
    # INJECTION_END_YR, then zero (monitoring). Implemented as a step LIST.
    def _pulse(ql, qg, xco2):
        rr = [f"    0.d0 {_fd(ql)} {_fd(qg)} 0.d0",
              f"    {_ft(INJECTION_END_YR)} 0.d0 0.d0 0.d0"]
        rx = [f"    0.d0 {xco2}", f"    {_ft(INJECTION_END_YR)} {xco2}"]
        return _list_block(rr, rx)

    # Injection-scenario settings. Every scenario injects for INJECTION_END_YR
    # with its characteristic schedule, then monitors (zero rate) to FINAL_TIME_YR.
    if scenario == 'dissolved':
        injection_spec = _pulse(QL0 * rate_mult, 0.0, "0.04d0")
        scenario_label = "S1_DISSOLVED"
    elif scenario == 'scco2':
        injection_spec = _pulse(0.0, QG0 * rate_mult, "0.99d0")
        scenario_label = "S2_SCCO2"
    elif scenario == 'swag':
        # Simultaneous co-injection: q_l = 5e-6, q_g = 2.5e-6, x_CO2 = 0.35
        injection_spec = _pulse(5.0e-6 * rate_mult, 2.5e-6 * rate_mult, "0.35d0")
        scenario_label = "S5_SWAG"
    elif scenario == 'wag6mo':
        injection_spec = _list_block(*_with_monitor(*_alt_schedule(0.5)))
        scenario_label = "S3_WAG6MO"
    elif scenario == 'wag3mo':
        injection_spec = _list_block(*_with_monitor(*_alt_schedule(0.25)))
        scenario_label = "S4_WAG3MO"
    elif scenario == 'adaptive':
        injection_spec = _list_block(*_with_monitor(*_adaptive_schedule()))
        scenario_label = "S6_ADAPTIVE"
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    deck = f"""# =====================================================================
# Da-Pe Disentangling Sweep
# Scenario: {scenario_label}
# κ (kinetic-rate multiplier): {kappa}
# D (molecular diffusion, m^2/s): {diffusion:.3e}
# Generated by generate_dape_decks.py (rebuilt from working baseline)
# =====================================================================

SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE GENERAL
      OPTIONS
        GAS_COMPONENT_FORMULA_WEIGHT 44.d0
      /
    /
    SUBSURFACE_TRANSPORT transport
      MODE GIRT
    /
  /
END

SUBSURFACE

NUMERICAL_METHODS FLOW
  TIMESTEPPER
    TS_ACCELERATION 8
    MAX_TS_CUTS 40
    MAX_STEPS 1000000
  /
  NEWTON_SOLVER
    ATOL 1.d-10
    RTOL 1.d-8
    STOL 1.d-30
    MAXIMUM_NUMBER_OF_ITERATIONS 40
    MAXF 200
  /
  LINEAR_SOLVER
    SOLVER ITERATIVE
  /
END

NUMERICAL_METHODS TRANSPORT
  TIMESTEPPER
    TS_ACCELERATION 8
    MAX_TS_CUTS 40
    MAX_STEPS 1000000
  /
  NEWTON_SOLVER
    ATOL 1.d-10
    RTOL 1.d-8
    STOL 1.d-30
    MAXIMUM_NUMBER_OF_ITERATIONS 40
  /
  LINEAR_SOLVER
    SOLVER ITERATIVE
  /
END

GRID
  TYPE STRUCTURED
  NXYZ {nxyz}
  DXYZ
{dxyz}
  /
END

MATERIAL_PROPERTY basalt
  ID 1
  POROSITY 0.15d0
  TORTUOSITY 0.3d0
  ROCK_DENSITY 2900.d0
  SPECIFIC_HEAT 900.d0
  THERMAL_CONDUCTIVITY_DRY 1.5d0
  THERMAL_CONDUCTIVITY_WET 2.0d0
  PERMEABILITY
    PERM_ISO 1.d-13
  /
  CHARACTERISTIC_CURVES cc1
/

CHARACTERISTIC_CURVES cc1
  SATURATION_FUNCTION VAN_GENUCHTEN
    ALPHA 1.d-4
    M 0.457d0
    LIQUID_RESIDUAL_SATURATION 0.15d0
    MAX_CAPILLARY_PRESSURE 1.d7
  /
  PERMEABILITY_FUNCTION MUALEM_VG_LIQ
    M 0.457d0
    LIQUID_RESIDUAL_SATURATION 0.15d0
  /
  PERMEABILITY_FUNCTION MUALEM_VG_GAS
    M 0.457d0
    GAS_RESIDUAL_SATURATION 0.10d0
    LIQUID_RESIDUAL_SATURATION 0.15d0
  /
END

# -----------------------------------------------------------------
# FLUID_PROPERTY — the diffusion coefficient is the Suite-B knob
# -----------------------------------------------------------------
FLUID_PROPERTY
  DIFFUSION_COEFFICIENT {diffusion:.6e}
/

CHEMISTRY
  PRIMARY_SPECIES
    H+
    Ca++
    Mg++
    Fe++
    Na+
    K+
    Al+++
    SiO2(aq)
    HCO3-
    SO4--
    Cl-
    O2(aq)
  /
  SECONDARY_SPECIES
    OH-
    CO2(aq)
    CO3--
    CaHCO3+
    CaCO3(aq)
    CaSO4(aq)
    MgHCO3+
    MgCO3(aq)
    Fe+++
    AlOH++
    Al(OH)2+
    Al(OH)3(aq)
    Al(OH)4-
    HSO4-
    H3SiO4-
  /
  PASSIVE_GAS_SPECIES
    CO2(g)
    O2(g)
  /
  MINERALS
    Forsterite
    Anorthite
    Diopside
    Calcite
    Magnesite
    Siderite
    Dolomite-ord
    Kaolinite
    SiO2(am)
  /
  MINERAL_KINETICS
    Forsterite
      ! Rimstidt et al. (2012): acidic pH<5.6 [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['forsterite_acid']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['forsterite_acid']['Ea']:.1f}
        PREFACTOR_SPECIES H+
          ALPHA {kin['forsterite_acid']['alpha']:.4f}d0
        /
      /
      ! Rimstidt et al. (2012): alkaline pH>5.6 [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['forsterite_alk']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['forsterite_alk']['Ea']:.1f}
        PREFACTOR_SPECIES H+
          ALPHA {kin['forsterite_alk']['alpha']:.4f}d0
        /
      /
    /
    Anorthite
      ! Palandri & Kharaka (2004) Table 13: neutral mechanism [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['anorthite_neutral']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['anorthite_neutral']['Ea']:.1f}
      /
      ! Palandri & Kharaka (2004) Table 13: acid mechanism [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['anorthite_acid']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['anorthite_acid']['Ea']:.1f}
        PREFACTOR_SPECIES H+
          ALPHA {kin['anorthite_acid']['alpha']:.4f}d0
        /
      /
    /
    Diopside
      ! Palandri & Kharaka (2004): neutral mechanism [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['diopside_neutral']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['diopside_neutral']['Ea']:.1f}
      /
      ! Palandri & Kharaka (2004): acid mechanism [scaled by κ={kappa}]
      PREFACTOR
        RATE_CONSTANT {kin['diopside_acid']['log_k']:.4f}d0 mol/m^2-sec
        ACTIVATION_ENERGY {kin['diopside_acid']['Ea']:.1f}
        PREFACTOR_SPECIES H+
          ALPHA {kin['diopside_acid']['alpha']:.4f}d0
        /
      /
    /
    Calcite
      RATE_CONSTANT {kin['calcite']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['calcite']['Ea']:.1f}
    /
    Magnesite
      RATE_CONSTANT {kin['magnesite']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['magnesite']['Ea']:.1f}
    /
    Siderite
      RATE_CONSTANT {kin['siderite']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['siderite']['Ea']:.1f}
    /
    Dolomite-ord
      RATE_CONSTANT {kin['dolomite_ord']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['dolomite_ord']['Ea']:.1f}
    /
    Kaolinite
      RATE_CONSTANT {kin['kaolinite']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['kaolinite']['Ea']:.1f}
    /
    SiO2(am)
      RATE_CONSTANT {kin['sio2_am']['log_k']:.4f}d0 mol/m^2-sec
      ACTIVATION_ENERGY {kin['sio2_am']['Ea']:.1f}
    /
  /
  DATABASE hanford.dat
  LOG_FORMULATION
  UPDATE_POROSITY
  UPDATE_PERMEABILITY
  OUTPUT
    ALL
    TOTAL
    PH
    MINERAL_VOLUME_FRACTION
  /
END

REGION all
  COORDINATES
    0.d0 0.d0 0.d0
    2000.d0 1.d0 100.d0
  /
END

REGION injection_well
  COORDINATES
    0.d0 0.d0 20.d0
    20.d0 1.d0 80.d0
  /
END

REGION outlet
  FACE EAST
  COORDINATES
    2000.d0 0.d0 0.d0
    2000.d0 1.d0 100.d0
  /
END

REGION obs_10m
  COORDINATE 10.d0 0.5d0 50.d0
/
REGION obs_50m
  COORDINATE 50.d0 0.5d0 50.d0
/
REGION obs_100m
  COORDINATE 100.d0 0.5d0 50.d0
/
REGION obs_250m
  COORDINATE 250.d0 0.5d0 50.d0
/
REGION obs_500m
  COORDINATE 500.d0 0.5d0 50.d0
/
REGION obs_1000m
  COORDINATE 1000.d0 0.5d0 50.d0
/

FLOW_CONDITION initial
  TYPE
    LIQUID_PRESSURE HYDROSTATIC
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  DATUM 0.d0 0.d0 100.d0
  LIQUID_PRESSURE 6.0d6
  MOLE_FRACTION 1.d-8
  TEMPERATURE 60.d0
/

FLOW_CONDITION outlet_bc
  TYPE
    LIQUID_PRESSURE HYDROSTATIC
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  DATUM 0.d0 0.d0 100.d0
  LIQUID_PRESSURE 6.0d6
  MOLE_FRACTION 1.d-8
  TEMPERATURE 60.d0
/

TRANSPORT_CONDITION initial_tc
  TYPE DIRICHLET_ZERO_GRADIENT
  CONSTRAINT_LIST
    0.d0 formation_water
  /
END

TRANSPORT_CONDITION injection_tc
  TYPE DIRICHLET
  CONSTRAINT_LIST
    0.d0 co2_water
  /
END

TRANSPORT_CONDITION outlet_tc
  TYPE DIRICHLET_ZERO_GRADIENT
  CONSTRAINT_LIST
    0.d0 formation_water
  /
END

CONSTRAINT formation_water
  CONCENTRATIONS
    H+        7.5d0   pH
    Ca++      2.0d-3  T
    Mg++      1.5d-3  T
    Fe++      5.0d-5  T
    Na+       1.0d-2  T
    K+        5.0d-4  T
    Al+++     1.0d-8  T
    SiO2(aq)  5.0d-4  T
    HCO3-     3.0d-3  T
    SO4--     1.0d-3  T
    Cl-       5.0d-3  T
    O2(aq)    1.0d-6  T
  /
  MINERALS
    Forsterite   0.05d0  5.d4 m^2/m^3
    Anorthite    0.30d0  5.d4 m^2/m^3
    Diopside     0.25d0  5.d4 m^2/m^3
    Kaolinite    0.25d0  1.d-2 m^2/m^3
    Calcite      {seed_carb_vf}  5.d3 m^2/m^3
    Magnesite    {seed_carb_vf}  5.d3 m^2/m^3
    Siderite     {seed_carb_vf}  5.d3 m^2/m^3
    Dolomite-ord {seed_carb_vf}  5.d3 m^2/m^3
    SiO2(am)     0.0d0   1.d2 m^2/m^3
  /
END

CONSTRAINT co2_water
  CONCENTRATIONS
    H+        4.5d0   pH
    Ca++      1.0d-5  T
    Mg++      1.0d-5  T
    Fe++      1.0d-7  T
    Na+       1.0d-3  T
    K+        1.0d-4  T
    Al+++     1.0d-10 T
    SiO2(aq)  1.0d-6  T
    HCO3-     5.0d-2  T
    SO4--     1.0d-4  T
    Cl-       1.0d-3  T
    O2(aq)    5.0d-5  T
  /
END

# -----------------------------------------------------------------
# INJECTION FLOW CONDITION — scenario-specific
#   dissolved: liquid-only flow (1e-5 m^3/s), MOLE_FRACTION 0.04
#   scco2:    gas-only flow (5e-6 m^3/s),  MOLE_FRACTION 0.99
# -----------------------------------------------------------------
FLOW_CONDITION injection
  TYPE
    RATE SCALED_VOLUMETRIC_RATE VOLUME
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  {injection_spec}
  TEMPERATURE 60.d0
/

STRATA
  MATERIAL basalt
  REGION all
END

INITIAL_CONDITION
  FLOW_CONDITION initial
  TRANSPORT_CONDITION initial_tc
  REGION all
END

BOUNDARY_CONDITION outlet_side
  FLOW_CONDITION outlet_bc
  TRANSPORT_CONDITION outlet_tc
  REGION outlet
END

SOURCE_SINK well
  FLOW_CONDITION injection
  TRANSPORT_CONDITION injection_tc
  REGION injection_well
END

OBSERVATION
  REGION obs_10m
/
OBSERVATION
  REGION obs_50m
/
OBSERVATION
  REGION obs_100m
/
OBSERVATION
  REGION obs_250m
/
OBSERVATION
  REGION obs_500m
/
OBSERVATION
  REGION obs_1000m
/

OUTPUT
  TIMES yr {_out_times}
  FORMAT HDF5
  VELOCITY_AT_CENTER
  VARIABLES
    LIQUID_PRESSURE
    GAS_SATURATION
    TEMPERATURE
    POROSITY
    PERMEABILITY_X
  /
  OBSERVATION_FILE
    PERIODIC TIME 30 d
  /
  MASS_BALANCE_FILE
    PERIODIC TIME 90 d
  /
  SNAPSHOT_FILE
    FORMAT HDF5
    PERIODIC TIME 365 d
  /
END

TIME
  FINAL_TIME {_ft_tag} yr
  INITIAL_TIMESTEP_SIZE 1.d-3 d
  MAXIMUM_TIMESTEP_SIZE 5.d0 d
/

END_SUBSURFACE
"""
    return deck


def main():
    out = Path('decks')
    out.mkdir(exist_ok=True)

    # --- Suite A: kinetic rate sweep (Pe held constant) ---
    print("=== SUITE A: kinetic rate sweep (constant Pe) ===")
    for scenario in ['dissolved', 'scco2']:
        for kappa in [0.1, 0.3, 1.0, 3.0, 10.0]:
            deck = build_deck(scenario, kappa=kappa, diffusion=1.0e-9)
            kappa_tag = str(kappa).replace('.', 'p')
            fname = out / f"suiteA_{scenario}_kappa{kappa_tag}.in"
            fname.write_text(deck)
            print(f"  ✓ {fname.name}  (κ={kappa}, D=1e-9)")

    # --- Suite B: diffusion sweep (Da held constant) ---
    print("\n=== SUITE B: diffusion coefficient sweep (constant Da) ===")
    diffusions = {
        '0p1':   1.0e-10,
        '1p0':   1.0e-9,
        '10p0':  1.0e-8,
        '100p0': 1.0e-7,
    }
    for scenario in ['dissolved', 'scco2']:
        for tag, D in diffusions.items():
            deck = build_deck(scenario, kappa=1.0, diffusion=D)
            fname = out / f"suiteB_{scenario}_D{tag}.in"
            fname.write_text(deck)
            print(f"  ✓ {fname.name}  (κ=1, D={D:.0e})")

    print(f"\nTotal decks generated: 18")
    print(f"Output directory: {out.absolute()}")


if __name__ == "__main__":
    main()
