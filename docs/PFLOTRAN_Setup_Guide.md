# Running the WAG CO₂ Mineralization Simulations with PFLOTRAN

## Complete Guide: Installation, Execution, and Post-Processing

**Reference:** Chen et al. (2026) *Earth-Science Reviews* — Carbon Mineralization in Mafic and Ultramafic Rocks

---

## 1. System Requirements

PFLOTRAN is a massively parallel reactive transport simulator written in Fortran 90. It requires:

- **Operating system:** Linux (Ubuntu 20.04+, CentOS 7+, Rocky 8+). macOS works but is less tested. Windows requires WSL2.
- **Compiler:** GNU Fortran (gfortran ≥ 9.0) or Intel Fortran (ifort ≥ 2021).
- **MPI:** OpenMPI ≥ 4.0 or MPICH ≥ 3.3.
- **PETSc:** Version 3.18–3.21 (PFLOTRAN's primary dependency — it provides the parallel linear algebra solvers).
- **HDF5:** Version 1.12+ with parallel I/O support (for output files).
- **RAM:** Minimum 8 GB; 32+ GB recommended for the full 100×50 grid.
- **Disk:** ~5 GB per scenario for HDF5 snapshots over 30 years.
- **Cores:** 4–16 cores for a single scenario; 32–128 for production runs.

If you have access to a university HPC cluster (e.g., Pawsey in Perth, NCI in Canberra), PFLOTRAN or PETSc may already be available as a module.

---

## 2. Installation

### Option A: Using an HPC Module (Fastest)

```bash
# Check if PFLOTRAN is available on your cluster
module spider pflotran
module spider petsc

# If available:
module load pflotran/4.0    # or latest version
which pflotran              # verify installation
```

### Option B: Install via PETSc + PFLOTRAN from Source

This is the standard method. The entire process takes 1–3 hours depending on your system.

#### Step 1: Install System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran gcc g++ \
  openmpi-bin libopenmpi-dev \
  cmake git wget curl python3 python3-pip \
  liblapack-dev libblas-dev

# CentOS/Rocky
sudo yum groupinstall "Development Tools"
sudo yum install gcc-gfortran openmpi openmpi-devel cmake git wget
export PATH=/usr/lib64/openmpi/bin:$PATH
```

#### Step 2: Set Environment Variables

```bash
# Create a working directory
export PFLOTRAN_DIR=$HOME/pflotran-build
mkdir -p $PFLOTRAN_DIR
cd $PFLOTRAN_DIR

# Set PETSc environment
export PETSC_DIR=$PFLOTRAN_DIR/petsc
export PETSC_ARCH=arch-linux-gnu-opt
```

#### Step 3: Download and Build PETSc

```bash
cd $PFLOTRAN_DIR
git clone -b release https://gitlab.com/petsc/petsc.git
cd petsc

# Configure PETSc with HDF5 (PFLOTRAN requirement)
# This will automatically download and build HDF5
./configure \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --with-debugging=0 \
  --with-shared-libraries=0 \
  --download-hdf5=1 \
  --download-hdf5-fortran-bindings=1 \
  --download-metis=1 \
  --download-parmetis=1 \
  COPTFLAGS="-O3" \
  CXXOPTFLAGS="-O3" \
  FOPTFLAGS="-O3"

# Build PETSc (takes 20-60 minutes)
make all

# Test the installation
make check
```

#### Step 4: Download and Build PFLOTRAN

```bash
cd $PFLOTRAN_DIR
git clone https://bitbucket.org/pflotran/pflotran.git
cd pflotran/src/pflotran

# Build PFLOTRAN
make pflotran

# Verify build
ls -la pflotran    # Should see the executable

# Add to PATH
export PATH=$PFLOTRAN_DIR/pflotran/src/pflotran:$PATH
echo 'export PATH=$HOME/pflotran-build/pflotran/src/pflotran:$PATH' >> ~/.bashrc
```

#### Step 5: Verify Installation

```bash
# Test with a simple example
cd $PFLOTRAN_DIR/pflotran/regression_tests/default/543
pflotran -input_prefix 543_flow
# Should run without errors and produce output files
```

### Option C: Docker Container (Easiest for Local Testing)

```bash
# Pull the PFLOTRAN Docker image
docker pull pflotran/pflotran:latest

# Run interactively
docker run -it -v $(pwd)/wag_simulations:/work pflotran/pflotran:latest /bin/bash

# Inside the container, PFLOTRAN is pre-installed
which pflotran
```

---

## 3. Preparing the Simulation Files

### Directory Structure

```
wag_simulations/
├── database/
│   └── hanford.dat              # Thermodynamic database
├── scenario1_dissolved/
│   ├── basalt_wag_base.in       # Base input deck
│   └── scenario1_dissolved.in   # Scenario-specific overrides
├── scenario2_scco2/
│   ├── basalt_wag_base.in
│   └── scenario2_scco2.in
├── scenario3_wag6mo/
│   ├── basalt_wag_base.in
│   └── scenario3_wag6mo.in
├── scenario4_wag3mo/
│   ├── basalt_wag_base.in
│   └── scenario4_wag3mo.in
├── scenario5_swag/
│   ├── basalt_wag_base.in
│   └── scenario5_swag.in
├── scenario6_adaptive/
│   ├── basalt_wag_base.in
│   └── scenario6_adaptive.in
└── postprocess/
    └── wag_optimizer.py
```

### Getting the Thermodynamic Database

The input deck references `hanford.dat`. PFLOTRAN ships with several databases:

```bash
# The database is included in the PFLOTRAN source
cp $PFLOTRAN_DIR/pflotran/database/hanford.dat ./database/

# Alternative databases (may need for specific minerals):
# - pflotran/database/thermodynamic_database.dat
# - pflotran/database/CrunchFlow/datacom.dat
```

**Important database note:** The provided input decks use mineral names that must match the database entries exactly. If you use a different database (e.g., `thermoddem.dat`), you may need to adjust mineral names. Common mismatches include `Basalt_glass` (not in all databases — may need to define as a custom kinetic mineral) and `Smectite-Na` (naming varies across databases). Check the database file for available mineral names:

```bash
grep -i "forsterite\|anorthite\|calcite\|magnesite\|siderite" database/hanford.dat
```

### Modifying the Input Deck for Your Database

If mineral names don't match, edit the `CHEMISTRY` block in `basalt_wag_base.in`:

```
# In the MINERALS section, replace names to match your database.
# For example, if your database uses "Olivine_Fo" instead of "Forsterite":
MINERALS
  Olivine_Fo          # was: Forsterite
  Anorthite
  ...
```

---

## 4. Running the Simulations

### Single Scenario (Serial or Small Parallel)

```bash
cd wag_simulations/scenario1_dissolved/

# Serial run (1 core) — good for testing
pflotran -input_prefix basalt_wag_base

# Parallel run (4 cores)
mpirun -np 4 pflotran -input_prefix basalt_wag_base

# Parallel run (16 cores, recommended for production)
mpirun -np 16 pflotran -input_prefix basalt_wag_base
```

**Runtime estimates** for the 100×1×50 grid over 30 years:

| Cores | Estimated Wall Time |
|-------|-------------------|
| 1     | 24–72 hours       |
| 4     | 8–24 hours        |
| 16    | 2–8 hours         |
| 64    | 0.5–2 hours       |

### Running All 6 Scenarios (Batch Script)

Create a bash script `run_all.sh`:

```bash
#!/bin/bash
#============================================================
# Run all 6 WAG scenarios sequentially or in parallel
#============================================================

PFLOTRAN_EXE=$(which pflotran)
NCORES=16
BASE_DIR=$(pwd)

scenarios=(
  "scenario1_dissolved"
  "scenario2_scco2"
  "scenario3_wag6mo"
  "scenario4_wag3mo"
  "scenario5_swag"
  "scenario6_adaptive"
)

for scenario in "${scenarios[@]}"; do
  echo "============================================"
  echo "Running: $scenario"
  echo "Start: $(date)"
  echo "============================================"
  
  cd $BASE_DIR/$scenario
  
  # Run PFLOTRAN
  mpirun -np $NCORES $PFLOTRAN_EXE -input_prefix basalt_wag_base \
    > run_${scenario}.log 2>&1
  
  # Check exit status
  if [ $? -eq 0 ]; then
    echo "  SUCCESS: $scenario completed at $(date)"
  else
    echo "  FAILED: $scenario — check run_${scenario}.log"
  fi
  
  cd $BASE_DIR
done

echo "All scenarios complete."
```

```bash
chmod +x run_all.sh
./run_all.sh
```

### HPC Job Submission (SLURM)

For Pawsey (Setonix), NCI (Gadi), or similar SLURM clusters:

```bash
#!/bin/bash
#SBATCH --job-name=wag_s1
#SBATCH --account=your_project_code
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

# Load modules (adjust for your cluster)
module load pflotran/4.0
# OR if building from source:
# module load gcc/12.2.0 openmpi/4.1.4 petsc/3.20.0
# export PATH=$HOME/pflotran-build/pflotran/src/pflotran:$PATH

cd $SLURM_SUBMIT_DIR

srun pflotran -input_prefix basalt_wag_base
```

Submit with:

```bash
# Submit all scenarios
for s in scenario{1..6}_*; do
  cd $s
  sbatch ../submit.slurm
  cd ..
done

# Monitor jobs
squeue -u $USER
```

### Pawsey-Specific Notes (Perth)

```bash
# On Setonix at Pawsey Supercomputing Centre
module load petsc/3.20.5
module load hdf5/1.14.3

# If PFLOTRAN is not a module, build from source using Pawsey's PETSc
cd $MYSOFTWARE
git clone https://bitbucket.org/pflotran/pflotran.git
cd pflotran/src/pflotran
make pflotran

# Submit to work queue
sbatch --account=your_pawsey_project submit.slurm
```

---

## 5. Monitoring Running Simulations

### Check Progress

```bash
# PFLOTRAN writes progress to stdout/log
tail -f run_scenario1.log

# Look for timestep information:
# "== GENERAL FLOW ===  Time= 3.65000E+02 [d]  Dt= 1.00000E+01 [d]"
# This tells you the simulation is at day 365 (year 1) with 10-day timesteps

# Check for convergence issues:
grep -i "WARNING\|ERROR\|CONVERGENCE\|CUT" run_scenario1.log | tail -20
```

### Common Issues and Fixes

**Issue: "Timestep cut" messages appearing frequently**

This means the solver is struggling to converge. Reduce the maximum timestep:

```
# In the input deck, change:
MAXIMUM_TIMESTEP_SIZE 30.d0 d
# To:
MAXIMUM_TIMESTEP_SIZE 10.d0 d
```

**Issue: "Negative concentration" or "Negative saturation" errors**

This usually means the chemistry is stiff. Add damping:

```
NEWTON_SOLVER
  RTOL 1.d-8
  ATOL 1.d-8
  STOL 1.d-30
  MAXIT 25
  MAX_ALLOW_REL_LIQ_PRES_CHANG_NI 1.d-1
/
```

**Issue: Database mineral not found**

Check that mineral names in your input deck match the database exactly:

```bash
grep "Forsterite" database/hanford.dat
# If not found, check available olivine names:
grep -i "oliv\|forst" database/hanford.dat
```

**Issue: Memory errors**

The 100×1×50 grid is moderate. If running on a laptop, reduce grid size:

```
NXYZ 50 1 25    # Half resolution — good for testing
```

---

## 6. Output Files and Post-Processing

### PFLOTRAN Output Structure

After a successful run, you'll find:

```
scenario1_dissolved/
├── basalt_wag_base.in
├── scenario1_dissolved_continuous-001.h5    # Snapshot at t=0.25 yr
├── scenario1_dissolved_continuous-002.h5    # Snapshot at t=0.5 yr
├── ...
├── scenario1_dissolved_continuous-obs.dat   # Observation point data
├── scenario1_dissolved_continuous-mas.dat   # Mass balance data
└── run_scenario1.log
```

### Reading HDF5 Output with Python

```bash
pip install h5py numpy matplotlib pandas
```

```python
#!/usr/bin/env python3
"""Read PFLOTRAN HDF5 output and extract key metrics."""

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Open a snapshot file
filename = "scenario1_dissolved_continuous-010.h5"  # Year 10
with h5py.File(filename, 'r') as f:
    # List available datasets
    print("Groups:", list(f.keys()))
    
    # Time information
    time_group = f['Time']
    print(f"Simulation time: {time_group.attrs['Time']} {time_group.attrs['Time Units']}")
    
    # Flow variables
    pressure = f['Liquid_Pressure'][:]          # Pa
    gas_sat = f['Gas_Saturation'][:]            # fraction
    temperature = f['Temperature'][:]            # °C
    porosity = f['Porosity'][:]                  # fraction
    perm_x = f['Permeability_X'][:]             # m²
    
    # Chemistry (if GIRT transport was used)
    if 'Total_H+' in f:
        pH = -np.log10(f['Total_H+'][:])
    
    # Mineral volume fractions
    if 'Volume_Fraction_Calcite' in f:
        calcite_vf = f['Volume_Fraction_Calcite'][:]
        magnesite_vf = f['Volume_Fraction_Magnesite'][:]
        siderite_vf = f['Volume_Fraction_Siderite'][:]
    
    print(f"Porosity range: {porosity.min():.4f} to {porosity.max():.4f}")
    print(f"Permeability range: {perm_x.min():.2e} to {perm_x.max():.2e} m²")
```

### Reading Observation Point Data

```python
import pandas as pd

# Observation data is in a text file
obs = pd.read_csv(
    "scenario1_dissolved_continuous-obs.dat",
    delim_whitespace=True,
    skiprows=1  # Skip header
)

print(obs.columns.tolist())
# Typical columns: Time, Liquid_Pressure, Gas_Saturation, Temperature,
# Total_H+, Total_Ca++, Total_Mg++, etc.
```

### Mass Balance Analysis

```python
# Mass balance file tracks total CO2 injected vs. mineralized
mass = pd.read_csv(
    "scenario1_dissolved_continuous-mas.dat",
    delim_whitespace=True,
    skiprows=1
)

# Calculate mineralization efficiency
co2_injected = mass['CO2(aq)_Injected']      # Total CO2 in (moles)
co2_mineral = mass['Calcite_Volume']          # Carbonate mineral volume
```

### Running the Full Post-Processing Pipeline

The `wag_optimizer.py` script I provided includes a semi-analytical model. Once you have real PFLOTRAN output, replace the semi-analytical engine with HDF5 readers:

```python
# In wag_optimizer.py, add this function to read real PFLOTRAN output:

def read_pflotran_results(scenario_dir, prefix):
    """Read PFLOTRAN HDF5 snapshots and observation data."""
    import h5py, glob, os
    
    # Find all snapshot files
    h5_files = sorted(glob.glob(os.path.join(scenario_dir, f"{prefix}-*.h5")))
    
    times = []
    porosities_nw = []
    porosities_ff = []
    perms_nw = []
    calcite_vols = []
    
    for h5f in h5_files:
        with h5py.File(h5f, 'r') as f:
            t = f['Time'].attrs.get('Time', 0)
            times.append(t)
            
            phi = f['Porosity'][:]
            k = f['Permeability_X'][:]
            
            # Near-wellbore = first 10 cells radially
            porosities_nw.append(np.mean(phi[:10, :, :]))
            porosities_ff.append(np.mean(phi[50:, :, :]))
            perms_nw.append(np.mean(k[:10, :, :]))
            
            if 'Volume_Fraction_Calcite' in f:
                calcite_vols.append(np.sum(f['Volume_Fraction_Calcite'][:]))
    
    return {
        'times': np.array(times),
        'porosity_nw': np.array(porosities_nw),
        'porosity_ff': np.array(porosities_ff),
        'perm_nw': np.array(perms_nw),
        'calcite_total': np.array(calcite_vols),
    }
```

---

## 7. Visualization with ParaView

PFLOTRAN HDF5 output can be visualized in ParaView:

```bash
# Install ParaView
sudo apt install paraview
# OR download from https://www.paraview.org/download/

# Open ParaView and use the XDMF reader:
# 1. File → Open → select the .h5 file
# 2. Choose "XDMF Reader" when prompted
# 3. Apply, then select variables to visualize
```

For batch visualization, use `pvpython`:

```python
# pvpython script for automated visualization
from paraview.simple import *
reader = XDMFReader(FileNames=['scenario1.h5'])
Show(reader)
ColorBy(reader, ('CELLS', 'Porosity'))
SaveScreenshot('porosity_year10.png')
```

---

## 8. Modifying Scenarios for Your Site

To adapt these simulations for a specific site (e.g., Deccan Traps, Paraná basalt):

### Change Rock Properties

```
MATERIAL_PROPERTY your_basalt
  POROSITY 0.08d0                    # Adjust to your formation
  PERMEABILITY
    PERM_ISO 5.d-14                  # Measured or estimated
  /
  ROCK_DENSITY 2950.d0              # Lab measurement
```

### Change Mineral Assemblage

Adjust the modal mineralogy in the `CONSTRAINT formation_water` block:

```
MINERALS
  Forsterite       0.10d0  1000.d0  # 10 vol% olivine (picritic basalt)
  Anorthite        0.20d0  2000.d0  # 20 vol% plagioclase
  Augite           0.25d0  1500.d0  # 25 vol% pyroxene
  Basalt_glass     0.10d0  5000.d0  # 10 vol% glass (crystalline basalt)
```

### Change Depth and Temperature

```
FLOW_CONDITION initial_formation
  LIQUID_PRESSURE 10.0d6            # 100 bar at 1000m
  TEMPERATURE 80.d0                  # Higher geothermal gradient
/
```

### Change Injection Rate

```
FLOW_CONDITION scenario1_injection
  RATE LIST
    0.0d0     2.0d-3     # Double injection rate
  /
```

---

## 9. Troubleshooting Quick Reference

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `pflotran: command not found` | Not in PATH | `export PATH=/path/to/pflotran/src/pflotran:$PATH` |
| PETSc version mismatch | PFLOTRAN requires specific PETSc | Check `pflotran/README` for compatible version |
| `Error: mineral not found in database` | Name mismatch | `grep -i "mineral_name" database/hanford.dat` |
| Simulation crashes at t=0 | Bad initial conditions | Start with coarser grid, verify pressure/temp |
| Very slow convergence | Chemistry too stiff | Reduce max timestep, increase Newton iterations |
| HDF5 files empty | Output times not reached | Check `TIME` block, ensure `FINAL_TIME` > first output |
| `MPI_ABORT` | Memory exceeded | Reduce grid or increase node memory |
| Unrealistic porosity | Precipitation rates too high | Reduce reactive surface area `A_s` by 10–100× |

---

## 10. Quick-Start Checklist

```
□  1. Install PETSc with HDF5 support
□  2. Build PFLOTRAN from source
□  3. Copy input decks to working directory
□  4. Copy thermodynamic database (hanford.dat)
□  5. Test with reduced grid (50×1×25) on 1 core
□  6. Verify output files are created
□  7. Scale to full grid (100×1×50) on multiple cores
□  8. Run all 6 scenarios
□  9. Post-process with Python (h5py + wag_optimizer.py)
□ 10. Visualize with ParaView or matplotlib
```

For questions about PFLOTRAN specifically, the documentation is at [https://www.pflotran.org](https://www.pflotran.org) and the user mailing list is active.
