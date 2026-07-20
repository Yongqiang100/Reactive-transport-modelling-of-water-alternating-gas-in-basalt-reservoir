# PFLOTRAN WAG Injection Simulation — Installation & Execution Guide

## Overview

This guide walks you through installing PFLOTRAN, running the 6 WAG injection scenarios for CO₂ mineralization in basalt, and post-processing the results. PFLOTRAN is a compiled Fortran/MPI reactive transport code that requires PETSc as its numerical backend.

**Estimated setup time:** 1–3 hours (mostly PETSc compilation)
**Estimated run time per scenario:** 2–48 hours depending on hardware and grid resolution

---

## 1. Prerequisites

You need a Linux or macOS system with:

- **Fortran compiler**: gfortran ≥ 9.0 (or Intel ifort/ifx)
- **C/C++ compilers**: gcc/g++ ≥ 9.0
- **Python** ≥ 3.6 (for PETSc configuration and post-processing)
- **Git**
- **Make / CMake**
- At least **8 GB RAM** and **20 GB disk space**

### Check your compilers

```bash
gcc --version
gfortran --version
python3 --version
git --version
```

On Ubuntu/Debian, install missing dependencies:

```bash
sudo apt update
sudo apt install -y build-essential gfortran cmake git python3 python3-pip \
    liblapack-dev libblas-dev zlib1g-dev
```

On macOS (with Homebrew):

```bash
brew install gcc cmake git python3
```

---

## 2. Install PETSc

PFLOTRAN depends on PETSc, which also downloads and builds MPI, HDF5, and other dependencies automatically.

### 2.1 Clone PETSc

```bash
# Choose a working directory
export PFLOTRAN_INSTALL=$HOME/pflotran_install
mkdir -p $PFLOTRAN_INSTALL && cd $PFLOTRAN_INSTALL

# Clone PETSc — use a version compatible with PFLOTRAN
git clone https://gitlab.com/petsc/petsc.git petsc
cd petsc

# PFLOTRAN v6.0 needs PETSc v3.21.x; master branch needs v3.24.x
# Check https://www.pflotran.org for the latest compatibility
git checkout v3.21.5
```

### 2.2 Configure PETSc

This single command downloads and compiles MPI, HDF5, BLAS/LAPACK, Metis, and ParMetis:

```bash
export PETSC_DIR=$PWD
export PETSC_ARCH=pflotran-opt

./configure \
  --PETSC_ARCH=$PETSC_ARCH \
  --COPTFLAGS='-O3' \
  --CXXOPTFLAGS='-O3' \
  --FOPTFLAGS='-O3 -Wno-unused-function' \
  --with-debugging=no \
  --download-mpich=yes \
  --download-hdf5=yes \
  --download-hdf5-fortran-bindings=yes \
  --download-fblaslapack=yes \
  --download-metis=yes \
  --download-parmetis=yes \
  --download-hdf5-configure-arguments="--with-zlib=yes"
```

**This takes 20–60 minutes.** If it succeeds, you'll see:

```
xxx=========================================================================xxx
 Configure stage complete. Now build PETSc libraries with:
   make PETSC_DIR=... PETSC_ARCH=... all
xxx=========================================================================xxx
```

### 2.3 Build PETSc

```bash
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all -j4
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH check
```

### 2.4 Set environment variables permanently

Add to your `~/.bashrc` or `~/.zshrc`:

```bash
export PETSC_DIR=$HOME/pflotran_install/petsc
export PETSC_ARCH=pflotran-opt
export PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH
```

Then reload: `source ~/.bashrc`

---

## 3. Install PFLOTRAN

### 3.1 Clone and build

```bash
cd $PFLOTRAN_INSTALL

# Clone PFLOTRAN
git clone https://bitbucket.org/pflotran/pflotran.git
cd pflotran

# Checkout a stable release (match your PETSc version)
git checkout v6.0

# Compile
cd src/pflotran
make -j4 pflotran
```

If successful, a `pflotran` executable appears in the current directory.

### 3.2 Verify installation

```bash
# Add to PATH
export PFLOTRAN_DIR=$PFLOTRAN_INSTALL/pflotran
export PATH=$PFLOTRAN_DIR/src/pflotran:$PATH

# Run a quick regression test
cd $PFLOTRAN_DIR/regression_tests/default/543
mpirun -n 1 pflotran -pflotranin 543_flow.in
```

If it runs without errors, PFLOTRAN is ready.

---

## 4. Prepare the WAG Simulation Files

### 4.1 Set up the project directory

```bash
mkdir -p $HOME/wag_basalt_mineralization
cd $HOME/wag_basalt_mineralization

# Copy all input decks from the provided files
# (assuming you downloaded them to ~/Downloads/)
cp ~/Downloads/pflotran_input_decks/*.in ./

# Copy the thermodynamic database
cp $PFLOTRAN_DIR/database/hanford.dat ./
```

### 4.2 Verify the database path

Each input deck references `DATABASE hanford.dat`. Make sure `hanford.dat` is in the same directory as the `.in` files, or update the path:

```bash
# Check that hanford.dat contains the minerals we need
grep -i "Forsterite\|Calcite\|Magnesite\|Ankerite\|Siderite" hanford.dat | head -20
```

### 4.3 Important: Input deck modifications before running

The provided input decks are structured as a base deck plus scenario overrides. For PFLOTRAN to run each scenario, you need to **merge the base deck with each scenario file** into a single standalone `.in` file. Here's how:

```bash
# Create a Python script to merge input decks
cat > merge_decks.py << 'EOF'
#!/usr/bin/env python3
"""Merge base input deck with scenario overrides."""
import sys

def merge(base_file, scenario_file, output_file):
    with open(base_file) as f:
        base = f.read()
    with open(scenario_file) as f:
        scenario = f.read()

    # Remove the INCLUDE line from scenario
    lines = scenario.split('\n')
    lines = [l for l in lines if not l.strip().startswith(':INCLUDE')]

    # Extract override blocks from scenario
    override_text = '\n'.join(lines)

    # Write merged file: base + overrides appended
    # (PFLOTRAN uses last-defined-wins for duplicate blocks)
    with open(output_file, 'w') as f:
        f.write(base)
        f.write('\n\n: === SCENARIO OVERRIDES ===\n')
        f.write(override_text)

    print(f"Created: {output_file}")

if __name__ == '__main__':
    merge(sys.argv[1], sys.argv[2], sys.argv[3])
EOF

# Merge each scenario
python3 merge_decks.py basalt_wag_base.in scenario1_dissolved_continuous.in run_S1.in
python3 merge_decks.py basalt_wag_base.in scenario2_scco2_continuous.in run_S2.in
python3 merge_decks.py basalt_wag_base.in scenario3_wag_6month.in run_S3.in
python3 merge_decks.py basalt_wag_base.in scenario4_wag_3month.in run_S4.in
python3 merge_decks.py basalt_wag_base.in scenario5_swag_continuous.in run_S5.in
python3 merge_decks.py basalt_wag_base.in scenario6_adaptive_wag.in run_S6.in
```

---

## 5. Run the Simulations

### 5.1 Single scenario (local machine)

```bash
cd $HOME/wag_basalt_mineralization

# Run Scenario 1 (dissolved CO2 baseline) on 4 cores
mpirun -n 4 pflotran -pflotranin run_S1.in

# Output files appear as:
#   scenario1_dissolved_continuous-###.h5  (HDF5 snapshots)
#   scenario1_dissolved_continuous-obs-0.tec  (observation point data)
#   scenario1_dissolved_continuous-mas.dat  (mass balance)
```

### 5.2 Run all 6 scenarios sequentially

```bash
#!/bin/bash
# run_all.sh — Run all WAG scenarios sequentially
NPROCS=4

for i in 1 2 3 4 5 6; do
    echo "=========================================="
    echo "  Running Scenario S${i}"
    echo "=========================================="
    mpirun -n $NPROCS pflotran -pflotranin run_S${i}.in \
        2>&1 | tee log_S${i}.txt
    echo "Scenario S${i} completed at $(date)"
    echo ""
done

echo "All scenarios complete."
```

Make executable and run:

```bash
chmod +x run_all.sh
./run_all.sh
```

### 5.3 Run on an HPC cluster (SLURM)

For clusters with SLURM job scheduling:

```bash
cat > submit_wag.slurm << 'SLURM'
#!/bin/bash
#SBATCH --job-name=wag_mineralization
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --output=pflotran_%j.out
#SBATCH --error=pflotran_%j.err

# Load modules (adjust for your cluster)
module load gcc/11.2.0
module load openmpi/4.1.4
module load hdf5/1.12.2

# Set environment
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=pflotran-opt
export PFLOTRAN_EXE=/path/to/pflotran

# Run each scenario
for SCENARIO in run_S1.in run_S2.in run_S3.in run_S4.in run_S5.in run_S6.in; do
    echo "Running $SCENARIO on $SLURM_NTASKS cores..."
    srun $PFLOTRAN_EXE -pflotranin $SCENARIO
    echo "$SCENARIO completed at $(date)"
done
SLURM

sbatch submit_wag.slurm
```

### 5.4 PBS/Torque clusters

```bash
cat > submit_wag.pbs << 'PBS'
#!/bin/bash
#PBS -N wag_mineralization
#PBS -l nodes=2:ppn=24
#PBS -l walltime=48:00:00
#PBS -q standard

cd $PBS_O_WORKDIR

module load gcc openmpi hdf5

NPROCS=$(wc -l < $PBS_NODEFILE)

for SCENARIO in run_S1.in run_S2.in run_S3.in run_S4.in run_S5.in run_S6.in; do
    mpirun -np $NPROCS pflotran -pflotranin $SCENARIO
done
PBS

qsub submit_wag.pbs
```

---

## 6. Monitor Progress

While PFLOTRAN is running:

```bash
# Watch the log output in real time
tail -f log_S1.txt

# Check for convergence issues
grep -i "WARNING\|ERROR\|CONVERGE\|CUT" log_S1.txt

# Check timestep progress
grep "Step" log_S1.txt | tail -5

# Monitor HDF5 output file sizes (growing = still running)
watch -n 10 'ls -lh *.h5'
```

### Common issues and fixes

| Problem | Symptom | Fix |
|---------|---------|-----|
| Timestep cuts | `TS_CUT` messages | Reduce `MAX_SATURATION_CHANGE` or initial timestep |
| Non-convergence | `NEWTON_ITERATION_EXCEEDED` | Reduce max timestep size; check mineral kinetics |
| Missing species in DB | `ERROR: Species not found` | Check `hanford.dat` has all minerals listed in CHEMISTRY block |
| Memory error | `SIGKILL` or `OOM` | Reduce grid resolution or use more nodes |
| Slow convergence | Very small timesteps | Simplify mineral kinetics (fewer minerals) for initial testing |

---

## 7. Post-Processing Results

### 7.1 Install Python dependencies

```bash
pip install numpy h5py matplotlib pandas scipy
```

### 7.2 Extract data from HDF5 output

PFLOTRAN writes results as HDF5 files. Here's how to read them:

```python
#!/usr/bin/env python3
"""Read PFLOTRAN HDF5 output and extract key metrics."""
import h5py
import numpy as np
import matplotlib.pyplot as plt

def read_pflotran_h5(filename):
    """Read a PFLOTRAN HDF5 snapshot file."""
    data = {}
    with h5py.File(filename, 'r') as f:
        # List available time groups
        print("Available groups:")
        for key in f.keys():
            print(f"  {key}")

        # Each time snapshot is a group like 'Time:  1.00000E+00 y'
        for time_key in sorted(f.keys()):
            if 'Time' not in time_key:
                continue

            group = f[time_key]
            t_data = {}

            # Read available datasets
            for dset_name in group.keys():
                t_data[dset_name] = group[dset_name][:]

            data[time_key] = t_data

    return data


def extract_observation_data(obs_file):
    """Read PFLOTRAN observation point .tec file."""
    # Tecplot format: header line, then columns of data
    with open(obs_file, 'r') as f:
        header = f.readline().strip()
        # Parse variable names from header
        variables = header.replace('"', '').split(',')

    data = np.loadtxt(obs_file, skiprows=1, delimiter=',')
    return variables, data


def plot_mineralization_comparison(scenario_files):
    """Plot cumulative CO2 mineralized across scenarios."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = ['#2563eb', '#dc2626', '#16a34a', '#9333ea', '#ea580c', '#0891b2']
    labels = ['S1: Dissolved', 'S2: scCO2', 'S3: WAG-6mo',
              'S4: WAG-3mo', 'S5: SWAG', 'S6: Adaptive']

    for i, (sfile, label, color) in enumerate(
            zip(scenario_files, labels, colors)):
        try:
            variables, data = extract_observation_data(sfile)
            time_col = 0  # First column is usually time
            time = data[:, time_col]

            # Find relevant columns (names vary by PFLOTRAN version)
            # Typical: pH, mineral volume fractions, porosity
            for j, var in enumerate(variables):
                if 'pH' in var:
                    axes[0, 0].plot(time, data[:, j], color=color, label=label)
                    axes[0, 0].set_ylabel('pH')
                if 'Porosity' in var or 'porosity' in var:
                    axes[0, 1].plot(time, data[:, j], color=color, label=label)
                    axes[0, 1].set_ylabel('Porosity')
                if 'Calcite' in var and 'VF' in var:
                    axes[1, 0].plot(time, data[:, j], color=color, label=label)
                    axes[1, 0].set_ylabel('Calcite Volume Fraction')
                if 'Permeability' in var or 'perm' in var.lower():
                    axes[1, 1].plot(time, data[:, j], color=color, label=label)
                    axes[1, 1].set_ylabel('Permeability (m²)')
                    axes[1, 1].set_yscale('log')
        except FileNotFoundError:
            print(f"Warning: {sfile} not found, skipping")

    for ax in axes.flat:
        ax.set_xlabel('Time (years)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('WAG Injection Strategy Comparison — PFLOTRAN Results',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('wag_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Saved: wag_comparison.png")


if __name__ == '__main__':
    # Update these filenames to match your PFLOTRAN output
    obs_files = [
        'scenario1_dissolved_continuous-obs-0.tec',
        'scenario2_scco2_continuous-obs-0.tec',
        'scenario3_wag_6month-obs-0.tec',
        'scenario4_wag_3month-obs-0.tec',
        'scenario5_swag_continuous-obs-0.tec',
        'scenario6_adaptive_wag-obs-0.tec',
    ]

    plot_mineralization_comparison(obs_files)
```

### 7.3 Use the provided analysis framework

The `wag_optimizer.py` script (provided with the input decks) can also process PFLOTRAN output. After running the simulations, update the script to read from HDF5 files instead of the semi-analytical model.

---

## 8. Troubleshooting

### PETSc won't configure

```bash
# Common fix: explicitly specify compilers
./configure \
  --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
  [... rest of options ...]
```

### PFLOTRAN won't compile

```bash
# Ensure environment is set
echo $PETSC_DIR    # Should show your PETSc directory
echo $PETSC_ARCH   # Should show your architecture name

# Clean and retry
make clean
make -j4 pflotran
```

### Database species not found

If PFLOTRAN reports minerals not found in `hanford.dat`, it's because some mineral names in the input deck may differ from the database spellings. Check:

```bash
grep -i "forsterite" hanford.dat
grep -i "basalt_glass" hanford.dat
```

The `hanford.dat` database may not include `Basalt_glass` as a named mineral. In that case, replace it with a proxy (e.g., `SiO2(am)` with adjusted kinetics) or use a custom database. The CarbFix project developed the `carbfix.dat` database specifically for basalt carbonation — it can be used with PHREEQC but would need conversion for PFLOTRAN format.

### Simulation runs very slowly

- Start with a **coarser grid** (e.g., 50×1×25 instead of 100×1×50) for testing
- Reduce the number of **secondary species** (keep only the most important ones)
- Increase `MAX_SATURATION_CHANGE` from 0.05 to 0.1
- Use fewer observation output times initially

---

## 9. Expected Output Files

After each scenario completes, you should see:

```
scenario*-001.h5     # HDF5 snapshot at year 0.25
scenario*-002.h5     # HDF5 snapshot at year 0.5
...
scenario*-011.h5     # HDF5 snapshot at year 30
scenario*-obs-0.tec  # Observation point time series (Tecplot format)
scenario*-mas.dat    # Mass balance over time
pflotran.out         # Main log file with convergence info
```

---

## 10. Quick-Start Summary

```bash
# 1. Install PETSc (one-time, ~1 hour)
git clone https://gitlab.com/petsc/petsc.git && cd petsc
git checkout v3.21.5
./configure --download-mpich=yes --download-hdf5=yes \
  --download-hdf5-fortran-bindings=yes --download-fblaslapack=yes \
  --download-metis=yes --download-parmetis=yes --with-debugging=no
make all && make check

# 2. Install PFLOTRAN (one-time, ~10 min)
git clone https://bitbucket.org/pflotran/pflotran.git
cd pflotran && git checkout v6.0
cd src/pflotran && make -j4 pflotran

# 3. Set up and run
cd ~/wag_basalt_mineralization
cp /path/to/hanford.dat .
python3 merge_decks.py basalt_wag_base.in scenario1_dissolved_continuous.in run_S1.in
mpirun -n 4 pflotran -pflotranin run_S1.in

# 4. Post-process
python3 wag_optimizer.py  # Semi-analytical comparison
python3 plot_results.py   # Plot PFLOTRAN HDF5 output
```

---

## References

- PFLOTRAN documentation: https://www.pflotran.org/documentation/
- PETSc installation: https://petsc.org/release/install/
- PFLOTRAN on Linux: https://www.pflotran.org/documentation/user_guide/how_to/installation/linux.html
- hanford.dat database: located in `pflotran/database/` after cloning
- CarbFix database (for PHREEQC): Voigt et al. (2018), Energy Procedia
