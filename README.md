# NIFREC: An Automated Ground-State Structure Optimization Workflow with No Imaginary Frequencies

NIFREC is an automated geometry optimization workflow for molecular ground states that incorporates an automated protocol for the elimination of imaginary frequencies. NIFREC enables sequential execution of conformer search and quantum chemical calculations across molecular datasets, automatically resolving imaginary frequencies at each stage of the workflow.  
NIFREC provides tools to generate conformers (RDKit), optimize geometries and analyze frequencies (xTB), run Gaussian opt+freq with robust imaginary-frequency remediation, and parse Gaussian results.

Version: 1.0.0

## Authors
- Hideya Tanaka @ Nara Institute of Science and Technology (Author)
- Tomoyuki Miyao @ Nara Institute of Science and Technology (Contributor)

## Requirements
- Python 3.12+
- See `environment.yml` for full list and version details.
- For Gaussian step: a working Gaussian16 installation (g16 available on PATH)

## Installation

1) Create environment from the repository and activate it:

```bash
curl -L https://raw.githubusercontent.com/tanaka-hideya/NIFREC/main/environment.yml -o environment.yml
conda env create -f environment.yml
conda activate nifrec
```

2) Install the package only (skip dependencies already provided by the environment):

```bash
pip install "nifrec @ git+https://github.com/tanaka-hideya/NIFREC.git@main" --no-deps
```

Explanation: --no-deps installs the nifrec package itself without pulling dependencies from PyPI.

No PYTHONPATH setup is required. Command-line entry points are provided via [project.scripts].

## Command-line tools
Installed scripts (see [pyproject.toml](https://github.com/tanaka-hideya/NIFREC/blob/main/pyproject.toml)):
- nifrec-rdkit — RDKit conformer generation
- nifrec-xtb — xTB optimization and vibrational analysis with automatic handling of small imaginary frequencies
- nifrec-gaussian-optfreq — Gaussian opt+freq runs with robust imaginary-frequency remediation
- nifrec-gaussian-parse — Parser for Gaussian logs and summary CSV aggregation

Tip: Append --help to any command for full options, e.g.:

```bash
nifrec-rdkit --help
```

Key nifrec-rdkit options (summary)
- --outfolder-rdkit PATH (required): output folder; created if missing
- --infile FILE (required): input CSV file containing SMILES
- --smicol NAME: SMILES column name (default: smiles)
- --idxcol N: zero-based index column to use as unique identifier (default: 0)

If any of the commands below fail to run, try adjusting your setuptools version.

## Quick start with the sample dataset
This section demonstrates the full pipeline in the current directory. Before starting, place the sample CSV in the working directory:

  ```bash
  curl -L https://raw.githubusercontent.com/tanaka-hideya/NIFREC/main/data/sample.csv -o sample.csv
  ```

The sample CSV uses columns: name, smi.

### Step 1 — RDKit conformer generation

Generate 3D conformers from SMILES. The sample file uses column "smi" for SMILES and the first column (index 0) as a unique identifier used for file naming.

Implementation note: Conformer generation leverages the MORFEUS library (ConformerEnsemble) on top of RDKit for ensemble creation, RMSD-based pruning, and sorting before export.

```bash
nifrec-rdkit --outfolder-rdkit rdkit --infile sample.csv --smicol smi
```

Outputs
- ./rdkit/xyz: conformer XYZ files named rdkit_idx_confid.xyz
- ./rdkit/sdf: 2D connectivity SDF files
- ./rdkit/rdkit_stats.csv: summary CSV

### Step 2 — xTB optimization and vibrational analysis

Optimize each RDKit conformer with xTB, iteratively handling small imaginary frequencies, and record the lowest-energy conformer per molecule.

```bash
nifrec-xtb --outfolder-xtb xtb --infolder-rdkit rdkit
```

xTB command line (default method: GFN2-xTB)

xtb input.xyz --ohess --chrg formal_charge --json

Details
- The formal charge is taken from the RDKit molecule.
- Frequencies are obtained via --ohess and parsed from the JSON output produced by --json.
- When significant imaginary modes remain, the same command is re-run on the distorted geometry file written by xTB (xtbhess.xyz) until convergence or the iteration limit is reached.

Outputs
- ./xtb/opt, ./xtb/imag_freq, ./xtb/working, ./xtb/worker: run artifacts and logs
- ./xtb/xtbopt_emin_xyz: XYZ files for the minimum-energy conformers (per molecule)
- ./xtb/xTB_stats_Emin.csv: per-molecule minima
- ./xtb/xTB_stats_all.csv: all conformers

### Step 3 — Gaussian opt+freq with imaginary-frequency remediation

Requires Gaussian16 (g16) in PATH. The route section is built as: "#p theory-level opt freq=noraman". Do not include opt/freq in --theory-level.

```bash
nifrec-gaussian-optfreq --outfolder-gaussian gaussian_optfreq_M062X_Def2TZVP --infolder-xtb xtb --suffix optfreq_M062X_Def2TZVP --theory-level M062X/Def2TZVP --nproc 8 --mem 32
```

Outputs
- ./gaussian_optfreq_M062X_Def2TZVP/gaussian_gjf_optfreq_M062X_Def2TZVP, ./gaussian_optfreq_M062X_Def2TZVP/gaussian_log_optfreq_M062X_Def2TZVP, ./gaussian_optfreq_M062X_Def2TZVP/gaussian_chk_optfreq_M062X_Def2TZVP: success artifacts
- ./gaussian_optfreq_M062X_Def2TZVP/gaussian_imagf_optfreq_M062X_Def2TZVP: runs that retained imaginary frequencies (files renamed with suffixes)
- ./gaussian_optfreq_M062X_Def2TZVP/gaussian_working_optfreq_M062X_Def2TZVP: working directory (contains only failed or still-problematic cases after completion)
- ./gaussian_optfreq_M062X_Def2TZVP/gaussian_optfreq_M062X_Def2TZVP_stats.csv: run summary

### Step 4 — Parse Gaussian results

Parse Gaussian logs to extract energies and (for restricted methods) HOMO/LUMO. Use a different output filename to avoid overwriting the Gaussian summary.

```bash
nifrec-gaussian-parse --infolder-gaussian gaussian_optfreq_M062X_Def2TZVP --infolder-gaussian-log gaussian_optfreq_M062X_Def2TZVP/gaussian_log_optfreq_M062X_Def2TZVP --infile gaussian_optfreq_M062X_Def2TZVP_stats.csv --outfile gaussian_optfreq_M062X_Def2TZVP_parse.csv
```

Notes
- For unrestricted (UHF) logs, append --no-homo-lumo to skip HOMO/LUMO extraction.

## Tips and troubleshooting
- Unique identifiers: The index column specified by --idxcol must uniquely identify molecules; it is used in filenames and CSV indices.
- Parallelism: Heavy steps support parallel execution. Use --njobs to control the number of workers (negative values use CPU cores - 1).
- Output folders: Each step creates its output folder with the exact path you specify; if the folder already exists, creation may fail. Use a fresh path or remove the existing directory before rerunning.

## Citation
If you use NIFREC in your work, please cite it. See [CITATION.cff](https://github.com/tanaka-hideya/NIFREC/blob/main/CITATION.cff) in this repository.

## License
This project is licensed under the terms of the MIT license. See [LICENSE](https://github.com/tanaka-hideya/NIFREC/blob/main/LICENSE) for details.



