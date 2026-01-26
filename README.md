# NIFREC: An Automated Ground-State Geometry Optimization Workflow with No Imaginary Frequencies

[![DOI](https://zenodo.org/badge/990414067.svg)](https://zenodo.org/badge/latestdoi/990414067)

NIFREC is an automated workflow for ground-state geometry optimization that includes an automated protocol to eliminate imaginary vibrational frequencies. NIFREC supports the sequential execution of conformer searches and quantum-chemical calculations across molecular datasets, automatically resolving imaginary frequencies at each stage of the workflow.  
NIFREC provides tools to generate conformers (RDKit); optimize geometries and analyze vibrational frequencies (xTB); run Gaussian optimization and frequency (opt+freq) jobs with robust imaginary-frequency remediation; and parse Gaussian results.

Version: 1.1.1

## Authors
- Hideya Tanaka @ Nara Institute of Science and Technology (Author)
- Tomoyuki Miyao @ Nara Institute of Science and Technology (Contributor)

## Requirements
- Python 3.12 or later
- See `environment.yml` for the complete dependency list and version details
- For the Gaussian step: a working installation of Gaussian 16. Either ensure `g16` is on your system `PATH`, or pass the full path to the Gaussian 16 executable via a CLI argument.

## Installation

1. Create the conda environment defined in the repository and activate it:

```bash
curl -L https://raw.githubusercontent.com/tanaka-hideya/NIFREC/main/environment.yml -o environment.yml
conda env create -f environment.yml
conda activate nifrec
```

2. Install the package only (skipping dependencies already provided by the environment):

```bash
pip install "nifrec @ git+https://github.com/tanaka-hideya/NIFREC.git@main" --no-deps
```

Explanation: --no-deps installs the nifrec package itself without pulling dependencies from PyPI.

No PYTHONPATH setup is required. Command-line entry points are provided via [project.scripts].

## Command-line tools
Installed scripts (see [pyproject.toml](https://github.com/tanaka-hideya/NIFREC/blob/main/pyproject.toml)):
- `nifrec-rdkit` — RDKit-based conformer generation
- `nifrec-xtb` — xTB geometry optimization and vibrational analysis with automatic handling of small imaginary frequencies
- `nifrec-gaussian-optfreq` — Gaussian opt+freq runs with robust imaginary-frequency remediation
- `nifrec-gaussian-parse` — Gaussian log-file parser with consolidated summary CSV output

Tip: Append --help to any command for full options, e.g.:

```bash
nifrec-rdkit --help
```

Then, the following output is shown:

```text
usage: nifrec_rdkit [-h] --outfolder-rdkit OUTFOLDER_RDKIT --infile INFILE [--smicol SMICOL] [--idxcol IDXCOL] [--outfile OUTFILE] [--nconfs NCONFS]
                    [--rmsd-thres RMSD_THRES] [--njobs NJOBS] [--backend BACKEND] [--random-seed RANDOM_SEED] [--force-field FORCE_FIELD]

Conformer generator using RDKit

options:
  -h, --help            show this help message and exit
  --outfolder-rdkit OUTFOLDER_RDKIT
                        Output folder to write RDKit results (XYZ and SDF files, logs). Accepts absolute or relative paths; '~' is expanded. The folder is
                        created.
  --infile INFILE       Path to the input CSV file. Accepts absolute or relative paths; '~' is expanded. Must contain a SMILES column specified by --smicol.
  --smicol SMICOL       Name of the column in the input CSV that contains SMILES strings (used for structure generation). (default: smiles)
  --idxcol IDXCOL       Zero-based index of the column in the input CSV to use as the unique molecule identifier (DataFrame index). Identifiers must be unique
                        per molecule and are used consistently across all outputs: the output CSV (--outfile) and the filenames of 3D structure files (XYZ/SDF).
                        Non-unique values may cause file overwrites and inconsistent results. (default: 0)
  --outfile OUTFILE     Name of the output CSV file to write summary (saved under --outfolder-rdkit). Note: If the SMILES column specified by --smicol is not
                        named 'smiles', it will be renamed to 'smiles' in the output CSV (unless a 'smiles' column already exists). Canonical smiles are stored
                        in the 'smiles' column. (default: rdkit_stats.csv)
  --nconfs NCONFS       Maximum number of RDKit conformers to generate per molecule. (default: 20)
  --rmsd-thres RMSD_THRES
                        RMSD pruning threshold between generated conformers. (default: 1)
  --njobs NJOBS         Number of parallel workers. If <= 0, uses (CPU cores - 1). (default: -1)
  --backend BACKEND     Parallel backend for joblib. One of: 'loky' (default), 'multiprocessing', 'threading'. (default: loky)
  --random-seed RANDOM_SEED
                        Random seed for reproducibility. (default: 42)
  --force-field FORCE_FIELD
                        Force field to use for RDKit conformer generation (MMFF94s, MMFF94, or UFF). (default: MMFF94s)
```

If any of the commands fail to run, try adjusting your setuptools version.

## Quick start with the sample dataset
This section demonstrates the full pipeline in the current directory. Before starting, place the sample CSV file in the working directory:

  ```bash
  curl -L https://raw.githubusercontent.com/tanaka-hideya/NIFREC/main/data/sample.csv -o sample.csv
  ```

The sample CSV file contains two columns: "name" and "smi".

### Step 1 — RDKit conformer generation

Generate 3D conformers from SMILES. The sample file uses the "smi" column for SMILES and the first column (index 0) as a unique identifier for file naming.

Implementation note: Conformer generation leverages the MORFEUS library (`ConformerEnsemble`) on top of RDKit for ensemble creation, RMSD-based pruning, and sorting before export.

```bash
nifrec-rdkit --outfolder-rdkit rdkit --infile sample.csv --smicol smi
```

Key outputs
- ./rdkit/xyz: conformer XYZ files named rdkit_idx_confid.xyz
- ./rdkit/rdkit_stats.csv: summary CSV

### Step 2 — xTB optimization and vibrational analysis

Optimize each RDKit conformer with xTB, iteratively addressing residual imaginary frequencies, and record the lowest-energy conformer for each molecule.

```bash
nifrec-xtb --outfolder-xtb xtb --infolder-rdkit rdkit
```

Internal xTB command (GFN2-xTB)

xtb input.xyz --ohess --chrg formal_charge --json

Details
- The formal charge is taken from the RDKit molecule.
- Frequencies are obtained via --ohess and parsed from the JSON output produced by --json.
- When significant imaginary modes remain, the same command is re-run on the distorted geometry file written by xTB (xtbhess.xyz) until convergence or the iteration limit is reached.

Key outputs
- ./xtb/xtbopt_emin_xyz: XYZ files for the minimum-energy conformers (per molecule)
- ./xtb/xTB_stats_Emin.csv: summary of the minimum-energy conformer for each molecule

### Step 3 — Gaussian opt+freq with imaginary-frequency remediation

Requires Gaussian 16. The route section is built as: "#p theory-level opt freq=noraman". Do not include opt/freq in --theory-level.

```bash
nifrec-gaussian-optfreq --outfolder-gaussian gaussian_optfreq_PM6 --infolder-xtb xtb --suffix optfreq_PM6 --theory-level PM6 --nproc 8 --mem 32
```

Key outputs
- ./gaussian_optfreq_PM6/gaussian_gjf_optfreq_PM6, ./gaussian_optfreq_PM6/gaussian_log_optfreq_PM6, ./gaussian_optfreq_PM6/gaussian_chk_optfreq_PM6: artifacts from successful runs
- ./gaussian_optfreq_PM6/gaussian_imagf_optfreq_PM6: runs that retained imaginary frequencies (files are renamed with suffixes)
- ./gaussian_optfreq_PM6/gaussian_working_optfreq_PM6: working directory (contains only failed cases after completion)
- ./gaussian_optfreq_PM6/gaussian_optfreq_PM6_stats.csv: summary CSV

### Step 4 — Parse Gaussian results

Parse Gaussian log files to extract energies and, for restricted methods, HOMO/LUMO orbital energies. Use a distinct output filename to avoid overwriting the Gaussian summary.

```bash
nifrec-gaussian-parse --infolder-gaussian gaussian_optfreq_PM6 --infolder-gaussian-log gaussian_optfreq_PM6/gaussian_log_optfreq_PM6 --infile gaussian_optfreq_PM6_stats.csv --outfile gaussian_optfreq_PM6_parse.csv
```

Key outputs
- ./gaussian_optfreq_PM6/gaussian_optfreq_PM6_parse.csv: summary CSV

Notes
- For unrestricted (UHF) calculations, append --no-homo-lumo to skip HOMO/LUMO extraction.

## Tips and troubleshooting
- Unique identifiers: The index column specified by --idxcol must uniquely identify molecules; it is used in filenames and CSV indices.
- Parallelism: Heavy steps support parallel execution. Use --njobs to control the number of workers (negative values use CPU cores - 1).
- Output folders: Each step creates its output folder with the exact path you specify; if the folder already exists, creation may fail. Use a fresh path or remove the existing directory before rerunning.

## Citation
If you use NIFREC in your work, please cite it. See [CITATION.cff](https://github.com/tanaka-hideya/NIFREC/blob/main/CITATION.cff) in this repository.

## License
This project is licensed under the terms of the MIT license. See [LICENSE](https://github.com/tanaka-hideya/NIFREC/blob/main/LICENSE) for details.
