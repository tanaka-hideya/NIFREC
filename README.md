# NIFREC: An Automated Local-Minimum Geometry Optimization Workflow with No Imaginary Frequencies

[![DOI](https://zenodo.org/badge/990414067.svg)](https://zenodo.org/badge/latestdoi/990414067)
[![tests](https://github.com/tanaka-hideya/NIFREC/actions/workflows/tests.yml/badge.svg)](https://github.com/tanaka-hideya/NIFREC/actions/workflows/tests.yml)

NIFREC is an automated workflow for local-minimum geometry optimization (no imaginary frequencies at the chosen level of theory) that includes an automated protocol to eliminate imaginary vibrational frequencies. NIFREC supports the sequential execution of conformer searches and quantum-chemical calculations across molecular datasets, automatically resolving imaginary frequencies at each stage of the workflow.  
NIFREC provides tools to generate conformers (RDKit); optimize geometries and analyze vibrational frequencies (xTB); run Gaussian optimization and frequency (opt+freq) jobs with robust imaginary-frequency remediation; and parse Gaussian results.

Note: the absence of imaginary frequencies establishes a local minimum at the chosen level of theory, not the global conformational minimum. By default, the Gaussian step starts from the lowest-energy xTB conformer; all xTB conformers without significant imaginary frequencies can be submitted instead (see Step 3).

Version: 2.0.0 (see [CHANGELOG.md](https://github.com/tanaka-hideya/NIFREC/blob/main/CHANGELOG.md) for the changes from v1.x)

## Authors
- Hideya Tanaka @ Nara Institute of Science and Technology (Author)
- Tomoyuki Miyao @ Nara Institute of Science and Technology (Contributor)

## Requirements
- See `environment.yml` for the complete dependency list and version details
- For the Gaussian step: a working installation of Gaussian 16. Either ensure `g16` is on your system `PATH`, or pass the full path to the Gaussian 16 executable via `nifrec-gaussian-optfreq`'s `--gcmd` argument.

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

### Installing xTB binary

Note: On macOS and Linux, xTB can be installed via conda, so the installation steps above are sufficient. If you want to use `nifrec-xtb` on Windows, follow the steps below. If you are using NIFREC without `nifrec-xtb` on Windows, the installation steps above are sufficient.

1. After downloading `environment.yml` in the steps above, remove only the xTB (`xtb`) entry from the YAML file, and then proceed with the same installation steps.

Alternatively, without creating an environment from the YAML file, you can install all Python dependencies except xTB into an existing environment via pip:

```bash
pip install "nifrec[chem] @ git+https://github.com/tanaka-hideya/NIFREC.git@main"
```

2. Download the xTB binary from the official [xTB GitHub repository](https://github.com/grimme-lab/xtb/releases/tag/v6.7.1). If `xtb` is not on your system `PATH`, pass the full path to the xTB executable via `nifrec-xtb`'s `--xcmd` argument.

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
- In this screening step, modes with |frequency| <= `--imagfreq-thres` (default: 5 cm^-1) are treated as numerical noise. In the Gaussian step, no imaginary frequency is tolerated (see Step 3).

Key outputs
- ./xtb/xtbopt_emin_xyz: XYZ files for the minimum-energy conformers (per molecule)
- ./xtb/xTB_stats_Emin.csv: summary of the minimum-energy conformer for each molecule

### Step 3 — Gaussian opt+freq with imaginary-frequency remediation

Requires Gaussian 16. The route section is built as: "#p theory-level opt freq=noraman". Do not include opt/freq in --theory-level.

```bash
nifrec-gaussian-optfreq --outfolder-gaussian gaussian_optfreq_PM6 --infolder-xtb xtb --suffix optfreq_PM6 --theory-level PM6 --nproc 8 --mem 32
```

To submit all xTB conformers without significant imaginary frequencies (not only the lowest-energy one), add `--infolder-xtb-xyz opt --infile xTB_stats_all.csv`.

Key outputs
- ./gaussian_optfreq_PM6/gaussian_gjf_optfreq_PM6, ./gaussian_optfreq_PM6/gaussian_log_optfreq_PM6, ./gaussian_optfreq_PM6/gaussian_chk_optfreq_PM6: artifacts from successful runs
- ./gaussian_optfreq_PM6/gaussian_imagf_optfreq_PM6: runs that retained imaginary frequencies (files are renamed with the suffixes _0, _1, and _2_k for Stage 0, Stage 1, and Stage 2 trial k)
- ./gaussian_optfreq_PM6/gaussian_working_optfreq_PM6: working directory (contains only failed cases after completion)
- ./gaussian_optfreq_PM6/gaussian_optfreq_PM6_stats.csv: summary CSV (see "Output columns of the Gaussian step")
- ./gaussian_optfreq_PM6/log_gaussian.txt: settings, NIFREC version, host name, and CPU

#### Imaginary-frequency remediation

Each molecule is processed in up to three stages. A structure is regarded as free of imaginary frequencies only if no vibrational frequency is negative; unlike the xTB screening, no tolerance threshold is applied, so that no case-by-case inspection of small imaginary frequencies is required.

- Stage 0: "opt freq" from the input structure.
- Stage 1: if imaginary frequencies are found, the optimization is restarted from the Stage 0 checkpoint file using the computed force constants (`opt=RCFC freq Guess=Read Geom=AllCheck`). Skipped with `--skip-stage1`.
- Stage 2: if imaginary frequencies persist, a normalized (unit-length, 3N-dimensional) Cartesian displacement vector is built from the imaginary modes of the preceding stage. Trial k (k = 0, 1, ..., `--max-repeat` - 1) displaces the structure of the preceding stage by `--base-disp` x (k + 1) angstrom along this vector and runs a fresh "opt freq" calculation (same route section as Stage 0). Every trial starts from the same structure with the same vector; the vector is not recomputed between trials. The loop stops as soon as no imaginary frequency remains.

Displacement vector
- Sign convention: the sign of a normal-mode displacement vector is arbitrary in the Gaussian output. NIFREC fixes it deterministically: the largest-magnitude Cartesian component of each mode vector (the first one in atom order in case of ties) is made positive.
- `--disp-vec sum` (default): component-wise sum of the sign-fixed vectors of all imaginary modes, normalized to unit length. `--disp-vec largest`: the sign-fixed vector of the most negative mode only. Both give the same vector when only one imaginary frequency is present.
- `--reverse-disp`: multiply the resulting vector by -1 (displace in the opposite direction).

Other options
- `--option-opt`: options added to the opt keyword in all stages (e.g., `MaxCycles=200`).
- `--option-opt-fc`: options used only in Stage 0 and Stage 2 because they must not be combined with RCFC (e.g., `CalcFC`).
- `--option-freq`: text appended to the freq keyword (default: `=noraman`).
- Gaussian's default optimization convergence criteria and numerical settings (e.g., the default integration grid for DFT) are used unless specified via `--theory-level`, `--option-opt`, or `--option-freq`.

#### Using Gaussian input files (.gjf/.com) instead of the xTB workflow

```bash
nifrec-gaussian-optfreq --outfolder-gaussian gaussian_from_gjf --infolder-gjf my_inputs --suffix b3lyp --theory-level "B3LYP/6-31G(d)"
```

- All .gjf and .com files in `--infolder-gjf` are processed in order of file name. The file name without the extension is the molecule identifier (index of the output CSV file) and must be unique and free of whitespace.
- The charge, multiplicity, and coordinates (Cartesian or Z-matrix, in angstrom) are read from each file with ASE (`ase.io.read(..., format='gaussian-in')`). The route section of the file is not used; it is built from the command-line options as usual. Input files that ASE cannot read (e.g., with freeze codes in the molecule specification) are reported with an error.
- Use this mode to set the charge and multiplicity explicitly (see "Charge and spin multiplicity").

#### Recalculating failed molecules

Molecules recorded with confid = 0 (abnormal termination, or imaginary frequencies remaining after Stage 2) can be recalculated, for example with different opt options, using the same input as the previous run and a new output folder:

```bash
nifrec-gaussian-optfreq --outfolder-gaussian gaussian_optfreq_PM6_recalc --infolder-xtb xtb --suffix optfreq_PM6_recalc --theory-level PM6 --option-opt-fc CalcFC --infile-gaussian-recalc gaussian_optfreq_PM6/gaussian_optfreq_PM6_stats.csv
```

The recalculation starts again from the input structures and applies the same Stage 0-2 procedure. The `fail_stage` column of the previous run shows where each failed molecule stopped.

#### Resuming an interrupted run

Molecules are processed sequentially from the top of the input CSV file, and the stats CSV file is updated after each molecule. To resume an interrupted run, delete the rows of the completed molecules from (a copy of) the input CSV file and run the step again with a new `--outfolder-gaussian`.

#### Charge and spin multiplicity

With `--infolder-xtb`, the charge is the formal charge of the RDKit molecule, and the spin multiplicity is the number of radical electrons in the RDKit molecule + 1 (i.e., the high-spin state is assumed; closed-shell molecules are singlets). This rule is not sufficient for all electronic states (e.g., open-shell singlets or low-spin states). For such systems, specify the charge and multiplicity explicitly in Gaussian input files and use `--infolder-gjf`. For unrestricted calculations, add `--no-homo-lumo` in Step 4.

#### Output columns of the Gaussian step

One row is written per molecule, and the stats CSV file is updated after each molecule; empty cells denote missing values. The same definitions are shown by `nifrec-gaussian-optfreq --help`. Stage 2 trial k (k = 0, 1, ..., `--max-repeat` - 1) displaces the structure of the preceding stage by `--base-disp` x (k + 1) angstrom along the normalized (unit-length) displacement vector. The final structure and the final energy of a Gaussian job are the last ones in its .log file.

| Column | Description |
| --- | --- |
| smiles, molid, total_energy_xTB | Copied from the input CSV file (total_energy_xTB: xTB total energy in hartree). With `--infolder-gjf`, smiles and total_energy_xTB are empty, and molid is the file name without the extension. |
| confid | Conformer ID from the input CSV file (1 with `--infolder-gjf`); 0 if no structure without imaginary frequencies was obtained. |
| charge, multiplicity | Charge and spin multiplicity used in the Gaussian jobs. |
| filepath | Name of the .log file of the successful Gaussian job. |
| success_stage | Stage (0, 1, or 2) at which a structure without imaginary frequencies was obtained. |
| success_disploop | Stage 2 trial k (0-based) at which a structure without imaginary frequencies was obtained; -1 if imaginary frequencies remained after `--max-repeat` trials. |
| fail_stage | Stage (0, 1, or 2) at which a Gaussian job terminated abnormally or its output could not be processed (empty otherwise). |
| n_imag_s0, n_imag_s1, n_imag_s2 | Number of imaginary (negative) frequencies after Stage 0, Stage 1, and the last Stage 2 trial. |
| imag_freqs_s0_per_cm, imag_freqs_s1_per_cm, imag_freqs_s2_per_cm | Imaginary frequencies (cm^-1) after the same stages, in ascending order and separated by ";" (empty if there are none). |
| rmsd_in_s0_angstrom | RMSD (angstrom) between the input structure and the final structure of Stage 0. |
| rmsd_s0_s1_angstrom, rmsd_s0_s2_angstrom | RMSD (angstrom) between the final structure of Stage 0 and that of Stage 1 or of the last Stage 2 trial. |
| dE_s0_s1_kJ_per_mol, dE_s0_s2_kJ_per_mol | Final energy of Stage 1 or of the last Stage 2 trial minus the final energy of Stage 0 (kJ/mol). |
| wall_time_seconds | Elapsed (wall-clock) time in seconds for the molecule, measured from the preparation of the Stage 0 input file to the end of the processing of its last stage, including the Gaussian jobs, the parsing of their .log files, and the file handling (also recorded for molecules that failed). |

Notes on the columns
- The columns n_imag, imag_freqs, rmsd, and dE of a stage are filled only if the Gaussian job of that stage (for Stage 2, the last trial) terminated normally and its .log file was parsed; they are empty if the stage was not run or failed.
- RMSDs are calculated over all atoms (same atom order, without mass weighting) after optimal superposition by translation and rotation with ASE (`ase.build.minimize_rotation_and_translation`); equivalent atoms are not permuted.
- Energies are the final energies parsed by cclib (the SCF energy, or the MP or CC energy if present); differences are converted to kJ/mol with `cclib.parser.utils.convertor`.
- wall_time_seconds depends on the hardware and the machine load; compare values only between runs under the same conditions. The host name and CPU are written to log_gaussian.txt.

### Step 4 — Parse Gaussian results

Parse Gaussian log files to extract energies and, for restricted methods, HOMO/LUMO orbital energies. Use a distinct output filename to avoid overwriting the Gaussian summary.

```bash
nifrec-gaussian-parse --infolder-gaussian gaussian_optfreq_PM6 --infolder-gaussian-log gaussian_optfreq_PM6/gaussian_log_optfreq_PM6 --infile gaussian_optfreq_PM6_stats.csv --outfile gaussian_optfreq_PM6_parse.csv
```

Key outputs
- ./gaussian_optfreq_PM6/gaussian_optfreq_PM6_parse.csv: summary CSV

Notes
- For unrestricted (UHF) calculations, append --no-homo-lumo to skip HOMO/LUMO extraction.
- All columns of the input CSV file (including the columns of Step 3) are carried over to the output CSV file.

## Testing

The test suite (pytest) covers the stage logic of the Gaussian step, the construction of the displacement vectors, failed-job handling, recalculation, the output records, and the parsing step. Gaussian is not required: the Gaussian runs are simulated in the tests. If xTB and MORFEUS are available, an end-to-end test of the RDKit, xTB, and Gaussian (simulated) steps is also run. The tests run automatically on GitHub Actions.

```bash
git clone https://github.com/tanaka-hideya/NIFREC.git
cd NIFREC
conda env create -f environment.yml
conda activate nifrec
pip install --no-deps -e .
pip install "pytest>=8"
python -m pytest
```

Optional: place real Gaussian "opt freq" log files in `tests/data/gaussian/` to check that they are parsed correctly by the installed cclib (`tests/test_gaussian_real_logs.py`).

## Tips and troubleshooting
- Unique identifiers: The index column specified by --idxcol must uniquely identify molecules; it is used in filenames and CSV indices.
- Parallelism: Heavy steps support parallel execution. Use --njobs to control the number of workers (negative values use CPU cores - 1).
- Output folders: Each step creates its output folder with the exact path you specify; if the folder already exists, creation may fail. Use a fresh path or remove the existing directory before rerunning.
- Dots in paths: Avoid `.` in `--outfolder-rdkit` (and in the current directory for relative paths). Per-conformer XYZ export may split the output path at the first dot and write files to an unexpected location.

## Citation
If you use NIFREC in your work, please cite it. See [CITATION.cff](https://github.com/tanaka-hideya/NIFREC/blob/main/CITATION.cff) in this repository.

## License
This project is licensed under the terms of the MIT license. See [LICENSE](https://github.com/tanaka-hideya/NIFREC/blob/main/LICENSE) for details.
