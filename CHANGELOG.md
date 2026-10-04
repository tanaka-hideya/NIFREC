# Changelog

## 2.0.0

### Changed (backward-incompatible)
- `nifrec-gaussian-optfreq`: `--imag-vec` was replaced by `--disp-vec {sum,largest}` (default: `sum`).
- Stage 2 uses a deterministic sign convention for the imaginary-mode displacement vectors: the largest-magnitude Cartesian component of each mode vector is made positive before the vector is used (for `sum`, before the vectors are added). Displacement directions may therefore differ from those of v1.x.
- The output .csv file of `nifrec-gaussian-optfreq` contains additional columns (see README.md).
- New dependencies: ASE (reading Gaussian input files, RMSD superposition) and py-cpuinfo (CPU information in the log).
- Terminology: the outputs are described as local-minimum geometries (no imaginary frequencies at the chosen level of theory).

### Added
- `--reverse-disp`: reverse the Stage 2 displacement vector.
- `--skip-stage1`: skip the RCFC re-optimization and go directly from Stage 0 to Stage 2.
- `--infolder-gjf`: read the charge, multiplicity, and coordinates (Cartesian or Z-matrix) from Gaussian input files (.gjf/.com) with ASE instead of using the xTB workflow.
- New output columns (the unit is part of the column name): `fail_stage`; the number and values of imaginary frequencies after each stage (`n_imag_s0`, `n_imag_s1`, `n_imag_s2`, `imag_freqs_s0_per_cm`, `imag_freqs_s1_per_cm`, `imag_freqs_s2_per_cm`); RMSDs and energy differences (`rmsd_in_s0_angstrom`, `rmsd_s0_s1_angstrom`, `dE_s0_s1_kJ_per_mol`, `rmsd_s0_s2_angstrom`, `dE_s0_s2_kJ_per_mol`); and the elapsed time per molecule (`wall_time_seconds`). The columns are defined in `nifrec-gaussian-optfreq --help` and README.md.
- The NIFREC version is written to the logs of all steps; the Gaussian log also contains the host name and CPU.
- `nifrec-gaussian-optfreq --help` shows the complete description of the module (identical to its docstring), including the definitions of the output columns.
- `nifrec-gaussian-parse` carries over all columns of the input .csv file and keeps integer columns as integers.
- Automated tests (pytest) and continuous integration (GitHub Actions).

### Fixed
- An error while processing one molecule in `nifrec-gaussian-optfreq` (e.g., a log file that cannot be parsed) no longer stops the batch; the molecule is recorded as failed (`fail_stage`).
- `nifrec-gaussian-parse`: the final SCF energy is converted from eV back to hartree with `cclib.parser.utils.convertor`, i.e., with the same factor that cclib used for the conversion to eV. A different hartree-eV factor was used before, which shifted the reported energies (and the derived zero-point-corrected energy) by about 4.4e-8 relative (e.g., 2.2e-5 hartree for -500 hartree).
- `--infile-gaussian-recalc`: molecule identifiers are matched as text, so that identifiers such as "001" are handled correctly.
