"""
File: nifrec_rdkit.py
Author: Hideya Tanaka at Nara Institute of Science and Technology
Supervised by: Tomoyuki Miyao at Nara Institute of Science and Technology
Description:
    This script performs the conformer generation using RDKit.
"""

import argparse
import os
import numpy as np
import pandas as pd
from joblib import delayed, Parallel, cpu_count
from morfeus.conformer import ConformerEnsemble
from rdkit import Chem
import sys
from pathlib import Path


def smiles_to_canonical(smi):
    if (smi is None) or (smi == ''):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    else:
        return Chem.MolToSmiles(mol) # canonical smiles are generated


def write_mols_to_sdf(mol_list, sdf_name):  
    w = Chem.SDWriter(sdf_name)
    for mol in mol_list:
        w.write(mol)
    w.close()


def process_rows_for_rdkit(outfd, file_path_input, file_path_output, n_confs, thres, njobs=-1, backend='loky', smicol='smiles', idxcol=0, random_seed=42, force_field='MMFF94s'):
    print(f'Conformer generation using RDKit')
    print('========== Settings ==========')
    print(f'outfolder-rdkit: {outfd}')
    print(f'infile: {file_path_input}')
    print(f'smicol: {smicol}')
    print(f'idxcol: {idxcol}')
    print(f'outfile: {file_path_output}')
    print(f'nconfs: {n_confs}')
    print(f'rmsd-thres: {thres}')
    print(f'njobs: {njobs}')
    print(f'backend: {backend}')
    print(f'random-seed: {random_seed}')
    print(f'force-field: {force_field}')
    print('------------------------------')
    
    original_df = pd.read_csv(file_path_input, index_col=idxcol)
    print(f'Loaded mols: {len(original_df)}')
    original_df[smicol] = original_df[smicol].apply(smiles_to_canonical)
    df = original_df.copy()
    
    # Check duplicates
    n_before = len(df)
    n_after = len(df.drop_duplicates(subset=smicol))
    if n_after != n_before:
        raise ValueError(f"Duplicate SMILES detected in column '{smicol}'")
    
    # Parallel calculation setting     
    njobs = cpu_count() -1 if njobs < 1 else njobs
    batch_df = np.array_split(df, njobs)
    
    # Both sdf and xyz are necessary to reproduce molecule after coordinate optimization.
    outfd_xyz = f'{outfd}/xyz'
    outfd_sdf = f'{outfd}/sdf'
    outfd_worker = f'{outfd}/worker'
    os.makedirs(outfd_xyz, exist_ok=False)
    os.makedirs(outfd_sdf, exist_ok=False)
    os.makedirs(outfd_worker, exist_ok=False)
    status_df_list = Parallel(n_jobs=njobs, backend=backend)([delayed(generate_conf_rdkit)(wid, batch[smicol], n_confs, random_seed, thres, outfd_xyz, outfd_sdf, outfd_worker, force_field) for wid, batch in enumerate(batch_df)])
    status_df = pd.concat(status_df_list)
    combined_df = pd.concat([original_df, status_df], ignore_index=False, axis=1)
    if combined_df[smicol].equals(combined_df['smiles_input_rdkit_confgen']):
        combined_df = combined_df.drop(columns='smiles_input_rdkit_confgen')
        print('RDKit conformer generation completed.')
    
    # Harmonize SMILES column name to 'smiles' for output, unless it already exists.
    # This avoids duplicate or conflicting columns in the final CSV.
    if smicol != 'smiles':
        if 'smiles' in combined_df.columns:
            print(f"Column 'smiles' already exists; skip renaming '{smicol}' to 'smiles'.")
        else:
            combined_df = combined_df.rename(columns={smicol: 'smiles'})
            print(f"Renamed column '{smicol}' to 'smiles' for output consistency.")
    combined_df.to_csv(f'{outfd}/{file_path_output}')


def generate_conf_rdkit(workerid, pds_smiles, nconfs, rseed, rmsdthres, outfd_xyz, outfd_sdf, outfd_worker, force_field):
    ntotal = len(pds_smiles)
    status_dict = dict()
    for idx, (number,smiles) in enumerate(pds_smiles.items()):
        with open(f'{outfd_worker}/log_rdkit_worker_{workerid}.txt', 'w') as f:
            print(f'Processing {idx+1}/{ntotal}, number {number}, smiles {smiles}', file=f)

        try:
            ce = ConformerEnsemble.from_rdkit(smiles, n_confs=nconfs, optimize=force_field, random_seed=rseed, rmsd_thres=0)
            ce.prune_rmsd(thres=rmsdthres)
            ce.sort()
            n_conformer = len(ce)
            ce.write_xyz(f'{outfd_xyz}/rdkit_{number}.xyz', ids=None, unit='hartree', relative=False, separate=True, n_decimals=10)
            
            # Save mol to connectivity information.
            rdmol = ce.mol
            rdmol.SetProp('_Name', str(number))
            write_mols_to_sdf([rdmol], f'{outfd_sdf}/rdkit2d_{number}.sdf')            
            status_dict[number] = {'smiles_input_rdkit_confgen': smiles, 'success_rdkit_confgen': True, 'nconfs_rdkit_confgen': n_conformer}
        except:
            print(f'Fail RDKit conformation generation process: {idx}: {smiles}')
            status_dict[number] = {'smiles_input_rdkit_confgen': smiles, 'success_rdkit_confgen': False, 'nconfs_rdkit_confgen': None}
    return pd.DataFrame.from_dict(status_dict, orient='index')


def _parse_cli_args(argv=None):

    parser = argparse.ArgumentParser(
                        prog='nifrec_rdkit',
                        description='Conformer generator using RDKit',
    )
    parser.add_argument('--outfolder-rdkit',
                     help="Output folder to write RDKit results (XYZ and SDF files, logs). Accepts absolute or relative paths; '~' is expanded. The folder is created.",
                     type=str,
                     required=True,
                     )
    parser.add_argument('--infile',
                     help="Path to the input CSV file. Accepts absolute or relative paths; '~' is expanded. Must contain a SMILES column specified by --smicol.",
                     type=str,
                     required=True,
                     )
    parser.add_argument('--smicol',
                     help='Name of the column in the input CSV that contains SMILES strings (used for structure generation). (default: smiles)',
                     type=str,
                     default='smiles',
                     )
    parser.add_argument('--idxcol',
               help=('Zero-based index of the column in the input CSV to use as the unique molecule identifier (DataFrame index). '
                   'Identifiers must be unique per molecule and are used consistently across all outputs: the output CSV (--outfile) '
                   'and the filenames of 3D structure files (XYZ/SDF). Non-unique values may cause file overwrites and inconsistent results. (default: 0)'),
                     type=int,
                     default=0,
                     )
    parser.add_argument('--outfile',
               help=("Name of the output CSV file to write summary (saved under --outfolder-rdkit). "
                   "Note: If the SMILES column specified by --smicol is not named 'smiles', it will be renamed to 'smiles' in the output CSV "
                   "(unless a 'smiles' column already exists). Canonical smiles are stored in the 'smiles' column. (default: rdkit_stats.csv)"),
                     type=str,
                     default='rdkit_stats.csv',
                     )
    parser.add_argument('--nconfs',
                     help='Maximum number of RDKit conformers to generate per molecule. (default: 20)',
                     type=int,
                     default=20,
                     )
    parser.add_argument('--rmsd-thres',
                     help='RMSD pruning threshold between generated conformers. (default: 1)',
                     type=float,
                     default=1,
                     )
    parser.add_argument('--njobs',
                     help='Number of parallel workers. If <= 0, uses (CPU cores - 1). (default: -1)',
                     type=int,
                     default=-1,
                     )
    parser.add_argument('--backend',
                     help="Parallel backend for joblib. One of: 'loky' (default), 'multiprocessing', 'threading'. (default: loky)",
                     type=str,
                     default='loky',
                     )
    parser.add_argument('--random-seed',
                     help="Random seed for reproducibility. (default: 42)",
                     type=int,
                     default=42,
                     )
    parser.add_argument('--force-field',
                     help="Force field to use for RDKit conformer generation (MMFF94s, MMFF94, or UFF). (default: MMFF94s)",
                     type=str,
                     default='MMFF94s',
                     )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_cli_args(argv)
    outfd = str(Path(args.outfolder_rdkit).expanduser().resolve())
    infile = str(Path(args.infile).expanduser().resolve())
    os.makedirs(outfd)
    sys.stdout = open(f'{outfd}/log_rdkit.txt', 'w')
    process_rows_for_rdkit(outfd, infile, args.outfile, args.nconfs, args.rmsd_thres, njobs=args.njobs, backend=args.backend, smicol=args.smicol, idxcol=args.idxcol, random_seed=args.random_seed, force_field=args.force_field)
    print('Finish')
    sys.stdout.close()
    
    
if __name__ == '__main__':
    main()