"""
File: 1_1_nifrec_rdkit.py
Author: Hideya Tanaka at Nara Institute of Science and Technology
Supervised by: Tomoyuki Miyao
Description:
    This script performs the conformer generation using RDKit.
"""
import argparse
import os
import numpy as np
import pandas as pd
from joblib import delayed, Parallel, cpu_count
from morfeus.conformer import ConformerEnsemble
from utility_opt import MakeFolder, SmilesToCanSmiles, WriteMolsSDF
import sys

def process_rows_for_rdkit(fd, file_path_input, file_path_output, n_confs, thres, njobs=-1,  backend='loky',smicol='smiles', idxcol=0):
    outfd       = MakeFolder(f'{fd}/rdkit', allow_override=True) 
    original_df = pd.read_csv(file_path_input, index_col=idxcol)
    print(f'Loaded mols: {len(original_df)}')
    df = original_df.copy()
    
    # Remove duplicates
    df[smicol] = original_df[smicol].apply(SmilesToCanSmiles)
    df.drop_duplicates(subset=smicol, inplace=True)
    print(f'Loaded mols after removing duplicates: {len(df)}')
    
    # Parallel calculation setting     
    njobs = cpu_count() -1 if njobs < 1 else njobs
    batch_df = np.array_split(df, njobs)

    random_seed = 42
    
    # Both sdf and xyz are necessary to reproduce molecule after coordinate optimization.
    outfd_xyz = MakeFolder(f'{outfd}/xyz', allow_override=True)
    outfd_sdf = MakeFolder(f'{outfd}/sdf', allow_override=True)
    status_df_list = Parallel(n_jobs=njobs, backend=backend)([delayed(generate_conf_rdkit)(wid, batch[smicol], n_confs, random_seed, thres, outfd_xyz, outfd_sdf, outfd) for wid, batch in enumerate(batch_df)])
    status_df = pd.concat(status_df_list)
    combined_df = pd.concat([original_df, status_df], ignore_index=False, axis=1)
    if combined_df[smicol].equals(combined_df['smiles_input_rdkit_confgen']):
        combined_df = combined_df.drop(columns='smiles_input_rdkit_confgen')
        print('success: combine')
    combined_df.to_csv(f'{fd}/{file_path_output}')

def generate_conf_rdkit(workerid, pds_smiles, nconfs, rseed, rmsdthres, outfd_xyz, outfd_sdf, outfd):
    ntotal = len(pds_smiles)
    status_dict = dict()
    for idx, (number,smiles) in enumerate(pds_smiles.items()):
        with open(f'{outfd}/log_rdkit_worker_{workerid}.txt', 'w') as f:
            print(f'Processing {idx+1}/{ntotal}, number {number}, smiles {smiles}', file=f)

        try:
            ce = ConformerEnsemble.from_rdkit(smiles, n_confs=nconfs, optimize='MMFF94', random_seed=rseed, rmsd_thres=0)
            ce.prune_rmsd(thres=rmsdthres)
            ce.sort()
            n_conformer = len(ce)
            ce.write_xyz(f'{outfd_xyz}/rdkit_{number}.xyz', ids=None, unit='hartree', relative=False, separate=True, n_decimals=10)
            
            # Save mol to connectivity information.
            rdmol = ce.mol
            rdmol.SetProp('_Name', str(number))
            WriteMolsSDF([rdmol], f'{outfd_sdf}/rdkit2d_{number}.sdf')            
            status_dict[number] = {'smiles_input_rdkit_confgen': smiles, 'success_rdkit_confgen': True, 'nconfs_rdkit_confgen': n_conformer}
        except:
            print(f'Fail RDkit conformation generation process: {idx}: {smiles}')
            status_dict[number] = {'smiles_input_rdkit_confgen': smiles, 'success_rdkit_confgen': False, 'nconfs_rdkit_confgen': None}
    return pd.DataFrame.from_dict(status_dict, orient='index')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                        prog='nifrec_rdkit',
                        description= 'Conformer generator using RDKit'
    )
    parser.add_argument('--infile',
                     help='Input csv filename',
                     type=str
                     )
    parser.add_argument('--smicol',
                     help='Smiles column for structure generation',
                     type=str,
                     default='smiles'
                     )
    parser.add_argument('--idxcol',
                     help='Index column idx (0 starts)',
                     type=int,
                     default=0
                     )
    parser.add_argument('--outfolder',
                     help='Output folder where conformers are stored',
                     type=str,
                     default='rdkit-confgen'
                     )
    parser.add_argument('--outfile',
                     help='Output csv filename without extension (i.e. csv)',
                     type=str,
                     default='rdkit_stats'
                     )
    parser.add_argument('--nconfs',
                     help='Number of maximum conformers to be generated',
                     type=int,
                     default=20,
                     )
    parser.add_argument('--rmsd-thres',
                     help='RMSD threshold between generated conformers',
                     type=float,
                     default=1,
                     )
    parser.add_argument('--njobs',
                     help='Number of parallel computation threads (-1 means all cores)',
                     type=int,
                     default=-1,
                     )
    parser.add_argument('--backend',
                     help='Parallel computation backend (loky or multithprocessing)',
                     type=str,
                     default='loky',
                     )
    p = parser.parse_args(sys.argv[1:])

    fd = p.outfolder
    fd = MakeFolder(fd, allow_override=True)
    sys.stdout = open(f'{fd}/log_rdkit.txt', 'w')
    process_rows_for_rdkit(fd, p.infile, p.outfile, p.nconfs, p.rmsd_thres, njobs=p.njobs, backend=p.backend,smicol=p.smicol)
    print('Finish')
    sys.stdout.close()