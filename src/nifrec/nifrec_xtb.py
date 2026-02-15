"""
File: nifrec_xtb.py
Author: Hideya Tanaka at Nara Institute of Science and Technology
Supervised by: Tomoyuki Miyao at Nara Institute of Science and Technology
Description:
    This script performs the conformer optimization and vibrational analysis using xTB, with automatic handling of small imaginary frequencies.
"""

import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem
import json
from joblib import cpu_count, delayed, Parallel
import sys
from pathlib import Path
import shutil
import shlex


def run_xtb_optimization_and_vibration_handle_imagfreq(file_path, charge, outfd, imagfreqoutfd, namespace=None, max_iterations=50, imagfreq_thres=5.0, xcmd='xtb', option_xtb=None):
    if namespace is None:
        namespace = ''

    output_xyz = f'{namespace}.xtbopt.xyz' # Output of optimized geometry by xTB
    output_json = f'{namespace}.xtbout.json'
    output_hessxyz = f'{namespace}.xtbhess.xyz' # Output of distorted structure by xTB in the case where imaginary frequencies are present.
    fname = os.path.basename(file_path)

    def run_xtb_command(file_path, charge, namespace, fname):
        fp = str(Path(file_path).resolve().as_posix())
        xc = str(Path(xcmd).resolve().as_posix())
        args = [xc, fp, '--ohess', '--chrg', str(charge), '--json']
        if namespace != '':
            args += ['--namespace', namespace]
        if option_xtb:
            args += shlex.split(option_xtb)
        try:
            subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f'{fname}, error during xTB run: {e}')
            return False
        return True
    
    # First run
    if not run_xtb_command(file_path, charge, namespace, fname):
        return None
    
    # Handle small imaginary frequencies
    for iteration in range(max_iterations):
        if os.path.exists(output_json):
            is_stable, total_energy = check_vibrational_frequencies_and_energy(output_json, imagfreq_thres)
            if not is_stable and os.path.exists(output_hessxyz):
                if not run_xtb_command(output_hessxyz, charge, namespace, fname):
                    return None
            else:
                break
        else:
            print(f'{output_json} does not exist')
            return None
        
    # Final output
    outfname = fname.replace('rdkit','xTB')
    
    if os.path.exists(output_json):
        is_stable, total_energy = check_vibrational_frequencies_and_energy(output_json, imagfreq_thres)
        if is_stable:
            new_xyz_path = f'{outfd}/{outfname}'
            new_json_path = new_xyz_path.replace('.xyz', '.json')
            new_hessxyz_path = f'{outfd}/{outfname}_hess'
        else:
            new_xyz_path = f'{imagfreqoutfd}/{outfname}'
            new_json_path = new_xyz_path.replace('.xyz', '.json')
            new_hessxyz_path = f'{imagfreqoutfd}/{outfname}_hess'

    if os.path.exists(output_xyz):
        os.rename(output_xyz, new_xyz_path)
    else:
        print(f'{output_xyz} does not exist')
        return None
    if os.path.exists(output_json):
        os.rename(output_json, new_json_path)
    else:
        print(f'{output_json} does not exist')
        return None
    if os.path.exists(output_hessxyz):
        os.rename(output_hessxyz, new_hessxyz_path)

    return new_json_path


def check_vibrational_frequencies_and_energy(json_file, imagfreq_thres=5.0):
    if json_file is None:
        return False, None
    
    with open(json_file, 'r') as file:
        try:
            data = json.load(file)
        except:
            print(f'loading json file failed: {json_file}')
            return False, None

    frequencies = data['vibrational frequencies / rcm']
    total_energy = data['total energy']

    # Check imaginary freq (ignore small freq)
    # The threshold of 5.0 was determined by referring to the official documentation of xTB.
    significant_frequencies = [freq for freq in frequencies if abs(freq) > imagfreq_thres]
    has_imaginary_frequencies = any(freq < 0 for freq in significant_frequencies)

    return not has_imaginary_frequencies, total_energy


def process_rows_for_xtb(outfd, infd, infile, file_path_output_all, file_path_output_min, max_nconf=20, max_repeat=50, imagfreq_thres=5.0, njobs=-1, backend='loky', xcmd='xtb', option_xtb=None):
    print(f'Conformer optimization and vibrational analysis using xTB')
    print('========== Settings ==========')
    print(f'outfolder-xtb: {outfd}')
    print(f'infolder-rdkit: {infd}')
    print(f'infile: {infile}')
    print(f'outfile: {file_path_output_min}')
    print(f'outfile-allresults: {file_path_output_all}')
    print(f'max-nconfs: {max_nconf}')
    print(f'max-repeat: {max_repeat}')
    print(f'imagfreq-thres: {imagfreq_thres}')
    print(f'njobs: {njobs}')
    print(f'backend: {backend}')
    print(f'xcmd: {xcmd}')
    print(f'option-xtb: {option_xtb}')
    print('------------------------------')
    
    file_path_input = f'{infd}/{infile}'
    rdkitxyzfd = f'{infd}/xyz'
    
    df = pd.read_csv(file_path_input, index_col=0)
    print(f'All molecules in the dataset (before RDKit optimization) {len(df)}')
    
    smicol = 'smiles'
    molcol = 'romol'
    original_df = df[df.success_rdkit_confgen]
    print(f'All molecules successfully processed by RDKit conformer generation {len(original_df)}')

    mols_df = original_df.copy()
    
    mols_df[molcol] = mols_df[smicol].apply(Chem.MolFromSmiles)
    mols_df.dropna(subset=molcol, inplace=True) # Drop fail molecules by RDKit
    mols_df['formal_charge'] = mols_df[molcol].apply(Chem.GetFormalCharge)

    print(f'Done. Cheking structures and constraints.\n xTB optimization for {len(mols_df)} mols.')

    # Parallel calculation setting
    njobs = cpu_count() -1 if njobs < 1 else njobs
    batches = np.array_split(mols_df, njobs)

    optoutfd = f'{outfd}/opt'
    imagfreqoutfd = f'{outfd}/imag_freq'
    workingfd = f'{outfd}/working'
    workerfd = f'{outfd}/worker'
    opt_emin_outfd = f'{outfd}/xtbopt_emin_xyz'
    os.makedirs(optoutfd, exist_ok=False)
    os.makedirs(imagfreqoutfd, exist_ok=False)
    os.makedirs(workingfd, exist_ok=False)
    os.makedirs(workerfd, exist_ok=False)
    os.makedirs(opt_emin_outfd, exist_ok=False)
    # xTB calculation
    status_all_df_list = []
    status_min_df_list = []
    results = Parallel(n_jobs=njobs, backend=backend)([delayed(run_xtb_multi_files)(wid, batch, max_nconf, max_repeat, rdkitxyzfd, imagfreq_thres, optoutfd, imagfreqoutfd, workingfd, workerfd, opt_emin_outfd, xcmd, option_xtb) for wid, batch in enumerate(batches)])
    for status_all, status_min in results:
        status_all_df_list.append(status_all)
        status_min_df_list.append(status_min)
        
    pd.concat(status_all_df_list).to_csv(f'{outfd}/{file_path_output_all}')
    status_min_df = pd.concat(status_min_df_list)

    combined_df = pd.concat([original_df, status_min_df], ignore_index=False, axis=1)
    if combined_df[smicol].equals(combined_df['smiles_xTB']):
        combined_df = combined_df.drop(columns='smiles_xTB')
        print('xTB optimization completed.')
    combined_df.to_csv(f'{outfd}/{file_path_output_min}')
    
    print(f'All molecules ({outfd}/{file_path_output_min}) {len(combined_df)}')
    combined_df = combined_df[combined_df['confid'] != 0]
    print(f'All molecules successfully processed by xTB optimization and vibrational analysis (no imaginary frequencies) {len(combined_df)}')
    print('confid = 0 means all conformers have imaginary frequencies.')


def run_xtb_multi_files(workerid, mols, nconfmax, max_repeat, xyz_infd, imagfreq_thres, optoutfd, imagfreqoutfd, workingfd, workerfd, opt_emin_outfd, xcmd, option_xtb):    
    pwd_prev = os.getcwd()
    os.chdir(workingfd) # Output files are stored in wd

    ntotal = len(mols)

    status_all_dict = dict()
    status_min_dict = dict()
    for cumnum, row in enumerate(mols.itertuples(index=True)):

        number = row.Index
        smiles = row.smiles
        nconfs = int(min(row.nconfs_rdkit_confgen, nconfmax))
        charge = row.formal_charge
        emin_cidx = 0 
        emin_val = np.inf
        emin_xyz_path = ''
        
        with open(f'{workerfd}/log_xTB_worker_{workerid}.txt', 'w') as f:
            print(f'Processing {cumnum+1}/{ntotal}, number {number}, smiles {smiles}', file=f)
        
        for confid in range(1, nconfs+1):
            xyz_file_path = f'{xyz_infd}/rdkit_{number}_{confid}.xyz' # 1-indexed
            if not os.path.exists(xyz_file_path):
                raise ValueError(f'xyz file: {xyz_file_path} does not exists. Raise exception.')
            # Run xTB
            json_path = run_xtb_optimization_and_vibration_handle_imagfreq(xyz_file_path, charge, optoutfd, imagfreqoutfd, f'id{workerid}', max_repeat, imagfreq_thres, xcmd, option_xtb)
            
            is_stable, total_energy = check_vibrational_frequencies_and_energy(json_path, imagfreq_thres)

            # Save generation information
            gen_xyz_path  = json_path.replace('.json', '.xyz') if json_path is not None else '' 
            status_all_dict[f'{number}_{confid}'] = {'smiles': smiles,
                                        'molid': number,
                                        'confid': confid if is_stable else 0,
                                        'total_energy_xTB': total_energy,
                                        'filepath': os.path.basename(gen_xyz_path)}

            if is_stable and (total_energy < emin_val):
                emin_cidx = confid
                emin_val = total_energy
                emin_xyz_path = gen_xyz_path

        status_min_dict[number] = {'smiles_xTB': smiles,
                                    'molid': number,
                                    'confid': emin_cidx,
                                    'total_energy_xTB': emin_val,
                                    'filepath': os.path.basename(emin_xyz_path)}
        
        if emin_cidx != 0:
            shutil.copy(f'{optoutfd}/{os.path.basename(emin_xyz_path)}', opt_emin_outfd)
        
    os.chdir(pwd_prev) # set back the previous folder
    return pd.DataFrame.from_dict(status_all_dict, orient='index'), pd.DataFrame.from_dict(status_min_dict, orient='index')


def _parse_cli_args(argv=None):

    parser = argparse.ArgumentParser(
                        prog='nifrec_xtb',
                        description='Conformer optimization and vibrational analysis using xTB, with automatic handling of small imaginary frequencies.',
    )
    parser.add_argument('--outfolder-xtb',
                     help="Output folder to write xTB results (XYZ files, logs). Accepts absolute or relative paths; '~' is expanded. The folder is created.",
                     type=str,
                     required=True,
                     )
    parser.add_argument('--infolder-rdkit',
                     help="Input folder for loading RDKit results (XYZ files). Accepts absolute or relative paths; '~' is expanded.",
                     type=str,
                     required=True,
                     )
    parser.add_argument('--infile',
                     help="Name of the input CSV file (RDKit summary) located under --infolder-rdkit. (default: rdkit_stats.csv)",
                     type=str,
                     default='rdkit_stats.csv',
                     )
    parser.add_argument('--outfile',
                     help='Name of the output CSV file to write summary (saved under --outfolder-xtb). (default: xTB_stats_Emin.csv)',
                     type=str,
                     default='xTB_stats_Emin.csv',
                     )
    parser.add_argument('--outfile-allresults',
                     help='Name of the output CSV file to write all results (saved under --outfolder-xtb). (default: xTB_stats_all.csv)',
                     type=str,
                     default='xTB_stats_all.csv',
                     )
    parser.add_argument('--max-nconfs',
                     help='Maximum number of conformers to be optimized with xTB per molecule. (default: 20)',
                     type=int,
                     default=20,
                     )
    parser.add_argument('--max-repeat',
                help=("Upper bound on iterations to resolve residual imaginary frequencies after xTB optimization "
                    "by following the distorted structure."
                    " Larger values allow more retries; set small to limit runtime (positive integer). (default: 50)"),
                     type=int,
                     default=50,
                     )
    parser.add_argument('--imagfreq-thres',
               help=("Threshold in cm^-1 to ignore near-zero vibrational modes when detecting imaginary frequencies. "
                   "Modes with |freq| <= threshold are treated as numerical noise. "
                   "Negative modes with |freq| > threshold are considered imaginary. "
                   "Applies both to the iterative re-optimization and final classification. "
                   "(default: 5.0)"),
                     type=float,
                     default=5.0,
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
    parser.add_argument('--xcmd',
                     help="Command to the xTB executable. 'xcmd file.xyz' (default: xtb)",
                     type=str,
                     default='xtb',
                     )
    parser.add_argument('--option-xtb',
                     help=("Additional command-line arguments passed to xTB. "
                           "Example: --option-xtb '--alpb water' (default: None)"),
                     type=str,
                     default=None,
                     )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_cli_args(argv)
    outfd = str(Path(args.outfolder_xtb).expanduser().resolve())
    infd = str(Path(args.infolder_rdkit).expanduser().resolve())
    os.makedirs(outfd)
    sys.stdout = open(f'{outfd}/log_xtb.txt', 'w')
    process_rows_for_xtb(outfd, infd, args.infile, args.outfile_allresults, args.outfile, args.max_nconfs, args.max_repeat, args.imagfreq_thres, args.njobs, args.backend, args.xcmd, args.option_xtb)
    print('Finish')
    sys.stdout.close()
    
    
if __name__ == '__main__':
    main()
