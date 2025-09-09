"""
File: nifrec_gaussian_parse.py
Author: Hideya Tanaka at Nara Institute of Science and Technology
Supervised by: Tomoyuki Miyao at Nara Institute of Science and Technology
Description:
    This script performs parse of Gaussian results.
"""

import argparse
import os
import pandas as pd
import sys
import cclib
from pathlib import Path

EV_PER_HARTREE = 27.211386245988  


def gaussian_analyze(infd, number, smiles, glogfd, filepath, no_homo_lumo):  
    results = {'is_success': False,
                'functional': None,
                'basis_set': None,
                'E_scf_final': None,
                'zpve': None,
                'Ezero': None,
                'H': None,
                'G': None,
                'T': None,
                'HOMO': None,
                'LUMO': None}
    
    glogpath = f'{glogfd}/{filepath}'
    if not os.path.exists(glogpath):
        print(f'{glogpath} does not exist')
        return results

    try:
        # cclib
        data = cclib.io.ccread(glogpath)
        results['functional'] = data.metadata['functional']
        results['basis_set'] = data.metadata['basis_set']
        
        # Check imaginary freq
        # vibfreqs: vibrational frequencies, 1/cm, array of rank 1 (cclib parsed data (version 1.8.1))
        vibfreqs = data.vibfreqs
        with open(f'{infd}/freq_gaussian_parse.txt', 'a') as f:
            print(f'number {number}, smiles {smiles}, freq {vibfreqs}', file=f)
        is_stable = any(freq < 0 for freq in vibfreqs)
        if is_stable:
            print(f'imaginary freq: number {number}, smiles {smiles}')
            return results

        # Check success of opt
        if data.metadata['success'] and data.optdone:
            results['is_success'] = True

        results['E_scf_final'] = data.scfenergies[-1] / EV_PER_HARTREE
        results['zpve'] = data.zpve
        results['Ezero'] = results['E_scf_final'] + results['zpve']
        results['H'] = data.enthalpy
        results['G'] = data.freeenergy
        results['T'] = data.temperature
        
        if no_homo_lumo:
            return results
        
        # restricted method version
        homoidx = data.homos[0]
        lumoidx = homoidx + 1
        homo_energy = data.moenergies[0][homoidx]
        lumo_energy = data.moenergies[0][lumoidx]
        results['HOMO'] = homo_energy
        results['LUMO'] = lumo_energy

    except Exception as e:
        print(f'parse error: {e}, number {number}, smiles {smiles}')
        results['is_success'] = False
        return results
    
    return results


def process_rows_for_gparse(infd, glogfd, infile, outfile, no_homo_lumo):
    print(f'Gaussian parse (opt freq)')
    print('========== Settings ==========')
    print(f'infolder-gaussian: {infd}')
    print(f'infolder-gaussian-log: {glogfd}')
    print(f'infile: {infile}')
    print(f'outfile: {outfile}')
    print(f'no-homo-lumo: {no_homo_lumo}')
    print('------------------------------')
    
    # Safety check: prevent accidental overwrite when input and output filenames are identical.
    # Using the same name would overwrite the original Gaussian summary.
    if infile == outfile:
        raise ValueError('The input file name (--infile) and output file name (--outfile) are identical. Use a different --outfile to avoid overwriting the source CSV.')

    file_path_input = f'{infd}/{infile}'
    file_path_output = f'{infd}/{outfile}'
            
    df = pd.read_csv(file_path_input, index_col=0)
    print(f'All molecules in the dataset {len(df)}')
    df = df[df['confid'] != 0]
    print(f"All molecules successfully processed by Gaussian 'opt freq' calculations (no imaginary frequencies) {len(df)}")

    ntotal = len(df)
    status_dict = dict()
    for cumnum, row in enumerate(df.itertuples(index=True)):
        
        number = row.Index
        smiles = row.smiles
        confid = row.confid
        charge = row.charge
        multiplicity = row.multiplicity
        energy_xtb = row.total_energy_xTB
        filepath  = row.filepath
        success_stage = row.success_stage
        success_disploop = row.success_disploop
        
        with open(f'{infd}/log_gaussian_worker_parse.txt', 'w') as f:
            print(f'Processing {cumnum+1}/{ntotal}, number {number}, smiles {smiles}', file=f)
            
        results_dict = gaussian_analyze(infd, number, smiles, glogfd, filepath, no_homo_lumo)

        if results_dict['is_success']:
            status_dict[number] = {'smiles': smiles,
                                    'confid': confid,
                                    'charge': charge,
                                    'multiplicity': multiplicity,
                                    'total_energy_xTB': energy_xtb,
                                    'filepath': filepath,
                                    'success_stage': success_stage,
                                    'success_disploop': success_disploop,
                                    'functional': results_dict['functional'],
                                    'basis_set': results_dict['basis_set'],
                                    'Final_SCF_Energy_hartree': results_dict['E_scf_final'],
                                    'Zero_point_correction_hartree': results_dict['zpve'],
                                    'Sum_of_electronic_and_zero_point_Energies_hartree': results_dict['Ezero'],
                                    'Sum_of_electronic_and_thermal_Enthalpies_hartree': results_dict['H'],
                                    'Sum_of_electronic_and_thermal_Free_Energies_hartree': results_dict['G'],
                                    'Temperature_K': results_dict['T'],
                                    'HOMO_eV': results_dict['HOMO'],
                                    'LUMO_eV': results_dict['LUMO']}
        else:
            status_dict[number] = {'smiles': smiles,
                                    'confid': 0,
                                    'charge': charge,
                                    'multiplicity': multiplicity,
                                    'total_energy_xTB': energy_xtb,
                                    'filepath': filepath,
                                    'success_stage': success_stage,
                                    'success_disploop': success_disploop,
                                    'functional': results_dict['functional'],
                                    'basis_set': results_dict['basis_set'],
                                    'Final_SCF_Energy_hartree': results_dict['E_scf_final'],
                                    'Zero_point_correction_hartree': results_dict['zpve'],
                                    'Sum_of_electronic_and_zero_point_Energies_hartree': results_dict['Ezero'],
                                    'Sum_of_electronic_and_thermal_Enthalpies_hartree': results_dict['H'],
                                    'Sum_of_electronic_and_thermal_Free_Energies_hartree': results_dict['G'],
                                    'Temperature_K': results_dict['T'],
                                    'HOMO_eV': results_dict['HOMO'],
                                    'LUMO_eV': results_dict['LUMO']}
            
        status_df = pd.DataFrame.from_dict(status_dict, orient='index')
        status_df.to_csv(file_path_output)
        
    status_df = status_df[status_df['confid'] != 0]
    print(f"All molecules successfully processed by Gaussian 'opt freq' calculations (no imaginary frequencies) {len(status_df)}")
    print("confid = 0 means the failed Gaussian 'opt freq' calculations or that imaginary frequencies could not be removed.")
    print("confid != 0: check 'optdone'")


def _parse_cli_args(argv=None):

    parser = argparse.ArgumentParser(
                        prog='nifrec_gaussian_parse',
                        description='This script parses Gaussian output files (opt freq).',
    )
    parser.add_argument('--infolder-gaussian',
               help=("Input folder for loading Gaussian results (CSV file). (Also the output folder to which parsed results are written.) "
                   "Accepts absolute or relative paths; '~' is expanded. "
                   "Advanced use case: You may manually re-run Gaussian for molecules that failed in nifrec_gaussian_optfreq before executing this parser. "
                   "If you adopt this workflow you MUST (1) edit the CSV given via --infile so each successfully re-computed molecule has a non-zero 'confid' and a correct 'filepath' pointing to the re-calculated structure; and (2) gather the corresponding successful Gaussian '.log' files into the directory passed via --infolder-gaussian-log. "
                   "For clarity and reproducibility, keep the original directory produced by nifrec_gaussian_optfreq unchanged; instead make a copy of that directory, perform the manual re-calculations there, and run this parsing step within the copied directory."),
                     type=str,
                     required=True,
                     )
    parser.add_argument('--infolder-gaussian-log',
                help=("Input folder for loading Gaussian .log files. Accepts absolute or relative paths; '~' is expanded."
                      "typically, it is located under --infolder-gaussian."),
                     type=str,
                     required=True,
                     )
    parser.add_argument('--infile',
                     help="Name of the input CSV file (Gaussian summary) located under --infolder-gaussian.",
                     type=str,
                     required=True,
                     )
    parser.add_argument('--outfile',
                     help='Name of the output CSV file to write parsed results (saved under --infolder-gaussian).',
                     type=str,
                     required=True,
                     )
    parser.add_argument('--no-homo-lumo',
               help=("Default: extract HOMO/LUMO energies from MO data (RESTRICTED METHODS ONLY). "
                   "Values are reported in eV. Unrestricted (UHF) logs are not supported."),
                     action='store_true',
                     )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_cli_args(argv)
    infd = str(Path(args.infolder_gaussian).expanduser().resolve())
    infd_log = str(Path(args.infolder_gaussian_log).expanduser().resolve())
    sys.stdout = open(f'{infd}/log_gaussian_parse.txt', 'w')
    process_rows_for_gparse(infd, infd_log, args.infile, args.outfile, args.no_homo_lumo)
    print('Finish')
    sys.stdout.close()
    
    
if __name__ == '__main__':
    main()

    