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

def gaussian_analyze(infd, number, smiles, glogfd, filepath):  
    results = {'is_success': False,
                'zero_corr': None,
                'E_corr': None,
                'H_corr': None,
                'G_corr': None,
                'Ezero': None,
                'Ethermal': None,
                'H': None,
                'G': None,
                'HOMO': None,
                'LUMO': None}
    
    glogpath = f'{glogfd}/{filepath}'
    if not os.path.exists(glogpath):
        print(f'{glogpath} does not exist')
        return results

    # Flag to check success of opt
    log_line = 0
    # Flag to manage pairing between Alpha occ. and Alpha virt.
    expect_alpha_virt = False
    current_alpha_occ = []

    try:
        # cclib
        data = cclib.io.ccread(glogpath)
        # Check imaginary freq
        # vibfreqs: vibrational frequencies, 1/cm, array of rank 1 (cclib parsed data (version 1.8.1))
        vibfreqs = data.vibfreqs
        with open(f'{infd}/freq_gaussian_parse.txt', 'a') as f:
            print(f'number {number}, smiles {smiles}, freq {vibfreqs}', file=f)
        is_stable = any(freq < 0 for freq in vibfreqs)
        if is_stable:
            print(f'imaginary freq: number {number}, smiles {smiles}')
            return results
        
        with open(glogpath, 'r') as f:
            for line in f:
                # Check success of opt
                if gchk in line and log_line == 0:
                    log_line = 1
                if gjf_key in line and log_line == 1:
                    log_line = 2
                if 'Optimized Parameters' in line and log_line == 2:
                    log_line = 3
                if 'Normal termination of Gaussian' in line and log_line == 3:
                    results['is_success'] = True
                    
                # Parse
                for key, pattern in compiled_patterns.items():
                    match = pattern.search(line)
                    if match:
                        if key == 'Alpha_occ':
                            energy_values = match.group(1).strip().split()
                            energy_values = [float(val) for val in energy_values]
                            current_alpha_occ = energy_values
                            # Set flag to expect the next Alpha virt. line
                            expect_alpha_virt = True
                        elif key == 'Alpha_virt':
                            if expect_alpha_virt:
                                energy_values = match.group(1).strip().split()
                                energy_values = [float(val) for val in energy_values]
                                if current_alpha_occ:
                                    results['HOMO'] = current_alpha_occ[-1]  # Last occupied is HOMO
                                if energy_values:
                                    results['LUMO'] = energy_values[0]   # First virtual is LUMO
                                # Reset the flag
                                expect_alpha_virt = False
                        else:
                            results[key] = float(match.group(1))
                            
        # Error check by using cclib
        if results['H'] != data.enthalpy or results['G'] != data.freeenergy:
            print(f'parse error based on cclib, number {number}, smiles {smiles}')
            results['is_success'] = False
            return results
                    
    except Exception as e:
        print(f'parse error: {e}, number {number}, smiles {smiles}')
        results['is_success'] = False
        return results
    
    return results

def process_rows_for_gparse(infd, glogfd, infile, outfile):
    print(f'Gaussian parse (opt freq)')
    print('========== Settings ==========')
    print(f'infolder-gaussian: {infd}')
    print(f'infolder-gaussian-log: {glogfd}')
    print(f'infile: {infile}')
    print(f'outfile: {outfile}')
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
            
        results_dict = gaussian_analyze(infd, number, smiles, glogfd, filepath)

        if results_dict['is_success']:
            status_dict[number] = {'smiles': smiles,
                                    'confid': confid,
                                    'charge': charge,
                                    'multiplicity': multiplicity,
                                    'total_energy_xTB': energy_xtb,
                                    'Zero_point_correction': results_dict['zero_corr'],
                                    'Thermal_correction_to_Energy': results_dict['E_corr'],
                                    'Thermal_correction_to_Enthalpy': results_dict['H_corr'],
                                    'Thermal_correction_to_Gibbs_Free_Energy': results_dict['G_corr'],
                                    'Sum_of_electronic_and_zero_point_Energies': results_dict['Ezero'],
                                    'Sum_of_electronic_and_thermal_Energies': results_dict['Ethermal'],
                                    'Sum_of_electronic_and_thermal_Enthalpies': results_dict['H'],
                                    'Sum_of_electronic_and_thermal_Free_Energies': results_dict['G'],
                                    'HOMO': results_dict['HOMO'],
                                    'LUMO': results_dict['LUMO'],
                                    'filepath': filepath,
                                    'success_stage': success_stage,
                                    'success_disploop': success_disploop}
        else:
            status_dict[number] = {'smiles': smiles,
                                    'confid': 0,
                                    'charge': charge,
                                    'multiplicity': multiplicity,
                                    'total_energy_xTB': energy_xtb,
                                    'Zero_point_correction': results_dict['zero_corr'],
                                    'Thermal_correction_to_Energy': results_dict['E_corr'],
                                    'Thermal_correction_to_Enthalpy': results_dict['H_corr'],
                                    'Thermal_correction_to_Gibbs_Free_Energy': results_dict['G_corr'],
                                    'Sum_of_electronic_and_zero_point_Energies': results_dict['Ezero'],
                                    'Sum_of_electronic_and_thermal_Energies': results_dict['Ethermal'],
                                    'Sum_of_electronic_and_thermal_Enthalpies': results_dict['H'],
                                    'Sum_of_electronic_and_thermal_Free_Energies': results_dict['G'],
                                    'HOMO': results_dict['HOMO'],
                                    'LUMO': results_dict['LUMO'],
                                    'filepath': filepath,
                                    'success_stage': success_stage,
                                    'success_disploop': success_disploop}
            
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
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_cli_args(argv)
    infd = str(Path(args.infolder_gaussian).expanduser().resolve())
    infd_log = str(Path(args.infolder_gaussian_log).expanduser().resolve())
    sys.stdout = open(f'{infd}/log_gaussian_parse.txt', 'w')
    process_rows_for_gparse(infd, infd_log, args.infile, args.outfile)
    print('Finish')
    sys.stdout.close()
    
    
if __name__ == '__main__':
    main()

    