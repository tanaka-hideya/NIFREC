"""
File: 2_2_nifrec_gaussian_parse.py
Author: Hideya Tanaka
Supervised by: Tomoyuki Miyao
Description:
    This script performs parse of Gaussian results.
"""

import os
import pandas as pd
import sys
import re
import cclib
from utility_opt import MakeFolder

def gaussian_analyze(number, smiles, glogfd, gchkfd, ggjffd, filepath, keyword, logoutfd, chkoutfd, gjfoutfd, gjf_key):
    
    patterns = {'zero_corr': r'Zero-point correction=\s+([-+]?\d*\.\d+|\d+)',
                'E_corr': r'Thermal correction to Energy=\s+([-+]?\d*\.\d+|\d+)',
                'H_corr': r'Thermal correction to Enthalpy=\s+([-+]?\d*\.\d+|\d+)',
                'G_corr': r'Thermal correction to Gibbs Free Energy=\s+([-+]?\d*\.\d+|\d+)',
                'Ezero': r'Sum of electronic and zero-point Energies=\s+([-+]?\d*\.\d+|\d+)',
                'Ethermal': r'Sum of electronic and thermal Energies=\s+([-+]?\d*\.\d+|\d+)',
                'H': r'Sum of electronic and thermal Enthalpies=\s+([-+]?\d*\.\d+|\d+)',
                'G': r'Sum of electronic and thermal Free Energies=\s+([-+]?\d*\.\d+|\d+)',
                'Alpha_occ': r'Alpha\s+occ\.\s+eigenvalues\s+--\s+([-+]?\d*\.\d+(?:\s+[-+]?\d*\.\d+)*)',
                'Alpha_virt': r'Alpha\s+virt\.\s+eigenvalues\s+--\s+([-+]?\d*\.\d+(?:\s+[-+]?\d*\.\d+)*)'}

    compiled_patterns = {key: re.compile(pattern) for key, pattern in patterns.items()}
  
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
    
    gchk = filepath.replace('.log', '.chk')
    gchkpath = f'{gchkfd}/{gchk}'
    if not os.path.exists(gchkpath):
        print(f'{gchkpath} does not exist')
        return results
    
    ggjf = filepath.replace('.log', '.gjf')
    ggjfpath = f'{ggjffd}/{ggjf}'
    if not os.path.exists(ggjfpath):
        print(f'{ggjfpath} does not exist')
        return results
        
    try:
        # cclib
        data = cclib.io.ccread(glogpath)
        # Check imaginary freq
        # vibfreqs: vibrational frequencies, 1/cm, array of rank 1 (cclib parsed data (version 1.8.1))
        vibfreqs = data.vibfreqs
        with open(f'{fd}/freq_gaussian_{keyword}_analyze.txt', 'a') as f:
            print(f'number {number}, smiles {smiles}, freq {vibfreqs}', file=f)
        is_stable = any(freq < 0 for freq in vibfreqs)
        if is_stable:
            print(f'imaginary freq: number {number}, smiles {smiles}')
            os.rename(f'{glogpath}', f'{logoutfd}/{filepath}')
            os.rename(f'{gchkpath}', f'{chkoutfd}/{gchk}')
            os.rename(f'{ggjfpath}', f'{gjfoutfd}/{ggjf}')
            return results
    except Exception as e:
        print(f'vibfreqs parse error: {e}, number {number}, smiles {smiles}')
        os.rename(f'{glogpath}', f'{logoutfd}/{filepath}')
        os.rename(f'{gchkpath}', f'{chkoutfd}/{gchk}')
        os.rename(f'{ggjfpath}', f'{gjfoutfd}/{ggjf}')
        return results
    
    # Flag to check success of opt
    log_line = 0
    # Flag to manage pairing between Alpha occ. and Alpha virt.
    expect_alpha_virt = False
    current_alpha_occ = []

    try:
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
            os.rename(f'{glogpath}', f'{logoutfd}/{filepath}')
            os.rename(f'{gchkpath}', f'{chkoutfd}/{gchk}')
            os.rename(f'{ggjfpath}', f'{gjfoutfd}/{ggjf}')
            return results
                    
    except Exception as e:
        print(f'parse error: {e}, number {number}, smiles {smiles}')
        results['is_success'] = False
        os.rename(f'{glogpath}', f'{logoutfd}/{filepath}')
        os.rename(f'{gchkpath}', f'{chkoutfd}/{gchk}')
        os.rename(f'{ggjfpath}', f'{gjfoutfd}/{ggjf}')
        return results
    
    return results

def process_rows_for_gparse(fd, file_path_input, glogfd, gchkfd, ggjffd, file_path_output, keyword, gjf_key):
    print(f'Gaussian parse (opt freq)')
    print('========== Settings ==========')
    print(f'keyword: {keyword}')
    print(f'gjf_route_section_input_for_check: {gjf_key}')
    print(f'file_path_input: {file_path_input}')
    print(f'glogfd: {glogfd}')
    print(f'gchkfd: {gchkfd}')
    print(f'ggjffd: {ggjffd}')
    print(f'file_path_output: {file_path_output}')
    print('------------------------------')
    
    outfd = MakeFolder(f'{fd}/gaussian_{keyword}_analyze_failure', allow_override=True)
    logoutfd = MakeFolder(f'{outfd}/gaussian_log_{keyword}_analyze_failure', allow_override=True)
    chkoutfd = MakeFolder(f'{outfd}/gaussian_chk_{keyword}_analyze_failure', allow_override=True)
    gjfoutfd = MakeFolder(f'{outfd}/gaussian_gjf_{keyword}_analyze_failure', allow_override=True)
            
    df = pd.read_csv(file_path_input, index_col=0)
    print(f'Total molecules {len(df)}')
    df = df[df['confid'] != 0]
    print(f'Total molecules (gaussian,success) {len(df)}')

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
        
        with open(f'{fd}/log_gaussian_worker_{keyword}_analyze.txt', 'w') as f:
            print(f'Processing {cumnum+1}/{ntotal}, number {number}, smiles {smiles}', file=f)
            
        results_dict = gaussian_analyze(number, smiles, glogfd, gchkfd, ggjffd, filepath, keyword, logoutfd, chkoutfd, gjfoutfd, gjf_key)

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
        
        status_df_failure = status_df[status_df['confid'] == 0]
        status_df_failure.to_csv(f'{outfd}/gaussian_{keyword}_analyze_failure.csv')
        
    print('success: check log')
        
if __name__ == '__main__':
    fd = os.path.dirname(os.path.abspath(__file__))
    
    # ========== Settings (CHANGE HERE) ==========
    # A name used to identify output directories and files in 2_1_nifrec_gaussian_optfreq.py (Modification is not mandatory)
    keyword = 'optfreq'
    
    # For check (e.g. 'M062X/Def2SVP')
    gjf_key = 'M062X/Def2SVP'
    
    # Enter the file and directory paths
    file_path_input = f'{fd}/gaussian_{keyword}_stats.csv'
    glogfd = f'{fd}/gaussian_{keyword}/gaussian_log_{keyword}'
    gchkfd = f'{fd}/gaussian_{keyword}/gaussian_chk_{keyword}'
    ggjffd = f'{fd}/gaussian_{keyword}/gaussian_gjf_{keyword}'
    file_path_output = f'{fd}/gaussian_{keyword}_parse.csv'
    # =============================================
    
    sys.stdout = open(f'{fd}/log_gaussian_{keyword}_parse.txt', 'w')
    process_rows_for_gparse(fd, file_path_input, glogfd, gchkfd, ggjffd, file_path_output, keyword, gjf_key)
    print('Finish')
    sys.stdout.close()
    
    