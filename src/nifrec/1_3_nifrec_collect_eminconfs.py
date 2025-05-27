"""
File: 1_3_nifrec_collect_eminconfs.py
Author: Hideya Tanaka
Supervised by: Tomoyuki Miyao
Description:
    For each molecule, the xyz files of the conformer with the lowest xTB 
    energy among all conformers without imaginary frequencies are collected.
"""

import os
import pandas as pd
import json
import pickle
import shutil
from rdkit.Geometry import Point3D
from utility_opt import MakeFolder, LoadMolsFromSDF, WriteMolsToSDF
import sys

def setcoords_from_xyz(fname, sdffd): 
	molid 	= os.path.basename(fname).split('_')[1]
	sdfpath = f'{sdffd}/rdkit2d_{molid}.sdf'
	mol = LoadMolsFromSDF(sdfpath, remove_Hs=False)[0]
	
	# Load coordinates
	coords 		= list()
	elements 	= list()
	with open(fname, 'r') as fp:
		for line_number,line in enumerate(fp):
			if line_number == 0:
				num_atoms = int(line)
			elif line_number == 1:
				comment = line 
			else:
				element, x, y, z = line.split()
				elements.append(element)
				coords.append([float(x),float(y),float(z)])		
	
	# Set the coordinates
	conf = mol.GetConformer()
	for idx in range(mol.GetNumAtoms()):
		x,y,z 	= coords[idx]
		element = elements[idx]
		if mol.GetAtomWithIdx(idx).GetSymbol() != element:
			raise ValueError(f'Atomic Order is not consistent. index: {molid}.')
		conf.SetAtomPosition(idx, Point3D(x,y,z))
	return mol

def collect_eminconfs(fd, file_path_input, xtbxyzfd, rdkitsdffd, file_path_output):
	outfd = MakeFolder(f'{fd}/xTBopt_Emin_xyz', allow_override=True)
	df = pd.read_csv(file_path_input, index_col=0)
	print(f'Total molecules {len(df)}')
	df = df[df['confid'] != 0]
	print(f'Total molecules (xTB,success) {len(df)}')
		
	qminf_all = dict()
	mols 	  = list()
	ntotal 	  = len(df)
	for cumnum, row in enumerate(df.itertuples(index=True)):
		print(f'processing {cumnum+1}/{ntotal}', end='\r')
		number = row.Index
		filepath  = row.filepath
		xtbpath = f'{xtbxyzfd}/{filepath}'
		jsonpath =  xtbpath.replace('.xyz', '.json') # If you need property of qm calculation
		qminf 	= json.load(open(jsonpath,'r'))
		qminf_all[number] = qminf
		mol = setcoords_from_xyz(xtbpath, rdkitsdffd)
		shutil.copy(xtbpath, outfd)
		mols.append(mol)
	pickle.dump(qminf_all, open(f'{fd}/xTBopt_Emin_jsons.pickle','wb'), protocol=pickle.DEFAULT_PROTOCOL)
	WriteMolsToSDF(mols, file_path_output)

if __name__ == '__main__':
	fd = os.path.dirname(os.path.abspath(__file__))
    
	# ========== Settings (CHANGE HERE) ==========
	# Enter the file and directory paths
	file_path_input = f'{fd}/xTB_stats_Emin.csv'
	xtbxyzfd = f'{fd}/xTB/opt'
	rdkitsdffd = f'{fd}/rdkit/sdf'
	file_path_output = f'{fd}/xTBopt_Emin_mols.sdf'
	# =============================================
 
	sys.stdout = open(f'{fd}/log_collect.txt', 'w')
	collect_eminconfs(fd, file_path_input, xtbxyzfd, rdkitsdffd, file_path_output)
	print('Finish')
	sys.stdout.close()
