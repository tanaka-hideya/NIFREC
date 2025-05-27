import os
import datetime
import time
from rdkit import Chem

def GetTime(return_second=False, return_date=False):
    t = datetime.datetime.fromtimestamp(time.time())
    if return_second:
        return '{}{}{}_{}{}{}'.format(t.year, t.month, t.day, t.hour, t.minute, t.second)
    elif return_date:
        return '{}{}{}'.format(t.year,t.month,t.day)
    else:
        return '{}{}{}_{}{}'.format(t.year, t.month, t.day, t.hour, t.minute)

def search_exist_suffix(f_path):
    dirname, basename = os.path.split(f_path)

    for i in range(1,1000):
        if '.' in basename: # file
            n_name = basename.replace('.', '_{}.'.format(i))
        else:
            n_name = basename + '_' + str(i) # folder
        
        new_name = os.path.join(dirname,n_name)
        if not os.path.exists(new_name):
            return new_name
        # could not find unused filename for 1000 loops
    ValueError('cannot find unused folder')

def MakeFolder(folder_path, allow_override=False, skip_create=False, time_stamp=False):
    """
    Make a folder
    """
    today = GetTime(return_date=True)
    if os.path.exists(folder_path) and skip_create:
        if time_stamp:
            return f'{folder_path}_{today}'
        return folder_path

    if os.path.exists(folder_path) and (not allow_override):
        Warning('Specified folder already exists. Create new one')
        folder_path = search_exist_suffix(folder_path)
        
    if time_stamp:
        folder_path = f'{folder_path}_{today}'
        
    os.makedirs(folder_path, exist_ok=allow_override)

    return folder_path

def SmilesToCanSmiles(smi: str):
    # Note, RDKit successfully translate a blank sentence to mol object
    if (smi is None) or (smi == ''):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    else:
        return Chem.MolToSmiles(mol) # canonical smiles are generated

def WriteMolsSDF(mols, fname):
    w = Chem.rdmolfiles.SDWriter(fname)
    for mol in mols: 
        w.write(mol)

def GetMolsFromStream(suppl):
    fail_mols = 0
    mols  = []
    for x in suppl:
        if x is None:
            fail_mols +=1
            continue
        mols.append(x)
    return mols, fail_mols

def LoadMolsFromSDF(fname: str, return_counts: bool=False, remove_Hs: bool=True, v3000=False, strict_parsing=True):
    suppl = Chem.SDMolSupplier(fname, removeHs=remove_Hs, strictParsing=strict_parsing) 
    mols, fail_mols = GetMolsFromStream(suppl)
   
    if return_counts:
        return mols, fail_mols
    else:
        return mols

def WriteMolsToSDF(mol_list: list, sdf_name: str):  
    w = Chem.SDWriter(sdf_name)
    for m in mol_list:
        w.write(m)
    w.close()