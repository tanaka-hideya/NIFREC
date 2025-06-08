# NIFREC library 
XXX Brief explanation of the libray 

## Authors
- Hideya Tanka @ Nara Institute of Science and Technology (Author)
- Tomoyuki Miyao @ Nara Institute of Science and Technology (Contributor)

# Installation
Conda creates virtual enviroment (nifrec) to run the codes. Gaussian16 is installed and paths are set approprietely. 

0. Update conda 
```bash 
conda update -n base -c conda-forge conda
```

1. Create virtual environment
```bash
git clone https://github.com/tanaka-hideya/NIFREC.git
cd NIFREC
conda env create -f environment.yml
```

2. Activate the `nifrec` environment and add the script folder to the `PYTHOPATH`
```bash
conda activate nifrec
NIFRECDIR=$(pwd)
export PYTHONPATH="$PYTHONPATH:$NIFRECDIR/src/nifrec"
```

3. Test whether modules can be loaded.
```bash
python -m 1_1_nifrec_rdkit --help
```
Then, following output shows

```
usage: nifrec_rdkit [-h] [--infile INFILE] [--smicol SMICOL] [--idxcol IDXCOL] [--outfolder OUTFOLDER] [--outfile OUTFILE]
                    [--nconfs NCONFS] [--rmsd-thres RMSD_THRES] [--njobs NJOBS] [--backend BACKEND]

Conformer generator using RDKit

options:
  -h, --help            show this help message and exit
  --infile INFILE       Input csv filename
  --smicol SMICOL       Smiles column for structure generation
  --idxcol IDXCOL       Index column idx (0 starts)
  --outfolder OUTFOLDER
                        Output folder where conformers are stored
  --outfile OUTFILE     Output csv filename without extension (i.e. csv)
  --nconfs NCONFS       Number of maximum conformers to be generated
  --rmsd-thres RMSD_THRES
                        RMSD threshold between generated conformers
  --njobs NJOBS         Number of parallel computation threads (-1 means all cores)
  --backend BACKEND     Parallel computation backend (loky or multiprocessing)
```

# Running geometry optimization using an example file.
From the SMILES strings in `data/sample.csv`, conformer generation and optimization will be conducted.

```bash
mkdir results
```

## Step 1: RDkit conformre generation
```
python -m 1_1_nifrec_rdkit --infile data/sample.csv --smicol smi --outfolder results/rdkit-confgen --backend multiprocessing
```
Then, you can find the `rdkit-confgen` folder under the `results` directory, where the conformation generation results are stored.





