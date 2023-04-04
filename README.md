# CombFold
The CombFold software tool is capable of combinatorially assembling multiple PDB modelscreated by AlphaFold2 Multimer(AFM) into a single complex. 
The software can be run in two different ways, as follows:
1. **Starting from PDBs**: The input for this mode of operation consists of subunit sequences and a folder containing PDBs of the subunits predicted by AFM. This mode of operation does not require a GPU and can be executed entirely on Google Colab.
2. **Starting from sequences**: The input for this mode of operation is only the subunit sequences. CombFold will first perform multiple AFM jobs (colabfold installation is required), and only then will it proceed to the assembly. To use this mode of operation, the user will need to checkout the repository and modify the configurable.py file to suit their specific environment.

For the purposes of this tutorial, we will only describe the "Starting from PDBs" mode of operation.

## Defining subunits
The first step is to divide the complex into subunits. Naively, each subbuit should simply be a complete chain in the complex.  In case a chain is long, it is required to cut it into several subunits. This can be done either naively, dividing the chain into same length subunits, or by using predictors for functional domains based on sequence. Another option is to predict disordered regions based on sequence(for example using https://iupred3.elte.hu/) and remove them and divide the sequence on these regions.


Subunit is defined by 4 fields:
- name: a unique name for the subunit
- sequence: the amino acid sequence of the subunit
- chain_names: a list of chain names representing also the stoichiometry of the subunit
- start_res: the index of the start residue of the sequence on the chain. Needed to set constraints on other subunits on the same chains.

for example:
```
{
  "name": "AD1",
  "chain_names": ["A", "B"],
  "start_res": 20,
  "sequence": "LTAAAQALDGLGDKFGRSIVDGNAILADVNPRMPQIRRDITGLANLGEVY"
}
```
which defines a subunit named AD1 with 50 amino acids (the sequence length) and that has 2 copies in the complex (on chains A and B). 


Notice that each unique sequence should appear in only a single subunit definition, which can be translated into multiple chains in the assembled complex, according to stochiometry.

## Running the Combinatorial Assembly
To run the combinatorial assembly, the user must use the Google Colab notebook provided at the following link: https://colab.research.google.com/drive/1r2T2f1IDyNaloYa0URsAmuslo5z8ahwU. The input data for the assembly must be placed in a folder in Google Drive containing a file named `subunits.json` and a folder named `pdbs` that contains all PDB files. The subunits file should be a dictionary where the keys represent subunit names, and the values are the subunit definitions. A sample subunits file is available in the `example` folder of this repository.

After installing the necessary dependencies by running the first cell of the notebook, the user must locate the input data folder on their Google Drive using the `Files` toolbar on the left of the Colab notebook. The user must then copy the path and enter it into the second cell. Finally, the third cell should be executed to generate the assembled complexes, which will be saved to the same folder under a new folder named "assembled". 


It should be noted that, unlike the AFM predictions, the Colab notebook does not require scores json files produced by ColabFold for each prediction, and so instead of using PAE-based scoring an interface plDDT-based scoring is used.
