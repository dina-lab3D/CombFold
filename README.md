# CombFold
The CombFold software tool is capable of combinatorially assembling multiple PDB models created by AlphaFold2-Multimer(AFM) into a single complex. 
The software can run in two different ways, as follows:
1. **Starting from PDBs**: The input for this mode of operation consists of subunit sequences and a folder containing PDBs of the subunits predicted by AFM. This mode of operation does not require a GPU and can also be executed entirely on Google Colab.
2. **Starting from sequences**: The input for this mode of operation is only the subunit sequences. CombFold will first perform multiple AFM jobs (ColabFold installation is required), and only then will it proceed to the assembly. To use this mode of operation, the user will need to checkout the repository and modify the configurable.py file to suit their specific environment.

# System Requirements
## Hardware requirements
For the combinatorial assembler, a standard computer with a standard CPU would suffice.

To run the complete method starting from sequences, a GPU is recommended for running the ColabFold pipeline.

## Software requirements
The Combinatorial assembler is written in C++, while the supporting scripts required for running are written in Python3. Therefore, both a g++ compiler and a Python3 interpreter are required.

### OS Requirements
`CombFold` is supported for *macOS* and *Linux*. The code has been tested on the following systems:
+ macOS: Ventura (13.3.1)
+ Linux: Debian 10

### C++ Dependencies 
The Combinatorial assembler depends on Boost which can be installed using:
```
# For linux
sudo apt-get install libboost-all-dev

# For MacOS
brew install boost
```
In case boost is not in the compiler default include path, modify the BOOST_INCLUDE and BOOST_LIB variables on top of the Makefile for the installation destination.

### Python Dependencies
The supporting scripts depend on these Python3 packages, which can be installed using pip:
```
numpy
biopython
scipy
```

# Installation Guide (locally):
(Not required if you plan to use the Colab notebook for assembly)
```
git clone https://github.com/dina-lab3D/CombFold.git
cd CombinatorialAssembler
make
```
Installation should take ~3 minutes

# Running Assembly starting from PDBs
In this use case, the user already has AlphaFold-Multimer PDB results for several different combinations of the subunits. For example for a complex with 27 chains (A9B9C9), there will be several PDBs with 2 chains (AA,AB,AC,BB,BC,CC) and optionally extended subcomplexes such as ABC,AAB,AABC, etc. The input for the combinatorial assembly will be those PDB files in addition to a `subunits.json` that defines the sequences of each **unique** subunit, along with the number of copies of this subunit. The assembly algorithm can run locally, or, after uploading all PDBs and `subunits.json` to Google Drive, run using Google Colab.

It should be noted that, unlike when starting from sequences, since scores JSON files produced by ColabFold are not required, instead of using PAE-based scoring an interface plDDT-based scoring is used.

## Demo
Using the Google Colab notebook provided at the following link: https://colab.research.google.com/github/dina-lab3D/CombFold/blob/master/CombFold.ipynb.
This demo Colab notebook run the assembly on the `example` folder in the repository. It should take ~10 minutes, including installation. The expected output can be found under `example/assembled`.


## Defining subunits
The first step is to divide the complex into subunits. Naively, each subunit should simply be a complete chain in the complex. In case a chain is long, it is required to cut it into several subunits. This can be done either naively, by dividing the chain into same-length subunits, or by using predictors for functional domains based on sequence. Another option is to predict disordered regions based on sequence(for example using https://iupred3.elte.hu/) and remove them and split the sequence on these regions.

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
which defines a subunit named AD1 with 50 amino acids (the sequence length) and that has 2 copies in the complex (chains labeled A and B). 


Notice that each unique sequence should appear in only a single subunit definition, which can be translated into multiple chains in the assembled complex, according to stochiometry.

## Running the Combinatorial Assembly
### Using local installation
By running `scripts/run_on_pdbs.py`:
```
python3 scripts/run_on_pdbs.py <path_to_subunits.json> <path_to_folder_of_pdbs> <path_to_empty_output_folder>
```

### Using Google Colab
By using the Demo Colab notebook, the user just need to create the input data for the thier complex in a folder in thier Google Drive. That is a folder with a file named `subunits.json` and a folder named `pdbs` that contains all PDB files. 

After installing the necessary dependencies by running the first cell of the notebook (and ignoring the `view example elements` cell), the user must locate the input data folder on their Google Drive using the `Files` toolbar on the left of the Colab notebook. The user must then copy the path and enter it into the `Run` cell and run it. Running this cell will generate the assembled complexes (in a PDB format), which will be saved to the base input folder under a new folder named "assembled". 


