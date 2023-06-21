# CombFold
The `CombFold` pipeline predicts the structure of large protein complexes starting from the sequences of chains in the complex (up to at least 18,000 amino acids and 32 subunits). 
The pipeline uses AlphaFold-Multimer (AFM) to predict structures of "possible subcomplexes" which are combinations of subunits from the target complex. The **CombFold Combinatorial Assembly** algorithm assembles those structures into a single large complex.

The pipeline has 4 stages:
1. Defining subunits in the complex
2. Predicting structures using AFM for all pairings of subunits
3. \[Optional\] Predicting structures using AFM for larger groups of subunits
4. Running the Combinatorial Assembly algorithm on all generated structures.


## Demo

To view a demonstration of inputs given to CombFold and a run of the assembly algorithm, use the [Demo Google Colab notebook](https://colab.research.google.com/github/dina-lab3D/CombFold/blob/master/CombFold.ipynb). 

This demo Colab notebook runs the assembly on the `example` folder in the repository. It should take ~10 minutes, including installation. The expected output can be found under `example/expected_assembled`.

If you already have some AFM predictions for possible subcomplexes for some target complex, you can also use the Colab Notebook to run the assembly for your target complex without any local installation.


# Installation
## System Requirements
### Hardware requirements
For the combinatorial assembler, a standard computer with a standard CPU would suffice.

For generating AlphaFold-Multimer predictions locally (which is likely required for heteromeric complexes), a GPU with at least 12GB of memory is recommended.

### Software requirements
The Combinatorial assembler is written in C++, while the supporting scripts required for running are written in Python3. Therefore, both a g++ compiler and a Python3 interpreter are required.

#### OS Requirements
`CombFold` is supported for *macOS* and *Linux*. The code has been tested on the following systems:
+ macOS: Ventura (13.3.1)
+ Linux: Debian 10

#### C++ Dependencies 
The Combinatorial assembler depends on Boost which can be installed using:
```
# For Linux
sudo apt-get install libboost-all-dev

# For MacOS
brew install boost
```
In case boost is not in the compiler default include path, modify the BOOST_INCLUDE and BOOST_LIB variables on top of the Makefile for the installation destination.

#### Python Dependencies
The supporting scripts depend on these Python3 packages, which can be installed using pip:
```
numpy
biopython
scipy
```

## Installation Guide (locally):
```
git clone https://github.com/dina-lab3D/CombFold.git
cd CombFold/CombinatorialAssembler
make
```
Installation should take ~3 minutes

# Running CombFold
## Stage 1 - Defining subunits
The first step is to divide the complex into subunits and create the `subunits.json` file that defines the complex. Subunits would not change their structure during the assembly (only their position relative to other subunits structures), so we would like to choose subunits that are a single structured domain.

**Naively, each subunit should simply be a complete chain in the complex**. In case a chain is long, it is required to cut it into several subunits. This can be done either naively, by dividing the chain into same-length subunits, or by using predictors for functional domains based on sequence. Another option is to predict disordered regions based on sequence(for example using [IUPred3](https://iupred3.elte.hu/)) and remove them and split the sequence on these regions.

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

**Notice** that each unique sequence should appear in only a single subunit definition, which can be translated into multiple chains in the assembled complex, according to stochiometry.

The `subunits.json` is a JSON dictionary of subunits, for example:
```
{
  "A0": {"name": "A0", "chain_names": ["A", "B"], "start_res": 1, "sequence": "MKDILEKLEERRAQARLGGGEKRLEAQHKRGKLTARERIELLLDHGSFEE"},
  "C0": {"name": "C0", "chain_names": ["C", "D"], "start_res": 1, "sequence": "MFDKILIANRGEIACRIIKTAQKMGIKTVAVYSDADRDAVHVAMADEAVH"},
  "E0": {"name": "E0", "chain_names": ["E"], "start_res": 1, "sequence": "MGDKIESKKAAAAAEVSTVPGFLGVIESPEHAVTIADEIGYPVMIKASAGA"},
  "E1": {"name": "E1", "chain_names": ["E"], "start_res": 51, "sequence": "GGGKGMRIAESADEVAEGFARAKSEASSSFGDDRVFVEKFITDPRHIEIQ"},
}
```
This describes a complex with 5 chains (A,B,C,D,E), where A & B are the same chains (length 50) and so is C & D (length 50). Additionally, in this example, the chain E (length 100) is divided into two subunits. This can happen, for example, if the complete E would have been too large to be predicted with other subunits in our GPU.

## Stage 2 - Predicting structures for pairs
In this stage, we will run AFM for every pairing of subunits.

Using the script from this repository:
```
python3 scripts/prepare_fastas.py subunits.json --stage pairs --output-fasta-folder <path_to_output_folder> --max-af-size 1800
```
will result in a folder with up to `((N+1)*N)/2` `.fasta` files. Each of these files can be used as an input for AFM.

To run AFM you can use [ColabFold](https://github.com/sokrypton/ColabFold) to run using a Google Colab Notebook. Alternatively, you can use [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) or [AlphaPullDown](https://github.com/KosinskiLab/AlphaPulldown) to run locally (which will require a GPU). 

Notice that the number of AFM predictions required is dependent on the number of subunits (*unique* chains in the complex), and not the number of chains in the complex. For example, a homooligomer with 10 chains, will only require a single AFM prediction at this stage. Therefore, for complexes with up to 3 subunits (3 unique chains), perhaps using a Colab Notebook can suffice. 

Notice that the command line defines `--max-af-size` which should be set to the maximal number of residues that can be predicted using your prediction environment (local GPU or GPU supplied by Google Colab).

If running ColabFold locally, you can use a command line similar to this for each fasta path:
```
colabfold_batch <fasta_path> <output_folder> --num-models 5
```

## Stage 3 - \[Optional\] Predicting structures for larger groups
In this stage, we will run AFM for larger groups of subunits (up to 6 subunits at a single prediction). To limit the number of required predictions, we will choose only larger groups that are more likely to give high-scored results, based on the scores of predictions of pairs. 

This stage is optional, as the assembly can be done using only the pairs predictions, however, this stage significantly improves the accuracy of generated results and the ability to assemble challenging complexes.

Using the same script, we will generate a folder with `.fasta` files of larger groups:
```
python3 scripts/prepare_fastas.py subunits.json  --stage groups --output-fasta-folder <path_to_output_folder>--max-af-size 1800 --input-pairs-results <path_to_AFM_pairs_results>
```

Here you will also need to supply path_to_AFM_pairs_results which will be a folder containing all `.pdb` files that were predicted by AFM in the previous stage.

While the script generates suggestions for groupings of subunits, to improve results, the user is also encouraged to use biological knowledge about the target complex to manually create `.fasta` files for known groupings of subunits not suggested by the script.

The generated and manually defined `.fasta` files are to be supplied to AFM as input, as done in the previous stage.

## Stage 4 - Combinatorial Assembly
In this use stage, the user already has AlphaFold-Multimer structure predictions (in the form of `.pdb` files) for several different combinations of the subunits. For example for a complex with 27 chains (A9B9C9), there will be all `.pdb` models with 2 chains (AA,AB,AC,BB,BC,CC) and optionally extended subcomplexes such as ABC,AAB,AABC, etc. The input for the combinatorial assembly will be those `.pdb` files in addition to the `subunits.json`.

Notice that unlike in the orginal paper, the scoring function for interactions used by the supplied scripts is based on interface-plDDT and not PAE.

The assembly can run either locally or in a Google Colab Notebook.

### Using Google Colab
By using the Demo Colab notebook, the user need to upload the input data for their complex into a folder in their Google Drive. That is a folder containing the file `subunits.json` and a folder named `pdbs` that contains all PDB files generated by AFM in the previous stages. 

To run the algorithm use the [Demo Google Colab notebook](https://colab.research.google.com/github/dina-lab3D/CombFold/blob/master/CombFold.ipynb).

After installing the necessary dependencies by running the first cell of the notebook (and ignoring the `view example elements` cell), the user must locate the input data folder on their Google Drive using the `Files` toolbar on the left of the Colab notebook. The user must then copy the path and enter it into the `Run` cell and run it. Running this cell will generate the assembled complexes (in a PDB or CIF format), which will be saved to the base input folder (containing `subunits.json` under a new folder named "assembled". 

### Using local installation
By running `scripts/run_on_pdbs.py`:
```
python3 scripts/run_on_pdbs.py <path_to_subunits.json> <path_to_folder_of_pdbs> <path_to_empty_output_folder>
```


