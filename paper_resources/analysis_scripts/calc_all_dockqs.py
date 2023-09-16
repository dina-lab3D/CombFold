import os
import subprocess
from typing import Tuple, List, Dict, Optional

import Bio.PDB
import Bio.SeqUtils
import numpy as np
import scipy.spatial

# prints all pairwise dockqs given a reference and a single model

INTERFACE_MIN_ATOM_DIST = 8

RMSD_PATH = "/cs/labs/dina/dina/projects/rmsd/rmsd3.linux"
# TODO: change this to your path
WORK_DIR = "~/tmp/check_dockq"


def get_pdb_model_readonly(pdb_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    return next(iter(pdb_struct))


def _get_chained_res_id(res) -> Tuple[str, int]:
    return res.parent.id, res.get_id()[1]


def _get_clean_res_ids(pdb_path) -> List[Tuple[str, int]]:
    return [_get_chained_res_id(res) for res in get_pdb_model_readonly(pdb_path).get_residues()
            if "CA" in res and Bio.SeqUtils.seq1(res.get_resname()) != "X"]


def extract_chains(pdb_path: str, output_path: str, chains_to_save: List[str]):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, "can't extract if more than one model"
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())
    for chain_id in chains_to_save:
        assert len([c for c in chains if c.id == chain_id]) == 1, f"Missing chain {chain_id} from {pdb_path}"
    for chain in chains:
        if chain.id not in chains_to_save:
            model.detach_child(chain.id)
    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def prepare_chains_pdbs(ref_path: str, sample_path: str, output_folder: str, chain_map: Optional[Dict[str, str]]) \
        -> Dict[str, Tuple[str, str]]:
    ref_pdb_model = get_pdb_model_readonly(ref_path)
    sample_pdb_model = get_pdb_model_readonly(sample_path)

    if chain_map is None:
        chain_map = {c.id: c.id for c in ref_pdb_model.get_chains() if c.id in sample_pdb_model}

    chain_to_files = {}
    for chain_name in chain_map.keys():
        print("Preparing", chain_name)
        output_pdb_path1 = os.path.join(output_folder, f"ref_{chain_name}.pdb")
        output_pdb_path2 = os.path.join(output_folder, f"sample_{chain_name}.pdb")

        chain_to_files[chain_name] = (output_pdb_path1, output_pdb_path2)

        chain_ref = ref_pdb_model[chain_name]
        chain_sample = sample_pdb_model[chain_map[chain_name]]

        chain_ref_res_ids = {res.get_id()[1] for res in chain_ref.get_residues() if "CA" in res}
        chain_sample_res_ids = {res.get_id()[1] for res in chain_sample.get_residues() if "CA" in res}
        joined_res_ids = chain_ref_res_ids.intersection(chain_sample_res_ids)

        ref_res_to_remove = []
        for res in chain_ref.get_residues():
            if res.get_id()[1] not in joined_res_ids:
                ref_res_to_remove.append(res)
        for res in ref_res_to_remove:
            res.parent.detach_child(res.id)

        sample_res_to_remove = []
        for res in chain_sample.get_residues():
            if res.get_id()[1] not in joined_res_ids:
                sample_res_to_remove.append(res)
        for res in sample_res_to_remove:
            res.parent.detach_child(res.id)

        io = Bio.PDB.PDBIO()
        io.set_structure(chain_ref)
        io.save(output_pdb_path1)

        io = Bio.PDB.PDBIO()
        io.set_structure(chain_sample)
        io.save(output_pdb_path2)

    return chain_to_files


def is_chains_coords_close(chain1_ca: np.ndarray, chain2_ca: np.ndarray):
    if len(chain1_ca) == 0 or len(chain2_ca) == 0:
        return False
    return np.min(scipy.spatial.distance.cdist(chain1_ca, chain2_ca)) < INTERFACE_MIN_ATOM_DIST


def get_interacting_chains(pdb_path) -> List[Tuple[str, str]]:
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())

    close_pairs = []
    for i in range(len(chains)):
        chain_i_ca = np.array([res['CA'].coord for res in chains[i].get_residues() if 'CA' in res])
        for j in range(i + 1, len(chains)):
            chain_j_ca = np.array([res['CA'].coord for res in chains[j].get_residues() if 'CA' in res])
            if is_chains_coords_close(chain_i_ca, chain_j_ca):
                close_pairs.append((chains[i].get_id(), chains[j].get_id()))
    return [(sorted(i)[0], sorted(i)[1]) for i in close_pairs]


def print_all_dockq(ref_path, sample_path, chain_map):
    tmpdir = os.path.join(WORK_DIR, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    chain_to_files = prepare_chains_pdbs(ref_path, sample_path, tmpdir, chain_map)
    interacting_chains = get_interacting_chains(ref_path)

    trans_path = os.path.join(tmpdir, "transformation.txt")
    open(trans_path, "w").write("1 0 0 0 0 0 0\n")

    tmp_output_path = os.path.join(tmpdir, "output.txt")

    scores = {}
    for c1, c2 in  interacting_chains:
        print("Running", c1, c2, chain_to_files[c1], chain_to_files[c2])

        subprocess.run(f"{RMSD_PATH} {chain_to_files[c1][0]} {chain_to_files[c1][1]} "
                       f"{chain_to_files[c2][0]} {chain_to_files[c2][1]} "
                       f"{trans_path} -o {tmp_output_path}", shell=True)
        dockqs = [float(line.split(" | ")[2].strip()) for line in open(tmp_output_path, "r").read().split("\n")[1:-3]
                  if len(line) > 0]
        print("DockQ:", dockqs)
        scores[(c1, c2)] = dockqs[0]

    print("Final scores:", scores)


if __name__ == '__main__':
    ref_path = "afm3_dataset/input_complexes/7ZKQ.pdb"

    sample_path = "outputClustered_6.pdb"
    chain_map = None

    # sample_path = "7ZKQ_22_AA_CC_TT_bb_unrelaxed_rank_003_alphafold2_multimer_v3_model_5_seed_000.pdb"
    # chain_map = {"2": "A", "A": "B", "C": "C", "T": "D", "b": "E"}

    print_all_dockq(ref_path, sample_path, chain_map)