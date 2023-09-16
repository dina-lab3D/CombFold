import os
import time
from collections import defaultdict
from typing import Dict, Tuple, List

import Bio.PDB
import Bio.SeqIO
import numpy as np
import scipy.spatial
from scipy.spatial import KDTree

INTERFACE_MIN_ATOM_DIST = 5


def get_chain_to_seq(pdb_path: str) -> Dict[str, str]:
    chain_to_seq = {str(record.id): str(record.seq) for record in Bio.SeqIO.parse(pdb_path, 'pdb-seqres')}
    if len(chain_to_seq) > 0:
        return chain_to_seq


def create_ident_chain_map(chain_to_seq: Dict[str, str]) -> Dict[str, str]:
    seq_to_chains = defaultdict(list)
    for chain_name, seq in chain_to_seq.items():
        seq_to_chains[seq].append(chain_name)

    ident_chain_map = {}
    for v in seq_to_chains.values():
        ident_name = sorted(v)[0]
        for chain_name in v:
            ident_chain_map[chain_name] = ident_name

    return ident_chain_map


def get_ident_chain_map_from_complex(pdb_path: str) -> Dict[str, str]:
    return create_ident_chain_map(get_chain_to_seq(pdb_path))


def get_contacts(pdb_path: str, ident_chains: Dict[str, str], sample_chain_map=None, offset: int = 0,
                 bfactor_threshold=0):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    pdb_model = next(iter(pdb_struct))

    coords = []
    identifiers: List[Tuple[str, int]] = []
    for res in pdb_model.get_residues():
        if "CA" in res:
            if res["CA"].get_bfactor() < bfactor_threshold:
                continue
            for atom in res:
                # if not heavy atom continue
                if atom.element == "H":
                    continue
                coords.append(atom.get_coord())
                if sample_chain_map is None:
                    # TODO: the offset is a hack because CombFold output is 0-based
                    identifiers.append((res.parent.id, res.id[1] + offset))
                else:
                    # TODO: this can cause problems for homodimer interfaces
                    identifiers.append((sample_chain_map[res.parent.id], res.id[1]))

    dists = scipy.spatial.distance.cdist(coords, coords)
    # turn off diagonal
    # np.fill_diagonal(dists, np.inf)

    # turn everything where i > j to np.inf
    # for i in range(len(dists)):
    #     dists[i][i:] = np.inf

    # prevent clashes between residues that are close in sequence
    for i in range(len(dists)):
        dists[i][max(i - 100, 0):] = np.inf

    clashes = np.argwhere(dists < 3)
    print("there are", len(clashes), "clashes")
    # TODO: try without clashes
    # dists[dists < 3] = np.inf

    # based on PMC5949145
    close_residues = np.argwhere(dists < INTERFACE_MIN_ATOM_DIST)
    inter_close_residues = [i for i in close_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]

    interface_ident_residues = set()
    for i, j in inter_close_residues:
        identifier1, identifier2 = identifiers[i], identifiers[j]
        ident_chain1, ident_chain2 = ident_chains[identifier1[0]], ident_chains[identifier2[0]]
        ident_identifier1 = (ident_chain1, identifier1[1])
        ident_identifier2 = (ident_chain2, identifier2[1])
        if ident_identifier1 > ident_identifier2:
            ident_identifier1, ident_identifier2 = ident_identifier2, ident_identifier1
        interface_ident_residues.add((ident_identifier1, ident_identifier2))

    print("there are", len(inter_close_residues), "interface pairs and", len(interface_ident_residues),
          "interface ident pairs")

    return list(interface_ident_residues)


def get_contacts_fast(pdb_path: str, ident_chains: Dict[str, str], sample_chain_map=None, offset: int = 0,
                      bfactor_threshold=0):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    pdb_model = next(iter(pdb_struct))

    coords = []
    identifiers: List[Tuple[str, int]] = []
    for res in pdb_model.get_residues():
        if "CA" in res:
            if res["CA"].get_bfactor() < bfactor_threshold:
                continue
            for atom in res:
                # if not heavy atom continue
                # if atom.element == "H":
                #     continue
                coords.append(atom.get_coord())
                if sample_chain_map is None:
                    # TODO: the offset is a hack because CombFold output is 0-based
                    identifiers.append((res.parent.id, res.id[1] + offset))
                else:
                    # TODO: this can cause problems for homodimer interfaces
                    identifiers.append((sample_chain_map[res.parent.id], res.id[1]))

    tree = KDTree(coords)
    close_residues = tree.query_pairs(INTERFACE_MIN_ATOM_DIST, p=2)
    inter_close_residues = [i for i in close_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]

    interface_ident_residues = set()
    for i, j in inter_close_residues:
        identifier1, identifier2 = identifiers[i], identifiers[j]
        ident_chain1, ident_chain2 = ident_chains[identifier1[0]], ident_chains[identifier2[0]]
        ident_identifier1 = (ident_chain1, identifier1[1])
        ident_identifier2 = (ident_chain2, identifier2[1])
        if ident_identifier1 > ident_identifier2:
            ident_identifier1, ident_identifier2 = ident_identifier2, ident_identifier1
        interface_ident_residues.add((ident_identifier1, ident_identifier2))

    print("there are", len(inter_close_residues), "interface pairs and", len(interface_ident_residues),
          "interface ident pairs")

    return list(interface_ident_residues)


def get_all_residues(pdb_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    pdb_model = next(iter(pdb_struct))

    identifiers: List[Tuple[str, int]] = []
    for res in pdb_model.get_residues():
        if "CA" in res:
            identifiers.append((res.parent.id, res.id[1]))

    return identifiers


def get_contact_score(target_pdb_path: str, sample_pdb_path: str, sample_chain_map: Dict[str, str] = None,
                      sample_offset: int = 0):
    # described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5949145/
    ident_chains = get_ident_chain_map_from_complex(target_pdb_path)
    # assume sample is based on afm

    # target_contacts = get_contacts(target_pdb_path, ident_chains)
    # sample_contacts = get_contacts(sample_pdb_path, ident_chains, sample_chain_map, sample_offset,
    #                                bfactor_threshold=50)

    target_contacts = get_contacts_fast(target_pdb_path, ident_chains)
    sample_contacts = get_contacts_fast(sample_pdb_path, ident_chains, sample_chain_map, sample_offset,
                                        bfactor_threshold=50)

    all_target_residues = get_all_residues(target_pdb_path)
    sample_contacts = [contact for contact in sample_contacts
                       if contact[0] in all_target_residues and contact[1] in all_target_residues]
    print("after filtering, there are", len(sample_contacts), "sample contacts")

    print(sorted(target_contacts)[:100])
    print(sorted(sample_contacts)[:100])

    precision = len(set(sample_contacts).intersection(target_contacts)) / len(sample_contacts)
    recall = len(set(sample_contacts).intersection(target_contacts)) / len(target_contacts)

    f1_score = 2 * precision * recall / (precision + recall)

    return f1_score


def calc_scores():
    target_folder = "input_complexes"

    target_paths = {}
    for target_pdb_name in os.listdir(target_folder):
        if not target_pdb_name.endswith(".pdb"):
            continue
        target_paths[target_pdb_name[:4]] = os.path.join(target_folder, target_pdb_name)

    # for combfold
    pdb_id_to_sample_chain_map = {}
    sample_folder = "top1_results/benchmark2"
    sample_paths = {}
    for sample_pdb_name in os.listdir(sample_folder):
        if not sample_pdb_name.endswith(".pdb"):
            continue
        sample_pdb_path = os.path.join(sample_folder, sample_pdb_name)
        sample_paths[sample_pdb_name[:4]] = sample_pdb_path

    # for afm3
    # sample_folder = "afm_results/benchmark2"
    # pdb_id_to_sample_pdb_name = {}
    # pdb_id_to_sample_chain_map = {}
    # for sample_pdb_name in os.listdir(sample_folder):
    #     if not sample_pdb_name.endswith(".pdb"):
    #         continue
    #     pdb_id = sample_pdb_name[:4]
    #     if pdb_id not in pdb_id_to_sample_pdb_name or sample_pdb_name < pdb_id_to_sample_pdb_name[pdb_id]:
    #         pdb_id_to_sample_pdb_name[pdb_id] = sample_pdb_name
    #         end_index = sample_pdb_name.split("_").index("unrelaxed")
    #         pdb_id_to_sample_chain_map[pdb_id] = {chr(ord("A") + i): v[1]
    #                                               for i, v in enumerate(sample_pdb_name.split("_")[1:end_index])}
    #
    # sample_paths = {pdb_id: os.path.join(sample_folder, sample_pdb_name)
    #                 for pdb_id, sample_pdb_name in pdb_id_to_sample_pdb_name.items()}
    # print(pdb_id_to_sample_pdb_name, pdb_id_to_sample_chain_map)

    ics_scores = {}
    for pdb_id, target_pdb_path in target_paths.items():
        print("calculating score for", pdb_id, "...")
        if pdb_id not in sample_paths:
            ics_scores[pdb_id] = 0
            continue
        sample_pdb_path = sample_paths[pdb_id]
        # this if is to check wether the input is simple AF or CombFold
        if pdb_id in pdb_id_to_sample_chain_map:
            ics_scores[pdb_id] = get_contact_score(target_pdb_path, sample_pdb_path, pdb_id_to_sample_chain_map[pdb_id])
        else:
            ics_scores[pdb_id] = get_contact_score(target_pdb_path, sample_pdb_path, sample_offset=1)
        print(pdb_id, ics_scores[pdb_id])
    print(ics_scores)


if __name__ == '__main__':
    calc_scores()
    exit(0)

