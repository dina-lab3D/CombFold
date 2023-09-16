import os
import shutil
import sys
from collections import defaultdict
from functools import lru_cache
from typing import Dict, Tuple, List

import Bio.PDB
import numpy as np
import scipy.spatial
import Bio.SeqIO


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


def generalize_xlinks(crosslinks: List[Tuple[Tuple[str, int], Tuple[str, int]]] , ident_chains: Dict[str, str]):
    # in case there are multiple crosslinks between identical chains, keep only one, as they are equivalent
    general_crosslinks = set()
    for (chain1, res1), (chain2, res2) in crosslinks:
        # before_size = len(general_crosslinks)
        to_add = tuple(sorted(((ident_chains[chain1], res1), (ident_chains[chain2], res2))))
        general_crosslinks.add(to_add)
        # if len(general_crosslinks) == before_size:
        #     print("not added", to_add, "based", ((chain1, res1), (chain2, res2)))
        # else:
        #     print("yes added", to_add, "based", ((chain1, res1), (chain2, res2)))
    return general_crosslinks


@lru_cache(maxsize=10)
def get_res_to_plddt(combfold_folder: str) -> Dict[Tuple[str, int], float]:
    subunit_names = [os.path.join(combfold_folder, i) for i in
                     open(os.path.join(combfold_folder, "chain.list")).read().split("\n") if i]
    res_to_plddt = {}
    for pdb_path in subunit_names:
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
        pdb_model = next(iter(pdb_struct))

        for res in pdb_model.get_residues():
            if "CA" not in res:
                continue
            res_to_plddt[(res.parent.id, res.id[1])] = res["CA"].get_bfactor()
    return res_to_plddt


def score_xlinks_based_on_plddt(corsslinks_path: str, combfold_folder: str, output_path: str):
    res_to_plddt = get_res_to_plddt(combfold_folder)

    crosslinks = [i.split() for i in open(corsslinks_path, "r").read().split("\n") if i]
    output_file = open(output_path, "w")
    for crosslink in crosslinks:
        if len(crosslink) == 7:
            output_file.write(" ".join(crosslink) + "\n")
            continue
        elif len(crosslink) == 6:
            res1, chains1, res2, chains2, min_dist, max_dist = crosslink
        elif len(crosslink) == 5:
            res1, chains1, res2, chains2, max_dist = crosslink
            min_dist = 0
        else:
            print("wrong crosslink", crosslink)
            continue
        res1, res2 = int(res1), int(res2)
        if (chains1[0], res1) not in res_to_plddt or (chains2[0], res2) not in res_to_plddt:
            print("missing res", res1, res2)
            continue
        plddt1 = res_to_plddt[(chains1[0], res1)]
        plddt2 = res_to_plddt[(chains2[0], res2)]
        score = round(((plddt1 + plddt2) / 2) / 100, 2)
        output_file.write(" ".join([str(i) for i in [res1, chains1, res2, chains2, min_dist, max_dist, score]]) + "\n")
    output_file.close()


def simulate_crosslinks(pdb_path: str, output_path: str):
    np.random.seed(0)
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    pdb_model = next(iter(pdb_struct))

    ident_chains = get_ident_chain_map_from_complex(pdb_path)

    coords = []
    identifiers = []
    for res in pdb_model.get_residues():
        if res.get_resname() == "LYS" and "CA" in res:
            coords.append(res["CA"].get_coord())
            identifiers.append((res.parent.id, res.id[1]))

    dists = scipy.spatial.distance.cdist(coords, coords)
    # turn off diagonal
    # np.fill_diagonal(dists, np.inf)

    # turn everything where i > j to np.inf
    for i in range(len(dists)):
        dists[i][i:] = np.inf

    close_residues = np.argwhere(dists < 30)
    far_residues = np.argwhere(dists > 40)
    inter_close_residues = [i for i in close_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]
    inter_far_residues = [i for i in far_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]

    print("there are", len(coords), "LYS and", len(close_residues), "close residues", len(inter_close_residues),
          "of them are inter")

    # filter crosslinks where the crosslinker will be disturbed
    filtered_res = []
    np_coords = np.array(coords)
    for res1, res2 in inter_close_residues:
        c1 = coords[res1]
        c2 = coords[res2]

        dir_vec = c2 - c1
        dir_size = np.linalg.norm(dir_vec)

        flag = True
        for i in range(3, int(dir_size) - 3, 2):
            checked_c = c1 + dir_vec * (i / dir_size)
            close_res = np.argwhere(np.linalg.norm(np_coords - checked_c, axis=1) < 1)
            if len(close_res) > 0:
                print("removing disturbed crosslink: res1", identifiers[res1], "res2", identifiers[res2], "dist",
                      dists[res1, res2], "i", i, "disturbed by", identifiers[close_res[0][0]])
                flag = False
                break
        if flag:
            filtered_res.append((res1, res2))
    print("there are", len(filtered_res), "filtered_res")

    generalized_crosslinks = generalize_xlinks([(identifiers[i[0]], identifiers[i[1]]) for i in filtered_res],
                                               ident_chains)
    print("there are", len(generalized_crosslinks), "generalized_crosslinks")

    # output all of them
    generalized_crosslinks = sorted(generalized_crosslinks)

    # output 10%
    output_file = open(output_path, "w")
    selected_crosslinks = np.array(generalized_crosslinks)
    np.random.shuffle(selected_crosslinks)
    selected_crosslinks = selected_crosslinks[:int(len(generalized_crosslinks) * 0.1)]
    print("-------------------")
    for (chain1, res1), (chain2, res2) in selected_crosslinks:
        chains1 = [k for k, v in ident_chains.items() if v == chain1]
        chains2 = [k for k, v in ident_chains.items() if v == chain2]
        # print(res1, "".join(chains1), res2, "".join(chains2), 30)
        output_file.write(f"{res1} {''.join(chains1)} {res2} {''.join(chains2)} 30\n")

    print("------------------- (2)")
    # output extra 5% False crosslinks
    # output_file = open(os.path.join(output_folder, "sampled_false_xlinks.txt"), "w")
    generalized_false_crosslinks = sorted(generalize_xlinks([(identifiers[i[0]], identifiers[i[1]])
                                                            for i in inter_far_residues],
                                                            ident_chains))
    generalized_false_crosslinks = np.array(generalized_false_crosslinks)
    np.random.shuffle(generalized_false_crosslinks)
    generalized_false_crosslinks = generalized_false_crosslinks[:int(len(selected_crosslinks) * 0.05)]
    for (chain1, res1), (chain2, res2) in generalized_false_crosslinks:
        chains1 = [k for k, v in ident_chains.items() if v == chain1]
        chains2 = [k for k, v in ident_chains.items() if v == chain2]
        # print(res1, "".join(chains1), res2, "".join(chains2), 30)
        output_file.write(f"{res1} {''.join(chains1)} {res2} {''.join(chains2)} 30\n")
    output_file.close()

    print("output selected", len(selected_crosslinks), "false", len(generalized_false_crosslinks))


def main(input_base_folder: str, output_base_folder: str, input_complexes_folder: str):
    os.makedirs(output_base_folder)
    for jobname in os.listdir(input_base_folder):
        if not os.path.isdir(os.path.join(input_base_folder, jobname)):
            continue
        shutil.copytree(os.path.join(input_base_folder, jobname), os.path.join(output_base_folder, jobname))
        input_complex_path = os.path.join(input_complexes_folder, jobname + ".pdb")
        combfold_path = os.path.join(output_base_folder, jobname, "assembly_output")
        xlinks_output_path = os.path.join(combfold_path, "xlink_consts.txt")
        simulate_crosslinks(input_complex_path, xlinks_output_path)
        score_xlinks_based_on_plddt(xlinks_output_path, combfold_path, xlinks_output_path)


if __name__ == '__main__':
    assert len(sys.argv) == 4, "Usage: <script> <input_base_folder> <output_base_folder> <input_complexes_folder>"
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]))
