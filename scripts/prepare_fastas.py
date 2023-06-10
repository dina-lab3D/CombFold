import argparse
import os
from collections import defaultdict
from typing import List, Optional, Dict, Tuple
import Bio.PDB
import Bio.SeqUtils
import numpy as np
import scipy.spatial.distance

from libs.utils_classes import SubunitsInfo, SubunitName, read_subunits_info, INTERFACE_MIN_ATOM_DIST


def save_fasta(subunit_names: List[str], subunits_info: SubunitsInfo, output_folder: str):
    output_path = os.path.join(output_folder, "_".join(subunit_names) + ".fasta")
    with open(output_path, "w") as f:
        f.write(f">{'_'.join(subunit_names)}\n")
        f.write(":".join([subunits_info[subunit_name].sequence for subunit_name in subunit_names]) + "\n")


def get_fastas_for_pairs(subunits_info: SubunitsInfo, output_folder: str, max_af_size: int):
    all_subunits_names = sorted(list(subunits_info.keys()))
    for i in range(len(all_subunits_names)):
        start_from = i if len(subunits_info[all_subunits_names[i]].chain_names) > 1 else i + 1
        for j in range(start_from, len(all_subunits_names)):
            subunit_i, subunit_j = subunits_info[all_subunits_names[i]], subunits_info[all_subunits_names[j]]
            assert len(subunit_i.sequence) + len(subunit_j.sequence) <= max_af_size, \
                f"Subunits too long to do all alphafold pairs, joined size of {subunit_i.name} {subunit_j.name} is " \
                f"larger than AlphaFold max size {max_af_size}. You can divide large subunits into smaller ones," \
                f"or, use a GPU with more memory, and increase --max-af-size."
            save_fasta([all_subunits_names[i], all_subunits_names[j]], subunits_info, output_folder)


def score_pdb_pair(pair_path: str, names_by_sequences: Dict[str, str]) -> Optional[Tuple[Tuple[str, str], float]]:
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("p", pair_path)
    if len(list(pdb_struct)) > 1:
        return None
    model = next(iter(pdb_struct))
    if len(list(model)) != 2:
        return None
    chains = list(model.get_chains())
    chain1_seq = "".join([Bio.SeqUtils.seq1(res.get_resname()) for res in chains[0].get_residues()])
    chain2_seq = "".join([Bio.SeqUtils.seq1(res.get_resname()) for res in chains[1].get_residues()])
    if chain1_seq not in names_by_sequences or chain2_seq not in names_by_sequences:
        return None
    chain1_name = names_by_sequences[chain1_seq]
    chain2_name = names_by_sequences[chain2_seq]
    print("found pair", chain1_name, chain2_name, os.path.basename(pair_path))

    chain1_res = list(chains[0].get_residues())
    chain2_res = list(chains[1].get_residues())
    chain1_ca = np.array([res["CA"].get_coord() for res in chain1_res])
    chain2_ca = np.array([res["CA"].get_coord() for res in chain2_res])

    close_residues = np.argwhere(scipy.spatial.distance.cdist(chain1_ca, chain2_ca) < INTERFACE_MIN_ATOM_DIST)
    if len(close_residues) == 0:
        return None

    chain1_interface, chain2_interface = set(), set()
    for i, j in close_residues:
        chain1_interface.add(i)
        chain2_interface.add(j)

    bfactors = [chain1_res[i]["CA"].get_bfactor() for i in chain1_interface] + \
               [chain2_res[i]["CA"].get_bfactor() for i in chain2_interface]
    return (chain1_name, chain2_name), sum(bfactors) / len(bfactors)


def get_job_length(subunit_names, subunits_info):
    return sum([len(subunits_info[subunit_name].sequence) for subunit_name in subunit_names])


def get_fastas_for_groups(subunits_info: SubunitsInfo, output_folder: str, max_af_size: int, pairs_folder: str):
    groups_jobs = set()

    names_by_sequences = {s.sequence: s.name for s in subunits_info.values()}

    best_for_subunit: Dict[SubunitName, Dict[SubunitName, float]] = defaultdict(dict)
    pairs_to_use = [os.path.join(pairs_folder, i) for i in os.listdir(pairs_folder) if i.endswith(".pdb")]
    for pair_path in pairs_to_use:
        score_result = score_pdb_pair(pair_path, names_by_sequences)
        if score_result is None:
            continue
        (subunit1, subunit2), score = score_result
        best_for_subunit[subunit1][subunit2] = max(best_for_subunit[subunit1].get(subunit2, 0), score)
        best_for_subunit[subunit2][subunit1] = max(best_for_subunit[subunit2].get(subunit1, 0), score)

    for subunit_name, subunit_info in subunits_info.items():
        sorted_best_for_subunit = sorted(best_for_subunit[subunit_name].keys(),
                                         key=lambda x: best_for_subunit[subunit_name][x], reverse=True)

        print("best for", subunit_name, best_for_subunit[subunit_name])
        # with repetitions
        job = [subunit_name, sorted_best_for_subunit[0]]
        for i in range(3, 6):
            for subunit_to_add in sorted_best_for_subunit:
                if job.count(subunit_to_add) < len(subunits_info[subunit_to_add].chain_names) \
                        and get_job_length(job + [subunit_to_add], subunits_info) <= max_af_size:
                    job.append(subunit_to_add)
                    break
            if len(job) != i:
                break
            groups_jobs.add(tuple(sorted(job)))

        # without repetitions
        job = [subunit_name, sorted_best_for_subunit[0]]
        for subunit_to_add in sorted_best_for_subunit[1:]:
            job.append(subunit_to_add)
            if get_job_length(job, subunits_info) > max_af_size:
                break
            if len(job) > 5:
                break
            groups_jobs.add(tuple(sorted(job)))

    for job in groups_jobs:
        save_fasta(list(job), subunits_info, output_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=str, choices=["pairs", "groups"], default="pairs")
    parser.add_argument("subunits_json", type=str)
    parser.add_argument("--output-fasta-folder", type=str)
    parser.add_argument("--max-af-size", type=int, default=1800)
    parser.add_argument("--input-pairs-results", type=str, default="")
    args = parser.parse_args()

    subunits_info = read_subunits_info(args.subunits_json)

    if os.path.exists(args.output_fasta_folder):
        print("Output folder already exists, exiting")
        return
    os.makedirs(args.output_fasta_folder, exist_ok=True)

    if args.stage == "pairs":
        get_fastas_for_pairs(subunits_info, args.output_fasta_folder, args.max_af_size)
    elif args.stage == "groups":
        assert args.input_pairs_results != "" and os.path.isdir(args.input_pairs_results), \
            "When running stage=groups, must supply pairs results folder with --input-pairs-results"
        get_fastas_for_groups(subunits_info, args.output_fasta_folder, args.max_af_size, args.input_pairs_results)


if __name__ == "__main__":
    main()
