import pickle
from functools import lru_cache
import os
from typing import Dict, List
import numpy as np

from automatic_pipeline.configurable import get_scores_from_path
from automatic_pipeline.libs.get_alphafold_jobs import get_residue_diff
from automatic_pipeline.libs.utils_classes import PdbPath, AlphaFoldJobInfo, AFResultScoredPair, SubunitsInfo, \
    SubunitPdbInfo, AFSubunitScores, AFInteractionScores
from automatic_pipeline.libs.utils_pdb import get_pdb_model_readonly, are_subunits_close_in_pdb, get_interface_res_ids

try:
    import ujson as json
except ModuleNotFoundError:
    print("using slow json!")
    import json


@lru_cache(5)
def get_alphafold_scores(pdb_path: str):
    return json.load(open(get_scores_from_path(pdb_path), "r"))


def parse_results_to_scored_pairs(alphafold_results: Dict[AlphaFoldJobInfo, List[PdbPath]],
                                  subunits_info: SubunitsInfo) -> List[AFResultScoredPair]:
    all_pairs = []
    for job_count, (af_job_info, pdb_paths) in enumerate(alphafold_results.items()):
        if job_count % 5 == 0:
            print("finding pairs", job_count, "out of", len(alphafold_results), "found", len(all_pairs))
        for pdb_path in pdb_paths:
            pdb_model = get_pdb_model_readonly(pdb_path)
            chain_names = sorted([c.id for c in pdb_model.get_chains()])

            subunits_per_chain = [[af_job_info.subunit_names[0]]]
            for subunit_name, should_merge in zip(af_job_info.subunit_names[1:], af_job_info.merged_subunits):
                if should_merge:
                    subunits_per_chain[-1].append(subunit_name)
                else:
                    subunits_per_chain.append([subunit_name])
            assert len(chain_names) == len(subunits_per_chain), f"Mismatch subunits and chains lengths " \
                                                                f"{len(chain_names)} {len(subunits_per_chain)}"

            subunits_pdb_infos = []
            pdb_start_res = 0
            for chain_name, subunit_list in zip(chain_names, subunits_per_chain):
                chain_start_res = 1

                for subunit_i in range(len(subunit_list)):
                    subunit_name = subunit_list[subunit_i]
                    subunit = subunits_info[subunit_name]
                    if subunit_i > 0:
                        prev_subunit = subunits_info[subunit_list[subunit_i - 1]]
                        chain_start_res += get_residue_diff(prev_subunit, subunit)
                    subunits_pdb_infos.append(SubunitPdbInfo(chain_id=chain_name, chain_residue_id=chain_start_res,
                                                             pdb_residue_id=pdb_start_res,
                                                             length=len(subunit.sequence)))
                    chain_start_res += len(subunit.sequence)
                    pdb_start_res += len(subunit.sequence)

            for i in range(len(subunits_pdb_infos)):
                for j in range(i + 1, len(subunits_pdb_infos)):
                    if are_subunits_close_in_pdb(pdb_path, subunits_pdb_infos[i], subunits_pdb_infos[j]):
                        all_pairs.append(AFResultScoredPair(pdb_path=pdb_path,
                                                            subunits_names=(af_job_info.subunit_names[i],
                                                                            af_job_info.subunit_names[j]),
                                                            subunit1_pdb_info=subunits_pdb_infos[i],
                                                            subunit2_pdb_info=subunits_pdb_infos[j]))

    return all_pairs


def _get_subunit_scores(plddt: np.ndarray, pae: np.ndarray, subunit_indexes: List[int], interface_indexes: List[int]) \
        -> AFSubunitScores:
    subunit_plddt = plddt[subunit_indexes]
    return AFSubunitScores(
        plddt_avg=float(np.mean(subunit_plddt)),
        plddt_percentile=list(np.percentile(subunit_plddt, np.arange(0, 101, 10))),
        plddt_interface_avg=float(np.mean(plddt[interface_indexes])),
        plddt_interface_percentile=list(np.percentile(plddt[interface_indexes], np.arange(0, 101, 10))),
        self_pae_avg=float(np.mean(pae[subunit_indexes][:, subunit_indexes])),
        self_pae_percentile=list(np.percentile(pae[subunit_indexes][:, subunit_indexes], np.arange(0, 101, 10)))
    )


def _get_interaction_scores(pae: np.ndarray, subunit_indexes1: List[int], subunit_indexes2: List[int],
                            interface_indexes1: List[int], interface_indexes2: List[int]) -> AFInteractionScores:
    all_pae = np.concatenate([pae[subunit_indexes1][:, subunit_indexes2].flatten(),
                              pae[subunit_indexes2][:, subunit_indexes1].flatten()])
    pae_joined_interface = np.concatenate([pae[interface_indexes1][:, interface_indexes2].flatten(),
                                           pae[interface_indexes2][:, interface_indexes1].flatten()])
    return AFInteractionScores(
        pae_avg=float(np.mean(all_pae)),
        pae_percentile=list(np.percentile(all_pae, np.arange(0, 101, 10))),
        pae_joined_interface_avg=float(np.mean(pae_joined_interface)),
        pae_joined_interface_percentile=list(np.percentile(pae_joined_interface, np.arange(0, 101, 10))),
        interface1_size=len(interface_indexes1),
        interface2_size=len(interface_indexes2),
    )


def score_af_results_as_pairs(alphafold_results: Dict[AlphaFoldJobInfo, List[PdbPath]],
                              output_path: str, subunits_info: SubunitsInfo) -> List[AFResultScoredPair]:
    if os.path.exists(output_path):
        with open(output_path, "rb") as f:
            all_pairs = pickle.load(f)
    else:
        all_pairs = parse_results_to_scored_pairs(alphafold_results, subunits_info)
        pickle.dump(all_pairs, open(output_path, "wb"))

    not_scored_results: List[AFResultScoredPair] = [i for i in all_pairs if i.subunit1_scores is None]

    # enrich with scores
    for count, pair in enumerate(not_scored_results):
        # pdb_info = AlphaFoldJobInfo.from_result_path(pair.pdb_path)

        job_scores = get_alphafold_scores(pair.pdb_path)
        plddt = np.array(job_scores["plddt"])
        pae = np.array(job_scores["pae"])

        # index are 0-based
        subunit1_indexes = [i for i in range(pair.subunit1_pdb_info.length)
                            if subunits_info[pair.subunits_names[0]].sequence[i] != "X"]
        subunit1_indexes = np.array(subunit1_indexes) + pair.subunit1_pdb_info.pdb_residue_id

        subunit2_indexes = [i for i in range(pair.subunit2_pdb_info.length)
                            if subunits_info[pair.subunits_names[1]].sequence[i] != "X"]
        subunit2_indexes = np.array(subunit2_indexes) + pair.subunit2_pdb_info.pdb_residue_id

        subunit1_interface, subunit2_interface = get_interface_res_ids(pair.pdb_path,
                                                                       pair.subunit1_pdb_info,
                                                                       pair.subunit2_pdb_info)
        subunit1_interface_indexes = np.array(list(subunit1_interface)) + pair.subunit1_pdb_info.pdb_residue_id
        subunit2_interface_indexes = np.array(list(subunit2_interface)) + pair.subunit2_pdb_info.pdb_residue_id

        pair.subunit1_scores = _get_subunit_scores(plddt, pae, subunit1_indexes, subunit1_interface_indexes)
        pair.subunit2_scores = _get_subunit_scores(plddt, pae, subunit2_indexes, subunit2_interface_indexes)
        pair.interaction_scores = _get_interaction_scores(pae, subunit1_indexes, subunit2_indexes,
                                                          subunit1_interface_indexes, subunit2_interface_indexes)
        if count % 100 == 0 or count == len(not_scored_results) - 1:
            print(f"Processed {count + 1}/{len(not_scored_results)} total pairs: {len(all_pairs)}")
            pickle.dump(all_pairs, open(output_path, "wb"))

    return all_pairs
