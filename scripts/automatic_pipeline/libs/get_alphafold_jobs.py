from collections import defaultdict
from typing import List, Optional, Dict

from automatic_pipeline.configurable import MAX_AF_SIZE
from automatic_pipeline.libs.utils_classes import SubunitsInfo, AlphaFoldJobInfo, SubunitInfo, AFResultScoredPair, \
    SubunitName


def get_residue_diff(subunit1: SubunitInfo, subunit2: SubunitInfo) -> int:
    return (subunit2.start_res - 1) - subunit1.get_end_res()


def get_merged_info(subunits_info: SubunitsInfo, subunit_names: List[str]) -> List[bool]:
    ret_merged = []
    residues_less_than_max = MAX_AF_SIZE - sum([len(subunits_info[s].sequence) for s in subunit_names])

    for i in range(len(subunit_names) - 1):
        subunit1 = subunits_info[subunit_names[i]]
        subunit2 = subunits_info[subunit_names[i + 1]]

        residue_diff = get_residue_diff(subunit1, subunit2)
        if subunit1.chain_names[0] == subunit2.chain_names[0] \
                and 0 <= residue_diff <= min(residues_less_than_max, 100) and subunit1.name != subunit2.name:
            residues_less_than_max -= residue_diff
            ret_merged.append(True)
        else:
            ret_merged.append(False)
    return ret_merged


def get_alphafold_job(subunits_info: SubunitsInfo, subunit_names: List[str]) -> Optional[AlphaFoldJobInfo]:
    assert tuple(subunit_names) == tuple(sorted(subunit_names)), f"subunit names must be sorted {subunit_names}"
    sequences_length = sum([len(subunits_info[subunit_name].sequence) for subunit_name in subunit_names])
    if sequences_length > MAX_AF_SIZE:
        print(f"Skipping because too long", subunit_names)
        return None

    merged_subunits = get_merged_info(subunits_info, subunit_names)

    sequences = []
    next_sequence = subunits_info[subunit_names[0]].sequence

    for i in range(1, len(subunit_names)):
        if merged_subunits[i - 1]:
            residue_diff = get_residue_diff(subunits_info[subunit_names[i - 1]], subunits_info[subunit_names[i]])
            next_sequence += "X" * residue_diff
            next_sequence += subunits_info[subunit_names[i]].sequence
        else:
            sequences.append(next_sequence)
            next_sequence = subunits_info[subunit_names[i]].sequence

    if next_sequence:
        sequences.append(next_sequence)

    return AlphaFoldJobInfo(subunit_names=subunit_names,
                            merged_subunits=get_merged_info(subunits_info, subunit_names),
                            sequences=sequences)


def get_af_jobs_for_pairs(subunits_info: SubunitsInfo) -> List[AlphaFoldJobInfo]:
    returned_jobs = []
    all_subunits_names = sorted(list(subunits_info.keys()))
    for i in range(len(all_subunits_names)):
        start_from = i if len(subunits_info[all_subunits_names[i]].chain_names) > 1 else i + 1
        for j in range(start_from, len(all_subunits_names)):
            af_job = get_alphafold_job(subunits_info, [all_subunits_names[i], all_subunits_names[j]])
            if af_job is not None:
                returned_jobs.append(af_job)

    return returned_jobs


def _score_result_transform(result: AFResultScoredPair) -> float:
    full_pae_avg: float = (result.interaction_scores.pae_avg + result.subunit1_scores.self_pae_avg
                           + result.subunit2_scores.self_pae_avg) / 3

    # when full_pae_avg=20, we will return 1. has quadric properties
    pae_based_score = max([1, 100 - (1/4) * (full_pae_avg ** 2)])
    iplddt_score = (result.subunit1_scores.plddt_interface_avg + result.subunit2_scores.plddt_interface_avg) / 2
    return pae_based_score + iplddt_score/100


def get_job_length(subunit_names, subunits_info):
    return sum([len(subunits_info[subunit_name].sequence) for subunit_name in subunit_names])


def get_af_jobs_for_groups(subunits_info: SubunitsInfo, parsed_pairs: List[AFResultScoredPair])\
        -> List[AlphaFoldJobInfo]:
    groups_jobs = set()

    best_for_subunit: Dict[SubunitName, Dict[SubunitName, float]] = defaultdict(dict)
    for pair in parsed_pairs:
        subunit1, subunit2 = pair.subunits_names
        score = _score_result_transform(pair)

        best_for_subunit[subunit1][subunit2] = max(best_for_subunit[subunit1].get(subunit2, 0), score)
        best_for_subunit[subunit2][subunit1] = max(best_for_subunit[subunit2].get(subunit1, 0), score)
    best_for_subunit = dict(best_for_subunit)

    for subunit_name, subunit_info in subunits_info.items():
        sorted_best_for_subunit = sorted(best_for_subunit[subunit_name].keys(),
                                         key=lambda x: best_for_subunit[subunit_name][x], reverse=True)

        print("best for", subunit_name, best_for_subunit[subunit_name], sorted_best_for_subunit)
        # with repetitions
        job = [subunit_name, sorted_best_for_subunit[0]]
        for i in range(3, 6):
            for subunit_to_add in sorted_best_for_subunit:
                if job.count(subunit_to_add) < len(subunits_info[subunit_to_add].chain_names) \
                        and get_job_length(job + [subunit_to_add], subunits_info) <= MAX_AF_SIZE:
                    job.append(subunit_to_add)
                    break
            if len(job) != i:
                break
            print("not adding yet, with repetitions", job)
            groups_jobs.add(tuple(sorted(job)))
        # if len(job) > 2:
        #     print("adding with repetitions", job)
        #     groups_jobs.add(tuple(sorted(job)))

        # without repetitions
        job = [subunit_name, sorted_best_for_subunit[0]]
        for subunit_to_add in sorted_best_for_subunit[1:]:
            job.append(subunit_to_add)
            if get_job_length(job, subunits_info) > MAX_AF_SIZE:
                break
            if len(job) > 5:
                break
            print("not adding yet, without repetitions", job)
            groups_jobs.add(tuple(sorted(job)))
        # job = job[:-1]
        # if len(job) > 2:
        #     print("adding without repetitions", job)
        #     groups_jobs.add(tuple(sorted(job)))

    return [get_alphafold_job(subunits_info, list(job)) for job in groups_jobs]
