import os
import shutil
from typing import List, Dict, Tuple

from automatic_pipeline.configurable import AF2TRANS_BIN_PATH, is_bash_running, COMB_ASSEMBLY_BIN_PATH, run_bash_file
from automatic_pipeline.libs.utils_classes import AFResultScoredPair, SubunitsInfo, SubunitName
from automatic_pipeline.libs.utils_pdb import extract_pdb_info, copy_pdb_rename_chain, copy_pdb_set_start_offset


def extract_ref_structs(af_scored_pairs: List[AFResultScoredPair], output_path: str, subunits_info: SubunitsInfo):
    all_subunit_names = sum([i.get_chained_names() for i in subunits_info.values()], [])
    if all([os.path.exists(os.path.join(output_path, f"{i}.pdb")) for i in all_subunit_names]):
        return
    os.makedirs(output_path, exist_ok=True)

    ref_structs: Dict[SubunitName, Tuple[float, AFResultScoredPair, int]] = {}
    for result in af_scored_pairs:
        score = result.subunit1_scores.plddt_avg
        if result.subunits_names[0] not in ref_structs or ref_structs[result.subunits_names[0]][0] < score:
            ref_structs[result.subunits_names[0]] = (score, result, 0)

        score = result.subunit2_scores.plddt_avg
        if result.subunits_names[1] not in ref_structs or ref_structs[result.subunits_names[1]][0] < score:
            ref_structs[result.subunits_names[1]] = (score, result, 1)

    assert all(i in ref_structs for i in subunits_info.keys()), f"Missing ref structs {list(ref_structs.keys())} " \
                                                                f"{list(subunits_info.keys())}"
    for subunit_name, (score, result, subunit_idx) in ref_structs.items():
        subunit_info = subunits_info[subunit_name]
        subunit_pdb_info = result.subunit1_pdb_info if subunit_idx == 0 else result.subunit2_pdb_info
        print(f"ref_struct {subunit_name} {score} {subunit_pdb_info}")

        ref_struct_path = os.path.join(output_path, f"{subunit_name}.pdb")
        extract_pdb_info(result.pdb_path, subunit_pdb_info, ref_struct_path)
        copy_pdb_set_start_offset(ref_struct_path, subunit_info.start_res - 1, ref_struct_path)
        for chain_name, ident_subunit_name in zip(subunit_info.chain_names, subunit_info.get_chained_names()):
            copy_pdb_rename_chain(ref_struct_path, chain_name,
                                  os.path.join(output_path, f"{ident_subunit_name}.pdb"))
        os.remove(ref_struct_path)


def _score_result_transform(result: AFResultScoredPair) -> float:
    subunit1_score = result.subunit1_scores
    subunit2_score = result.subunit2_scores
    interaction_score = result.interaction_scores
    full_pae_avg: float = (interaction_score.pae_avg + subunit1_score.self_pae_avg + subunit2_score.self_pae_avg) / 3

    # when full_pae_avg=20, we will return 1. has quadric properties
    return max([1, 100 - (1/4) * (full_pae_avg ** 2)])


def create_transformations(af_scored_pairs: List[AFResultScoredPair], ref_structs_folder: str,
                           subunits_info: SubunitsInfo, transformations_output_folder: str):
    os.makedirs(transformations_output_folder, exist_ok=True)
    os.chdir(transformations_output_folder)

    temp_dirname = "temp_pdbs_dir"
    temp_dir = os.path.join(transformations_output_folder, temp_dirname)

    total_pairs_used = 0
    ordered_subunits = sorted(list(subunits_info.keys()))
    for i in range(len(ordered_subunits)):
        for j in range(i, len(ordered_subunits)):
            os.makedirs(temp_dir, exist_ok=True)
            subunit_name1 = ordered_subunits[i]
            subunit_name2 = ordered_subunits[j]
            chained_subunit_name1 = subunits_info[subunit_name1].get_chained_names()[0]
            chained_subunit_name2 = subunits_info[subunit_name2].get_chained_names()[0]

            if i == j and len(subunits_info[subunit_name1].chain_names) == 1:
                continue

            output_file_path = os.path.join(transformations_output_folder,
                                            f"{subunit_name1}_plus_{subunit_name2}")
            if all([os.path.exists(os.path.join(transformations_output_folder,
                                                f"{chained_subunit_name1}_plus_{chained_subunit_name2}"))
                    for chained_subunit_name1 in subunits_info[subunit_name1].get_chained_names()
                    for chained_subunit_name2 in subunits_info[subunit_name2].get_chained_names()]):
                continue

            results_for_pair = [result for result in af_scored_pairs
                                if result.subunits_names == (subunit_name1, subunit_name2)]
            if len(results_for_pair) == 0:
                print("No results for pair", subunit_name1, subunit_name2)
                continue
            total_pairs_used += len(results_for_pair)

            cmd = AF2TRANS_BIN_PATH + " "
            ref_subunit1_path = os.path.join(ref_structs_folder, chained_subunit_name1 + ".pdb")
            ref_subunit2_path = os.path.join(ref_structs_folder, chained_subunit_name2 + ".pdb")
            cmd += f"{ref_subunit1_path} {ref_subunit2_path} "

            results_for_pair.sort(key=_score_result_transform, reverse=True)
            scores = []
            for counter, result in enumerate(results_for_pair):
                chain1_pdb = os.path.join(temp_dir, f"{counter}_0_{os.path.basename(result.pdb_path)}")
                chain2_pdb = os.path.join(temp_dir, f"{counter}_1_{os.path.basename(result.pdb_path)}")
                extract_pdb_info(result.pdb_path, result.subunit1_pdb_info, chain1_pdb)
                extract_pdb_info(result.pdb_path, result.subunit2_pdb_info, chain2_pdb)

                cmd += f"{chain1_pdb} {chain2_pdb} "
                scores.append(_score_result_transform(result))
            print(f"result count for {subunit_name1} {subunit_name2} {len(results_for_pair)}")
            os.system(f"{cmd} > {output_file_path}")
            print(f"Saved transformations to {output_file_path}")

            # Calculate and set score based on PAE for transformations
            content_lines = open(output_file_path, "r").read().split("\n")
            scored_lines = []
            assert len(content_lines) == len(scores) + 1
            for line, score in zip(content_lines, scores):
                if len(line) == 0:
                    continue
                splitted = line.split(" | ")
                scored_lines.append(" | ".join([splitted[0], str(score)] + splitted[2:]))
            open(output_file_path, "w").write("\n".join(scored_lines) + "\n")

            # copy results for ident chains
            aliases_c1 = subunits_info[subunit_name1].get_chained_names()
            aliases_c2 = subunits_info[subunit_name2].get_chained_names()

            for c1 in range(len(aliases_c1)):
                start_from = c1 + 1 if subunit_name1 == subunit_name2 else 0
                for c2 in range(start_from, len(aliases_c2)):
                    alias_output_path = os.path.join(transformations_output_folder,
                                                     f"{aliases_c1[c1]}_plus_{aliases_c2[c2]}")
                    shutil.copy(output_file_path, alias_output_path)
            os.remove(output_file_path)
            shutil.rmtree(temp_dir)
    print("used a total of", total_pairs_used, "pairs out of ", len(af_scored_pairs))


def run_assembly(af_scored_pairs: List[AFResultScoredPair], subunits_info: SubunitsInfo, output_path: str) -> str:
    assembly_output_path = os.path.join(output_path, "assembly_output")
    transformations_folder = os.path.join(output_path, "transformations")

    os.makedirs(assembly_output_path, exist_ok=True)
    os.makedirs(transformations_folder, exist_ok=True)
    run_file_path = os.path.join(assembly_output_path, "run.sh")

    done_file_path = os.path.join(assembly_output_path, "assembly_done")
    if os.path.exists(done_file_path):
        return "success"

    if is_bash_running(run_file_path):
        return "running"

    # extract ref structs
    extract_ref_structs(af_scored_pairs, assembly_output_path, subunits_info)

    # extract transformations
    create_transformations(af_scored_pairs, assembly_output_path, subunits_info, transformations_folder)

    # prepare hasp input files
    with open(os.path.join(assembly_output_path, "chain.list"), "w") as f:
        sorted_all_subunits = sorted(sum([i.get_chained_names() for i in subunits_info.values()], []))
        for chained_subunit_name in sorted_all_subunits:
            f.write(f"{chained_subunit_name}.pdb\n")
    os.chdir(assembly_output_path)
    open("xlink_consts.txt", "w").close()

    with open(run_file_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"cd {assembly_output_path}\n")
        f.write(f"{COMB_ASSEMBLY_BIN_PATH} chain.list {transformations_folder}/ 900 100 xlink_consts.txt "
                f"-b 0.05 -t 80\n")
        f.write(f"touch {done_file_path}\n")
    run_bash_file(run_file_path)
    print("Started assembly")
    return "running"


def get_assembly_results(output_path: str) -> float:
    clusters_path = os.path.join(output_path, "assembly_output", "output_clustered.res")
    if not os.path.exists(clusters_path):
        return 0

    results_as_str = [i for i in open(clusters_path, "r").read().split("\n") if len(i) > 0]
    scores = []
    for i, result_as_str in enumerate(results_as_str):
        splitted_result = result_as_str.split(" ")
        # score = float(splitted_result[splitted_result.index("transScore_") + 1])
        scores.append(float(splitted_result[splitted_result.index("weightedTransScore") + 1]))
    return max(scores, default=0)
