import os
import sys
from collections import Counter

from libs.get_alphafold_jobs import get_af_jobs_for_pairs, get_af_jobs_for_groups
from libs.parse_alphafold_jobs import score_af_results_as_pairs
from libs.run_alphafold_jobs import run_alphafold_on_jobs
from libs.run_assembly import run_assembly, get_assembly_results
from libs.utils_classes import RunAlphaFoldResult, read_subunits_info


THIS_SCRIPT_PATH = os.path.abspath(__file__)
BASE_PATH = os.path.dirname(THIS_SCRIPT_PATH)
SUBUNITS_INFO_PATH = os.path.join(BASE_PATH, "subunits_info")
MAIN_OUTPUT_PATH = os.path.join(BASE_PATH, "output")


def run_on_subunits_json(subunits_json_path: str, output_path: str):
    subunits_info = read_subunits_info(subunits_json_path)

    # run alphafold pairs
    print("-- Pairs run AlphaFold")
    pairs_af_output_path = os.path.join(output_path, "alphafold_pairs")
    pairs_fasta_path = os.path.join(output_path, "fastas_pairs")

    pairs_alphafold_jobs = get_af_jobs_for_pairs(subunits_info)
    alphafold_results = run_alphafold_on_jobs(pairs_alphafold_jobs, pairs_fasta_path, pairs_af_output_path)

    results_counter = Counter([i[0] for i in alphafold_results.values()])
    if results_counter[RunAlphaFoldResult.RUNNING] > 0:
        print(f"still running", results_counter)
        return "af_pairs"

    # parse pairs
    print("-- Pairs Parse AlphaFold")
    parsed_pairs_path = os.path.join(output_path, "parsed_pairs.pkl")
    parsed_pairs = score_af_results_as_pairs({af_job: pdb_paths for af_job, (_, pdb_paths) in
                                              alphafold_results.items()}, parsed_pairs_path, subunits_info)

    # run alphafold groups
    print("-- Groups run AlphaFold")
    groups_af_output_path = os.path.join(output_path, "alphafold_groups")
    groups_fasta_path = os.path.join(output_path, "fastas_groups")

    groups_alphafold_jobs = get_af_jobs_for_groups(subunits_info, parsed_pairs)
    print("Created groups", [i.subunit_names for i in groups_alphafold_jobs])

    return "tmp_no_groups"

    alphafold_results = run_alphafold_on_jobs(groups_alphafold_jobs, groups_fasta_path, groups_af_output_path)

    results_counter = Counter([i[0] for i in alphafold_results.values()])
    if results_counter[RunAlphaFoldResult.RUNNING] > 0:
        print(f"still running", results_counter)
        return "af_groups"

    # parse groups
    print("-- Groups Parse AlphaFold")
    parsed_groups_path = os.path.join(output_path, "parsed_groups.pkl")
    parsed_groups = score_af_results_as_pairs({af_job: pdb_paths for af_job, (_, pdb_paths) in
                                               alphafold_results.items()}, parsed_groups_path, subunits_info)

    # run Hierarchical Assembly of Predicted Subunits (HAPS)
    print("-- Assembly")
    assembly_path = os.path.join(output_path, "assembly")
    assembly_result = run_assembly(parsed_pairs + parsed_groups, subunits_info, assembly_path)

    if assembly_result == "running":
        print(f"still running")
        return "assembly"
    elif assembly_result == "success":
        print("done_" + str(get_assembly_results(assembly_path)))
        return "done_" + str(get_assembly_results(assembly_path))

    return "done_0"


def run_on_folders():
    progression = {"done": [], "af_pairs": [], "af_groups": [], "assembly": []}
    scores = []
    for subunits_info_json_path in os.listdir(SUBUNITS_INFO_PATH):
        jobname = subunits_info_json_path.split(".")[0]
        print("--------", jobname, "--------")
        stage = run_on_subunits_json(os.path.join(SUBUNITS_INFO_PATH, subunits_info_json_path),
                                     os.path.join(MAIN_OUTPUT_PATH, jobname))
        if stage.startswith("done_"):
            scores.append((jobname, float(stage.split("_")[1])))
            stage = "done"

        progression[stage].append(jobname)
    print({f"{k}_{len(v)}": v for k, v in progression.items()})
    print(scores)


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == "all":
        run_on_folders()
    elif len(sys.argv) != 3:
        print("usage: <script> subunits_info output_path")
    else:
        run_on_subunits_json(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]))
