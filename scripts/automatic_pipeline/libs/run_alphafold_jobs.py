from typing import List, Dict, Tuple

from automatic_pipeline.configurable import run_alphafold_on_fasta
from automatic_pipeline.libs.utils_classes import AlphaFoldJobInfo, RunAlphaFoldResult
import os


def run_alphafold_on_jobs(alphafold_jobs: List[AlphaFoldJobInfo], fastas_folder: str, af_output_folder: str) \
        -> Dict[AlphaFoldJobInfo, Tuple[RunAlphaFoldResult, List[str]]]:
    jobs_to_results = {}
    print("running", len(alphafold_jobs), "jobs")
    for af_job in alphafold_jobs:
        af_jobname = "_".join(af_job.subunit_names)

        os.makedirs(fastas_folder, exist_ok=True)
        os.makedirs(af_output_folder, exist_ok=True)
        fasta_path = os.path.join(fastas_folder, af_jobname + ".fasta")
        if not os.path.exists(fasta_path):
            open(fasta_path, "w").write(f">{af_job.get_jobname()}\n" + ":".join(af_job.sequences) + "\n")

        jobs_to_results[af_job] = run_alphafold_on_fasta(fasta_path, af_output_folder)
        print("doing job", af_job.subunit_names, af_job.merged_subunits)
    return jobs_to_results
