import re
from functools import lru_cache
import subprocess
from typing import List, Optional, Dict, Tuple
import os

from libs.utils_classes import AlphaFoldJobInfo, RunAlphaFoldResult

# Constants
MAX_AF_SIZE = 3000
AF2TRANS_BIN_PATH = "/cs/labs/dina/bshor/repos/CombFold/CombinatorialAssembler/AF2trans.out"
COMB_ASSEMBLY_BIN_PATH = "/cs/labs/dina/bshor/repos/CombFold/CombinatorialAssembler/CombinatorialAssembler.out"


# slurm utils
@lru_cache(1)
def get_squeue_output():
    return subprocess.check_output(["squeue"])


def get_slurm_ids_from_log(slurms_log: str) -> List[str]:
    if not os.path.exists(slurms_log):
        return []
    return re.findall(r"Submitted batch job (\d+)", open(slurms_log, "r").read())


def check_for_active_slurm(slurms_log: str) -> bool:
    slurm_ids = get_slurm_ids_from_log(slurms_log)
    if len(slurm_ids) == 0:
        return False
    squeue = get_squeue_output()
    for slurm_id in slurm_ids:
        for line in squeue.split(b"\n"):
            if slurm_id.encode() in line:
                return True
    return False


def get_last_slurm_output(output_folder: str, slurms_log: str) -> Optional[str]:
    slurm_ids = get_slurm_ids_from_log(slurms_log)
    if len(slurm_ids) == 0:
        return None
    if not os.path.exists(os.path.join(output_folder, f"slurm-{slurm_ids[-1]}.out")):
        return None
    return open(os.path.join(output_folder, f"slurm-{slurm_ids[-1]}.out"), "r").read()


def is_bash_running(bash_path: str):
    output_path = os.path.dirname(bash_path)
    slurm_ids_log_path = os.path.join(output_path, f"slurm_logs_{os.path.basename(bash_path)}.log")
    if check_for_active_slurm(slurm_ids_log_path):
        return True
    return False


def run_bash_file(bash_path: str):
    output_path = os.path.dirname(bash_path)
    slurm_ids_log_path = os.path.join(output_path, f"slurm_logs_{os.path.basename(bash_path)}.log")
    subprocess.run(f"sbatch -c2 --mem=10G --time=3:0:0 {bash_path} >> {slurm_ids_log_path}", shell=True)


# alphafold run
def get_scores_from_path(af_pdb_path: str):
    return af_pdb_path.replace("unrelaxed", "scores").replace("pdb", "json")


def run_alphafold_on_fasta(fasta_path: str, output_folder: str) -> Tuple[RunAlphaFoldResult, List[str]]:
    """
    This function gets a path to a fasta and an output folder and runs alphafold on the fasta.
    The function return a tuple of the status of the alphafold job(running/completed/failed) and if completed,
    the second part of the returned tuple should be a list of paths to the alphafold predicted pdbs.
    """
    af_jobname = [i for i in open(fasta_path, "r").read().split("\n") if i.startswith(">")][0][1:]

    os.chdir(output_folder)
    cluster_run_path = os.path.join(output_folder, f"cluster_run_{af_jobname}.sh")
    slurm_ids_log_path = os.path.join(output_folder, f"slurm_logs_{af_jobname}.sh")

    # check if done
    result_files = [os.path.join(output_folder, filename) for filename in os.listdir(output_folder)
                    if filename.endswith(".pdb") and filename.startswith(af_jobname + "_unrelaxed")
                    and os.path.exists(get_scores_from_path(os.path.join(output_folder, filename)))]
    if result_files and not check_for_active_slurm(slurm_ids_log_path):
        return RunAlphaFoldResult.SUCCESS, result_files

    # for phoenix
    # running_line = f"sbatch --gres=gpu:1,vmem:20g --exclude=gsm-03 --time=24:0:0 --mem=20G " \
    #                f"--killable {cluster_run_path} >> {slurm_ids_log_path}"
    # for moriah
    running_line = f"sbatch --gres=gpu:a100-1-10 --time=24:0:0 --mem=30G {cluster_run_path} >> {slurm_ids_log_path}"
    # running_line = f"sbatch --gres=gpu:a30 --time=24:0:0 --mem=30G {cluster_run_path} >> {slurm_ids_log_path}"

    if os.path.exists(cluster_run_path):
        # check if running
        if check_for_active_slurm(slurm_ids_log_path):
            return RunAlphaFoldResult.RUNNING, []

        # check if failed
        last_output = get_last_slurm_output(output_folder, slurm_ids_log_path)
        if last_output is None:
            print("Some error with cluster run, no log of canceled run, rerunning")
        elif "DUE TO PREEMPTION" in last_output or "CANCELLED AT" in last_output:
            print(f"rerunning due to preemption {os.path.basename(fasta_path)}")
        elif ("[Errno 110] Connection timed out" in last_output or "Connection reset by peer" in last_output
              or "MMseqs2 API is giving errors" in last_output):
            print(f"rerunning due to connection timeout/reset {os.path.basename(fasta_path)}")
        elif "CUDA_ERROR_ILLEGAL_ADDRESS" in last_output:
            print(f"rerunning due to CUDA error {os.path.basename(fasta_path)}")
        elif "Out of memory while trying to allocate" in last_output:
            running_line = f"sbatch --gres=gpu:a30 --time=24:0:0 --mem=30G {cluster_run_path} >> {slurm_ids_log_path}"
            if running_line in open(cluster_run_path, "r").read():
                print(f"skipping because failed {af_jobname} on memory len: {len(open(fasta_path, 'r').read())} ")
                return RunAlphaFoldResult.FAILED, []
            print(f"rerunning due to out of memory {os.path.basename(fasta_path)}")
        else:
            print(f"skipping because failed {af_jobname} len: {len(open(fasta_path, 'r').read())} ")
            return RunAlphaFoldResult.FAILED, []

    logname = "log_" + os.path.basename(fasta_path).split(".")[0] + ".txt"
    with open(cluster_run_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write('export PATH="/sci/labs/dina/bshor/projects/colabfold/phoenix/202302/colabfold_batch/bin:$PATH"\n')
        f.write("module load cuda/11.1\n")
        f.write("module load cudnn/8.0.5\n")
        f.write(f"colabfold_batch {fasta_path} {output_folder} --logname {logname} --num-models 5"
                f" --num-recycle 3\n")

    subprocess.run(running_line, shell=True)
    print("Started Alphafold", af_jobname)
    return RunAlphaFoldResult.RUNNING, []
