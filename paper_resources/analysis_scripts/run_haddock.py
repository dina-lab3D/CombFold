import os

"""
This assumes you already in haddock environment
bash
conda activate haddock2.4
source {PROJECT_PATH}/tools/haddock/haddock2.4-2023-08/haddock_configure.sh
"""

# TODO: define project path
PROJECT_PATH = ""
HADDOCK_RUN = f"{PROJECT_PATH}/tools/conda_install/miniconda3/envs/haddock2.4/bin/python2.7" \
              f" {PROJECT_PATH}/tools/haddock/haddock2.4-2023-08/Haddock/RunHaddock.py"


def run_haddock(assembly_output_folder: str, haddock_output_folder: str):
    chain_filenames = [i for i in open(os.path.join(assembly_output_folder, "chain.list"), "r").read().split("\n") if i]
    if len(chain_filenames) > 19:
        print("ERROR: too many chains")
        return

    os.makedirs(haddock_output_folder, exist_ok=True)
    os.chdir(haddock_output_folder)

    # copy chain files
    for chain_filename in chain_filenames:
        os.system(f"cp {os.path.join(assembly_output_folder, chain_filename)} {haddock_output_folder}")

    # create xlinks
    original_xlinks = open(os.path.join(assembly_output_folder, "xlink_consts.txt"), "r").read().split("\n")
    new_xlinks = ""
    for line in original_xlinks:
        if not line:
            continue
        res1, chains1, res2, chains2, min_dist, max_dist, score = line.split(" ")
        # new_xlinks += f"{chain1} {res1} CA {chain2} {res2} CA 0.0 30.0\n"
        for chain1 in chains1:
            for chain2 in chains2:
                new_xlinks += f"assign (name ca and segid {chain1} and resi {res1}) " \
                              f"(name ca and segid {chain2} and resi {res2}) 20.0 15.0 15.0\n"
    open(os.path.join(haddock_output_folder, "xlinks.tbl"), "w").write(new_xlinks)


    # create params file
    pararm_content = f"HADDOCK_DIR={PROJECT_PATH}/tools/haddock/haddock2.4-2023-08\n"
    pararm_content += f"N_COMP={len(chain_filenames)}\n"
    for i, chain_filename in enumerate(chain_filenames):
        pararm_content += f"PDB_FILE{i+1}={chain_filename}\n"
    pararm_content += "PROJECT_DIR=./\n"
    pararm_content += "RUN_NUMBER=1\n"
    pararm_content += "UNAMBIG_TBL=./xlinks.tbl\n"

    open(os.path.join(haddock_output_folder, "run.param"), "w").write(pararm_content)

    os.system(HADDOCK_RUN)

    # optional = reduce number of structures
    os.chdir("run1")
    run_cns_content = open("run.cns", "r").read()
    run_cns_content = run_cns_content.replace("structures_0=1000", "structures_0=50")
    run_cns_content = run_cns_content.replace("structures_1=200", "structures_1=10")
    run_cns_content = run_cns_content.replace("anastruc_1=200", "anastruc_1=10")

    open("run.cns", "w").write(run_cns_content)

    with open("haddock_tmp_run.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(HADDOCK_RUN)
    os.system("sbatch --time=24:00:00 --mem=16G haddock_tmp_run.sh")

    os.chdir("../..")


def main():
    for pdb_id in os.listdir(f"{PROJECT_PATH}/afm3/base_xlinks25"):
        print("running", pdb_id, "cwd", os.getcwd())
        input("continue?")
        if pdb_id in ("7QRU", "7T3B"):
            continue
        run_haddock(f"{PROJECT_PATH}/afm3/base_xlinks25/{pdb_id}/assembly_output",
                    f"{PROJECT_PATH}/afm3/haddock/{pdb_id}")


if __name__ == '__main__':
    main()
