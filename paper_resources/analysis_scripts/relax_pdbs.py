# Notice! this script relays on running in a colabfold environment. faster if using GPU. if not - mark use_gpu=False
from alphafold.common import protein, residue_constants
from alphafold.relax import relax
from pathlib import Path
import sys
import os
import time
import Bio.PDB


def relax_pdb_path(pdb_path, output_path, use_gpu=False):
    pdb_lines = Path(pdb_path).read_text()
    pdb_obj = protein.from_pdb_string(pdb_lines)

    amber_relaxer = relax.AmberRelaxation(
        max_iterations=0,
        tolerance=2.39,
        stiffness=10.0,
        exclude_residues=[],
        max_outer_iterations=3,
        use_gpu=use_gpu)

    start_time = time.time()
    print("Relaxing", pdb_path, flush=True)
    relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=pdb_obj)
    open(output_path, "w").write(relaxed_pdb_lines)
    print("Relax took", time.time() - start_time, flush=True)

    return relaxed_pdb_lines


def relax_folder(input_folder: str, output_folder: str, use_gpu=False):
    os.makedirs(output_folder, exist_ok=True)
    input_names = [i for i in os.listdir(input_folder) if i.endswith(".pdb")]
    # sort by size
    input_names.sort(key=lambda x: os.path.getsize(os.path.join(input_folder, x)))
    print("Relaxing", len(input_names), "pdbs", input_names, flush=True)

    for pdb_name in input_names:
        pdb_path = os.path.join(input_folder, pdb_name)
        output_path = os.path.join(output_folder, pdb_name[:-4] + "_relaxed.pdb")
        # if os.path.exists(output_path):
        #     print("Skipping", pdb_path, "already exists", flush=True)
        #     continue
        # relax_pdb_path(pdb_path, output_path, use_gpu)
        if os.path.exists(output_path):
            print(f"--fixing chain ids {pdb_name}", flush=True)
            fix_chain_ids(pdb_path, output_path, output_path)


def fix_chain_ids(original_pdb_path: str, relaxed_pdb_path: str, output_path: str):
    # get original chain ids
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    original_pdb_model = next(iter(pdb_parser.get_structure("o_pdb", original_pdb_path)))
    original_chain_ids = [c.id for c in original_pdb_model.get_chains()]

    relaxed_pdb_struct = pdb_parser.get_structure("r_pdb", relaxed_pdb_path)
    relaxed_pdb_model = next(iter(relaxed_pdb_struct))
    relaxed_chain_ids = [c.id for c in relaxed_pdb_model.get_chains()]

    PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    available_chain_ids = [c for c in PDB_CHAIN_IDS if (c not in relaxed_chain_ids and c not in original_chain_ids)]
    assert len(original_chain_ids) == len(relaxed_chain_ids)
    chain_id_to_chain = {c.id: c for c in relaxed_pdb_model.get_chains()}

    for original_name, new_name in zip(original_chain_ids, relaxed_chain_ids):
        if original_name in relaxed_pdb_model:
            tmp_name = available_chain_ids.pop()
            relaxed_pdb_model[original_name].id = tmp_name
            print("TempRenaming", original_name, "to", tmp_name, flush=True)
        chain_id_to_chain[new_name].id = original_name
        print("Renaming", new_name, "to", original_name, flush=True)

    io = Bio.PDB.PDBIO()
    io.set_structure(relaxed_pdb_struct)
    io.save(output_path)


if __name__ == "__main__":
    assert len(sys.argv) == 3, "Usage: python af_relax.py <input_pdb_path> <output_pdb_path>"

    path1, path2 = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    if os.path.isdir(path1):
        relax_folder(path1, path2, True)
    else:
        relax_pdb_path(path1, path2, True)
