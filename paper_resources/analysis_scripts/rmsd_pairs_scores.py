import json
import os
import pickle
import sys
import Bio.PDB
import Bio.SeqUtils


BENCHMARK_PATH = "runs/20230330_afm3_benchmark_rec3"
INPUT_COMPLEXES_PATH = "afm3/input_complexes/"

OUTPUT_PATH = os.path.join(BENCHMARK_PATH, "rmsd_validate.json")
sys.path.append(BENCHMARK_PATH)


def get_scores_from_path(af_pdb_path: str):
    return af_pdb_path.replace("unrelaxed", "scores").replace("pdb", "json")


def get_pdb_model_readonly(pdb_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    return next(iter(pdb_struct))


def get_rmsd(input_model, chain_id1, chain_id2, sample_model):
    chain1 = input_model[chain_id1]
    chain2 = input_model[chain_id2]

    ref_atoms = []
    ref_ids = []
    sample_atoms = []

    for res in list(chain1.get_residues()) + list(chain2.get_residues()):
        if "CA" not in res or Bio.SeqUtils.seq1(res.get_resname()) == "X":
            continue
        ref_atoms.append(res["CA"])
        ref_ids.append((res.parent.id, res.get_id()[1]))

    # print("ref_ids", ref_ids)

    for res in sample_model["A"].get_residues():
        if (chain_id1, res.get_id()[1]) not in ref_ids:
            continue
        sample_atoms.append(res["CA"])

    for res in sample_model["B"].get_residues():
        if (chain_id2, res.get_id()[1]) not in ref_ids:
            continue
        sample_atoms.append(res["CA"])

    assert len(ref_atoms) == len(sample_atoms), f"len(ref_atoms)={len(ref_atoms)} len(sample_atoms)={len(sample_atoms)}"

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    return super_imposer.rms


def main():
    scores_by_pdb = {}
    for pdb_id in os.listdir(os.path.join(BENCHMARK_PATH, "output")):
        # only those that are completely asymetric, to avoid confusion of multiple interfaces
        if pdb_id not in ["7ZKQ", "8A5O", "7USE", "7PKN", "7T3B", "8HIL", "7ARC", "7UIC", "8E9G", "7WFF", "7OBA",
                          "8F5O"]:
            continue

        print("Doing", pdb_id)
        pairs_pickle_path = os.path.join(BENCHMARK_PATH, "output", pdb_id, "parsed_pairs.pkl")
        if not os.path.exists(pairs_pickle_path):
            continue
        with open(pairs_pickle_path, "rb") as f:
            all_pairs = pickle.load(f)
        merged_scores = []

        input_complex_path = os.path.join(INPUT_COMPLEXES_PATH, f"{pdb_id}.pdb")
        input_model = get_pdb_model_readonly(input_complex_path)

        for result in all_pairs:
            json_scores = json.load(open(get_scores_from_path(result.pdb_path), "r"))
            if "iptm" not in json_scores or json_scores["iptm"] is None:
                print(f"{os.path.basename(result.pdb_path)} no iptm")
                continue
            iptm = json_scores["iptm"]

            sample_model = get_pdb_model_readonly(result.pdb_path)
            chain1, chain2 = result.subunits_names[0][0], result.subunits_names[1][0]
            if result.subunits_names[0][2] != "F" or result.subunits_names[1][2] != "F":
                print(f"{os.path.basename(result.pdb_path)} not full")
                continue


            rmsd = get_rmsd(input_model, chain1, chain2, sample_model)

            merged_scores.append(
                {"pdb": os.path.basename(result.pdb_path),
                 "iptm": iptm,
                 "rmsd": rmsd,
                 "s1_plddt": result.subunit1_scores.plddt_avg,
                 "s1_iplddt": result.subunit1_scores.plddt_interface_avg,
                 "s1_pae": result.subunit1_scores.self_pae_avg,
                 "s2_plddt": result.subunit2_scores.plddt_avg,
                 "s2_iplddt": result.subunit2_scores.plddt_interface_avg,
                 "s2_pae": result.subunit2_scores.self_pae_avg,
                 "full_pae": result.interaction_scores.pae_avg,
                 "interface_pae": result.interaction_scores.pae_joined_interface_avg,
                 "interface1_size": result.interaction_scores.interface1_size,
                 "interface2_size": result.interaction_scores.interface2_size,
                 }
            )
            print(merged_scores[-1])
        scores_by_pdb[pdb_id] = merged_scores
        json.dump(scores_by_pdb, open(OUTPUT_PATH, "w"))


if __name__ == '__main__':
    main()