from functools import lru_cache
from typing import Tuple, List

import Bio.PDB
import numpy as np
import scipy.spatial

from automatic_pipeline.libs.utils_classes import INTERFACE_MIN_ATOM_DIST, SubunitPdbInfo


@lru_cache(5)
def get_pdb_model_readonly(pdb_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, f"Too many models! {pdb_path}"
    return next(iter(pdb_struct))


def get_pdb_model_no_cache(pdb_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, f"Too many models! {pdb_path}"
    return next(iter(pdb_struct))


def extract_chains(pdb_path: str, output_path: str, chains_to_save: List[str]):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, "can't extract if more than one model"
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())
    for chain_id in chains_to_save:
        assert len([c for c in chains if c.id == chain_id]) == 1, f"Missing chain {chain_id} from {pdb_path}"
    for chain in chains:
        if chain.id not in chains_to_save:
            model.detach_child(chain.id)
    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def extract_pdb_info(pdb_path: str, pdb_info: SubunitPdbInfo, output_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, "can't extract if more than one model"
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())
    assert len([c for c in chains if c.id == pdb_info.chain_id]) == 1, f"Missing chain {pdb_info.chain_id}"
    for chain in chains:
        if chain.id != pdb_info.chain_id:
            model.detach_child(chain.id)

    res_to_remove = []
    for res in model.get_residues():
        if not (pdb_info.chain_residue_id <= res.id[1] < pdb_info.chain_residue_id + pdb_info.length):
            res_to_remove.append(res)
    for res in res_to_remove:
        res.parent.detach_child(res.id)

    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def copy_pdb_rename_chain(original_pdb_path: str, new_chain_name: str, output_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", original_pdb_path)
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())
    assert len(chains) == 1
    assert len(new_chain_name) == 1, f"Bad chain name provided {new_chain_name}"
    chains[0].id = new_chain_name

    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def copy_pdb_add_offset(original_pdb_path: str, offset: int, output_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", original_pdb_path)
    model = next(iter(pdb_struct))
    chains = list(model.get_chains())
    assert len(chains) == 1

    # prevent "already used for a sibling of this entity" error
    for res in chains[0].get_residues():
        res.id = (res.id[0], res.id[1] + 10**4, res.id[2])

    for res in chains[0].get_residues():
        res.id = (res.id[0], res.id[1] - 10**4 + offset, res.id[2])

    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def copy_pdb_set_start_offset(original_pdb_path: str, offset: int, output_path: str):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", original_pdb_path)
    model = next(iter(pdb_struct))
    current_offset = min([res.id[1] for res in model.get_residues()])
    copy_pdb_add_offset(original_pdb_path, offset - current_offset, output_path)


def is_chains_coords_close(chain1_ca: np.ndarray, chain2_ca: np.ndarray):
    if len(chain1_ca) == 0 or len(chain2_ca) == 0:
        return False
    return np.min(scipy.spatial.distance.cdist(chain1_ca, chain2_ca)) < INTERFACE_MIN_ATOM_DIST


def are_subunits_close_in_pdb(pdb_path: str, subunit1: SubunitPdbInfo, subunit2: SubunitPdbInfo):
    model = get_pdb_model_readonly(pdb_path)
    chains = list(model.get_chains())

    chain1 = [i for i in chains if i.get_id() == subunit1.chain_id][0]
    chain2 = [i for i in chains if i.get_id() == subunit2.chain_id][0]

    chain1_ca = np.array([res['CA'].coord for res in chain1.get_residues() if 'CA' in res
                          and subunit1.chain_residue_id <= res.get_id()[1]
                          <= subunit1.chain_residue_id + subunit1.length])
    chain2_ca = np.array([res['CA'].coord for res in chain2.get_residues() if 'CA' in res
                          and subunit2.chain_residue_id <= res.get_id()[1]
                          <= subunit2.chain_residue_id + subunit2.length])

    assert len(chain1_ca) > 0 and len(chain2_ca) > 0, f"No CA atoms found {subunit1} {subunit2} {pdb_path}"

    return is_chains_coords_close(chain1_ca, chain2_ca)


def get_interface_res_ids(pdb_path: str, subunit1: SubunitPdbInfo, subunit2: SubunitPdbInfo) \
        -> Tuple[List[int], List[int]]:
    model = get_pdb_model_readonly(pdb_path)
    chains = list(model.get_chains())

    chain1 = [i for i in chains if i.get_id() == subunit1.chain_id][0]
    chain2 = [i for i in chains if i.get_id() == subunit2.chain_id][0]

    chain1_res = [res for res in chain1.get_residues() if 'CA' in res and subunit1.chain_residue_id <= res.get_id()[1]
                  <= subunit1.chain_residue_id + subunit1.length]
    chain2_res = [res for res in chain2.get_residues() if 'CA' in res and subunit2.chain_residue_id <= res.get_id()[1]
                  <= subunit2.chain_residue_id + subunit2.length]

    chain1_ca = np.array([res['CA'].coord for res in chain1_res])
    chain2_ca = np.array([res['CA'].coord for res in chain2_res])

    close_residues = np.argwhere(scipy.spatial.distance.cdist(chain1_ca, chain2_ca) < INTERFACE_MIN_ATOM_DIST)

    chain1_interface, chain2_interface = set(), set()
    for i, j in close_residues:
        chain1_interface.add(chain1_res[i].id[1] - subunit1.chain_residue_id)
        chain2_interface.add(chain2_res[j].id[1] - subunit2.chain_residue_id)

    return sorted(list(chain1_interface)), sorted(list(chain2_interface))
