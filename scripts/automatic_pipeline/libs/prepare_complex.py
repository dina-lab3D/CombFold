import os
import shutil
import sys

import numpy as np
from typing import List, Optional
from math import cos, sin

import Bio
import Bio.PDB, Bio.SeqIO, Bio.SeqUtils
import Bio.PDB.Residue


def read_model_path(pdb_path: str):
    if pdb_path.endswith(".cif"):
        return Bio.PDB.MMCIFParser().get_structure("s_cif", pdb_path)
    return Bio.PDB.PDBParser(QUIET=True).get_structure("s_pdb", pdb_path)


def _merge_models(model_path1: str, model_path2: str, output_path: str, output_cif: bool = False):
    # print(f"merging {model_path1} and {model_path2} to {output_path}")
    model_struct1 = read_model_path(model_path1)
    model_struct2 = read_model_path(model_path2)

    model_1 = next(iter(model_struct1))
    model_2 = next(iter(model_struct2))

    chains_1 = {i.get_id(): i for i in model_1.get_chains()}
    chains_2 = {i.get_id(): i for i in model_2.get_chains()}

    for chain in chains_2.values():
        if chain.id in chains_1:
            res_to_move = [i for i in chain.get_residues()]
            chains_1[chain.id].child_list += res_to_move
            chains_1[chain.id].child_list.sort(key=lambda x: x.id[1])
            for res in res_to_move:
                chain.detach_child(res.id)
                res.detach_parent()
                res.parent = chains_1[chain.id]

            chains_1[chain.id].child_list += chain.child_list
        else:
            model_1.child_list.append(chain)
            chain.parent = model_1
    # save the new structure to output_pdb
    if output_cif:
        io = Bio.PDB.MMCIFIO()
    else:
        io = Bio.PDB.PDBIO()
    io.set_structure(model_struct1)
    io.save(output_path)


def _rotate_atom(coord, euler_rotation_tuple):
    # based on ** gamb::Matrix3, gamb::RigidTrans **
    x, y, z = euler_rotation_tuple
    cx, cy, cz = cos(x), cos(y), cos(z)
    sx, sy, sz = sin(x), sin(y), sin(z)
    v0 = np.array([cz*cy, -sy*sx*cz - sz*cx, -sy*cx*cz + sz*sx])
    v1 = np.array([sz*cy, -sy*sx*sz + cx*cz, -sy*cx*sz - sx*cz])
    v2 = np.array([sy, cy*sx, cy*cx])

    return np.array([
        v0[0] * coord[0] + v0[1] * coord[1] + v0[2] * coord[2],
        v1[0] * coord[0] + v1[1] * coord[1] + v1[2] * coord[2],
        v2[0] * coord[0] + v2[1] * coord[1] + v2[2] * coord[2],
        ])


def apply_transform(pdb_path: str, output_path: str, transform_numbers: List[float]):
    assert len(transform_numbers) == 6

    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = pdb_parser.get_structure("original_pdb", pdb_path)
    assert len(list(pdb_struct)) == 1, f"Too many models! {pdb_path}"

    rotation_tuple, translation_tuple = transform_numbers[:3], transform_numbers[3:]
    for atom in next(iter(pdb_struct)).get_atoms():
        coord = atom.get_coord()
        new_coord = _rotate_atom(coord, rotation_tuple)
        new_coord += np.array(translation_tuple)
        atom.set_coord(new_coord)

    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def create_transformation_pdb(assembly_path: str, transforms_str: str, output_path: str, output_cif: bool = False):
    # example: transforms_str = "0(-0.00760431 -1.00941 -2.38469 21.6924 -87.7727 -10.4274),1(0 0 0 0 0 0),
    # 2(-2.96419 0.240678 -0.964002 65.6247 -63.7248 -42.0541)"
    tmp_pdb_path = os.path.join(assembly_path, "tmp_comb_trans.pdb")
    chains_path = os.path.join(assembly_path, "chain.list")
    if not os.path.exists(chains_path):
        chains_path = os.path.join(assembly_path, "chains.txt")
    subunit_filenames = [i.split(" ")[0] for i in open(chains_path, "r").read().split("\n") if len(i) > 0]

    if os.path.exists(output_path):
        os.remove(output_path)
    for transform_str in transforms_str.split(","):
        chain_ind = int(transform_str.split("(")[0])
        subunit_path = os.path.join(assembly_path, subunit_filenames[chain_ind])
        transform_numbers_as_str = transform_str.split("(")[1][:-1]
        transform_numbers = list(map(float, transform_numbers_as_str.split(" ")))

        apply_transform(subunit_path, tmp_pdb_path, transform_numbers)
        if os.path.exists(output_path):
            _merge_models(output_path, tmp_pdb_path, output_path, output_cif=output_cif)
        else:
            if output_cif:
                io = Bio.PDB.MMCIFIO()
                io.set_structure(read_model_path(tmp_pdb_path))
                io.save(output_path)
            else:
                shutil.copy(tmp_pdb_path, output_path)
    if output_cif:
        # fix _atom_site in cif
        cif_lines = open(output_path, "r").read().split("\n")
        cif_lines = [i + " " if i.startswith("_atom_site") else i for i in cif_lines]
        open(output_path, "w").write("\n".join(cif_lines) + "\n")
    os.remove(tmp_pdb_path)


def create_complexes(result_path: str, first_result: Optional[int] = None, last_result: Optional[int] = None,
                     output_folder: Optional[str] = None, output_cif: bool = False) \
        -> List[str]:
    output_folder = output_folder or os.path.dirname(result_path)
    os.makedirs(output_folder, exist_ok=True)
    assembly_path = os.path.dirname(result_path)
    transforms_strs = [i[i.index("[") + 1:i.index("]")] for i in open(result_path, "r").read().split("\n") if i]

    if first_result is None:
        first_result = 0
    if last_result is None or last_result > len(transforms_strs):
        last_result = len(transforms_strs)

    output_files = []
    for i in range(first_result, last_result):
        file_ext = ".cif" if output_cif else ".pdb"
        output_path = os.path.join(output_folder,
                                   os.path.basename(result_path).split(".")[0] + "_" + str(i) + file_ext)
        create_transformation_pdb(assembly_path, transforms_strs[i], output_path=output_path, output_cif=output_cif)
        output_files.append(output_path)
    return output_files


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise Exception("Usage: <script> result_path first_result_num last_result_num")
    _result_path = os.path.abspath(sys.argv[1])
    _first_result = int(sys.argv[2]) - 1
    _last_result = int(sys.argv[3])

    create_complexes(_result_path, _first_result, _last_result)
