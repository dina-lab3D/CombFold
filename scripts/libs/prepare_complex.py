import os
import shutil
import sys

import numpy as np
from typing import List, Optional
from math import cos, sin

import Bio
import Bio.PDB, Bio.SeqIO, Bio.SeqUtils
import Bio.PDB.Residue


def _merge_pdbs(pdb1_path: str, pdb2_path: str, output_path: str, override_chain_names: bool = False,
                override_res_nums: bool = False):
    pdb1_lines = open(pdb1_path, "rb").read().split(b"\n")
    pdb2_lines = open(pdb2_path, "rb").read().split(b"\n")

    final_lines = [line for line in pdb1_lines if line[:6] == b"ATOM  "]

    pdb1_chains = set()
    pdb1_last_res = pdb1_last_atom = 0
    for line in pdb1_lines:
        if line[:6] != b"ATOM  ":
            continue
        pdb1_last_atom = max(pdb1_last_atom, int(line[6:11]))
        pdb1_last_res = max(pdb1_last_res, int(line[22:26]))
        pdb1_chains.add(line[21:22])
    pdb2_new_chains = {}
    pdb2_new_res = {}
    for line in pdb2_lines:
        if line[:6] != b"ATOM  ":
            continue
        new_atom_num = str((pdb1_last_atom + 1) % 100000).rjust(5).encode("ascii")
        pdb1_last_atom += 1
        if line[22:26] not in pdb2_new_res:
            if override_res_nums:
                pdb2_new_res[line[22:26]] = str(pdb1_last_res + 1).rjust(4).encode("ascii")
            else:
                pdb2_new_res[line[22:26]] = line[22:26]
            pdb1_last_res += 1
        new_res_num = pdb2_new_res[line[22:26]]

        if override_chain_names and line[21:22] not in pdb2_new_chains:
            available_chars = [bytes([c]) for c in b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                               if bytes([c]) not in pdb1_chains and bytes([c]) not in pdb2_new_chains.values()]
            pdb2_new_chains[line[21:22]] = available_chars[0]
            print(line[21:22], available_chars, pdb2_new_chains)

        final_lines.append(line[:6] + new_atom_num + line[11:21] + pdb2_new_chains.get(line[21:22], line[21:22])
                           + new_res_num + line[26:])
    open(output_path, "wb").write(b"\n".join(final_lines))


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


def fix_residue_order(pdb_path: str, output_path: str):
    parser = Bio.PDB.PDBParser(QUIET=True)
    pdb_struct = parser.get_structure("original_pdb", pdb_path)

    for chain in pdb_struct.get_chains():
        residues = [r for r in chain]
        chain.child_list = sorted(residues, key=lambda x: x.full_id)

    io = Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(output_path)


def create_transformation_pdb(assembly_path: str, transforms_str: str, output_path: str):
    # based on /cs/labs/dina/bshor/repos/PenCombDock/bin/prepareComplex.pl chain.list clusters.res 1 10 combdock_complex
    # transforms_str = "0(-0.00760431 -1.00941 -2.38469 21.6924 -87.7727 -10.4274),1(0 0 0 0 0 0),
    # 2(-2.96419 0.240678 -0.964002 65.6247 -63.7248 -42.0541)"
    tmp_pdb_path = os.path.join(assembly_path, "bshor_tmp_comb_trans.pdb")
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
            _merge_pdbs(output_path, tmp_pdb_path, output_path)
        else:
            shutil.copy(tmp_pdb_path, output_path)
    fix_residue_order(output_path, output_path)

    os.remove(tmp_pdb_path)


def create_complexes(result_path: str, first_result: Optional[int] = None, last_result: Optional[int] = None) \
        -> List[str]:
    assembly_path = os.path.dirname(result_path)
    transforms_strs = [i[i.index("[") + 1:i.index("]")] for i in open(result_path, "r").read().split("\n") if i]

    if first_result is None:
        first_result = 0
    if last_result is None:
        last_result = len(transforms_strs)

    pdb_files = []
    for i in range(first_result, min(len(transforms_strs), last_result)):
        output_path = os.path.join(assembly_path,
                                   os.path.basename(result_path).split(".")[0] + "_" + str(i) + ".pdb")
        create_transformation_pdb(assembly_path, transforms_strs[i], output_path=output_path)
        pdb_files.append(output_path)
    return pdb_files


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise Exception("Usage: <script> result_path first_result_num last_result_num")
    _result_path = os.path.abspath(sys.argv[1])
    _first_result = int(sys.argv[2]) - 1
    _last_result = int(sys.argv[3])

    create_complexes(_result_path, _first_result, _last_result)
