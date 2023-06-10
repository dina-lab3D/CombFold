import dataclasses
import json
from enum import Enum
from typing import List, Dict, Tuple, Optional

# Constants
INTERFACE_MIN_ATOM_DIST = 8.0


SubunitName = str
ChainedSubunitName = str
ChainName = str
PdbPath = str


@dataclasses.dataclass
class SubunitInfo:
    name: SubunitName
    chain_names: List[str]
    start_res: int  # inclusive, should be greater than 0
    sequence: str

    def get_unstructured_res_ids(self) -> List[int]:
        return [self.start_res + i for i, resname in enumerate(self.sequence) if resname == "X"]

    def get_end_res(self) -> int:
        return self.start_res + len(self.sequence) - 1

    def get_active_res_ids(self) -> List[int]:
        return [i for i in range(self.start_res, self.get_end_res() + 1) if i not in self.get_unstructured_res_ids()]

    def get_relative_active_res_ids(self) -> List[int]:
        return [i - self.start_res for i in range(self.start_res, self.get_end_res() + 1)
                if i not in self.get_unstructured_res_ids()]

    def get_chained_names(self) -> List[str]:
        return [f"{self.name}_{chain_name}" for chain_name in self.chain_names]

    def to_dict(self):
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


SubunitsInfo = Dict[SubunitName, SubunitInfo]


def save_subunits_info(subunits_info: SubunitsInfo, output_path: str):
    json_data = {domain_name: domain_info.to_dict() for domain_name, domain_info in subunits_info.items()}
    json.dump(json_data, open(output_path, "w"), indent=2)


def read_subunits_info(output_path: str) -> SubunitsInfo:
    json_data = json.load(open(output_path))
    subunits_info = {domain_name: SubunitInfo.from_dict({k: v for k, v in domain_info.items() if k != "end_res"})
                     for domain_name, domain_info in json_data.items()}
    for domain_name, domain_info in subunits_info.items():
        assert domain_name == domain_info.name, "Domain name and domain info name must match"
    return subunits_info


@dataclasses.dataclass
class SubunitPdbInfo:
    chain_id: ChainName
    chain_residue_id: int
    pdb_residue_id: int
    length: int
