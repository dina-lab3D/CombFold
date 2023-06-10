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
    # return {domain_name: SubunitInfo.from_dict(domain_info) for domain_name, domain_info in json_data.items()}
    return {domain_name: SubunitInfo.from_dict({k: v for k, v in domain_info.items() if k != "end_res"})
            for domain_name, domain_info in json_data.items()}


@dataclasses.dataclass
class AlphaFoldJobInfo:
    subunit_names: List[SubunitName]
    merged_subunits: List[bool]  # if True, then the subunit is merged with the next one
    sequences: List[str]

    def get_jobname(self):
        return "_".join(self.subunit_names)

    def get_as_fasta(self):
        return f">{self.get_jobname()}\n" + ":".join(self.sequences) + "\n"

    def __hash__(self):
        return hash(tuple([*self.subunit_names, *self.merged_subunits]))


class RunAlphaFoldResult(Enum):
    SUCCESS = 0
    NOT_STARTED = 1
    RUNNING = 2
    FAILED = 3
    SHOULD_RERUN = 4


@dataclasses.dataclass
class AFSubunitScores:
    plddt_avg: float
    plddt_percentile: List[float]  # list of 11 numbers, 0-100 percentile in skips of 10
    plddt_interface_avg: float
    plddt_interface_percentile: List[float]  # list of 11 numbers, 0-100 percentile in skips of 10

    self_pae_avg: float
    self_pae_percentile: List[float]  # list of 11 numbers, 0-100 percentile in skips of 10


@dataclasses.dataclass
class AFInteractionScores:
    pae_avg: float
    pae_percentile: List[float]  # list of 11 numbers, 0-100 percentile in skips of 10
    pae_joined_interface_avg: float
    pae_joined_interface_percentile: List[float]  # list of 11 numbers, 0-100 percentile in skips of 10

    interface1_size: int
    interface2_size: int


@dataclasses.dataclass
class SubunitPdbInfo:
    chain_id: ChainName
    chain_residue_id: int
    pdb_residue_id: int
    length: int


@dataclasses.dataclass
class AFResultScoredPair:
    subunits_names: Tuple[SubunitName, SubunitName]
    pdb_path: str

    subunit1_pdb_info: SubunitPdbInfo
    subunit2_pdb_info: SubunitPdbInfo

    # chains_in_pdb: Tuple[ChainName, ChainName]

    subunit1_scores: Optional[AFSubunitScores] = None
    subunit2_scores: Optional[AFSubunitScores] = None
    interaction_scores: Optional[AFInteractionScores] = None
