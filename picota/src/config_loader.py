from dataclasses import dataclass
from typing import Dict, List, Optional
import yaml

# -------------------------
# Tolerance config
# -------------------------
@dataclass
class ToleranceConfig:
    tolerance: int = 1000
    small_overlap: int = 300
    boosted_tolerance: int = 1000
    max_boost: int = 5000
    mode_gap_tolerance: int = 200
    min_large_fraction: float = 0.25
    absolute_small_gap: int = 200
    max_reasonable_gap: int = 3000  

# -------------------------
# Plotting config
# -------------------------
@dataclass
class PlottingConfig:
    outdir: str = "figures"
    prefix: str = "global"
    hist_bins: int = 30
    scatter_point_size: int = 30
    scatter_alpha: float = 0.6

# -------------------------
# Logging config
# -------------------------
@dataclass
class LoggingConfig:
    log_file: str = "picota/picota.log"
    level: str = "INFO"
    logger_name: str = "picota_analysis"

# -------------------------
# Paths config
# -------------------------
@dataclass
class PathsConfig:
    outdir: str
    sra_id_file: str
    fastq_dump: str
    assembly_threads: int
    assembly_k_mer_list: str
    quiet: bool
    keep_temp_files: bool
    path_of_spades: str
    path_of_fastp: str
    skip_filtering: bool
    assembler_type: str
    path_of_megahit: str
    gfa_tools_path: str
    path_of_bandage: str
    path_to_antibiotics: str
    path_to_xenobiotics: str
    path_to_ises: str
    path_to_TNs: str
    find_all_path: bool
    path_limit: int

# -------------------------
# Options config
# -------------------------
@dataclass
class OptionsConfig:
    delete_fastq_files: bool
    min_size_of_cycle: int
    max_size_of_cycle: int
    name_prefix_cycle: str
    min_component_number: int
    max_component_number: int
    k_mer_sim: int
    threshold_sim: int
    mapping_threads: int
    mean_of_CompTns: int
    std_of_CompTns: int
    total_score_type: int
    threshold_final_score: int
    max_z: int
    dist_type: int
    path_of_prodigal: str
    path_of_blastn: str
    path_of_makeblastdb: str
    path_of_blastx: str
    path_of_blastp: str
    split_min_score: int


# -------------------------
# Main config
# -------------------------
@dataclass
class Config:
    paths: PathsConfig
    options: OptionsConfig
    patterns: Dict[str, List[str]]
    tolerances: ToleranceConfig
    plotting: PlottingConfig
    logging: LoggingConfig

# -------------------------
# Loader
# -------------------------
def load_config(yaml_path: str) -> Config:
    with open(yaml_path, "r") as fh:
        data = yaml.safe_load(fh)

    return Config(
        paths=PathsConfig(**data["paths"]),
        options=OptionsConfig(**data["options"]),
        patterns=data.get("patterns", {}),
        tolerances=ToleranceConfig(**data.get("tolerances", {})),
        plotting=PlottingConfig(**data.get("plotting", {})),
        logging=LoggingConfig(**data.get("logging", {}))
    )