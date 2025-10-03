# config_loader.py
from dataclasses import dataclass
from typing import Dict, List
import yaml

@dataclass
class ToleranceConfig:
    """
    Gap/overlap related parameters.
    """
    tolerance: int = 1000              # general gap tolerance used for coarse filtering
    small_overlap: int = 300           # overlaps <= this are tolerated/ignored
    boosted_tolerance: int = 1000      # baseline boosted tolerance for asymmetry boosting
    max_boost: int = 5000              # hard cap to prevent huge boosted tolerance
    mode_gap_tolerance: int = 200      # how close to mode is considered 'around mode'
    min_large_fraction: float = 0.25   # fraction: if large_gap_count/total_count < this -> ignore singletons
    absolute_small_gap: int = 200      # absolute bp threshold below which gap is always acceptable

@dataclass
class PlottingConfig:
    outdir: str = "figures"
    prefix: str = "global"
    hist_bins: int = 30
    scatter_point_size: int = 30
    scatter_alpha: float = 0.6

@dataclass
class Config:
    patterns: Dict[str, List[str]]
    tolerances: ToleranceConfig
    plotting: PlottingConfig

def load_config(yaml_path: str) -> Config:
    """
    Load YAML configuration into Config dataclass structure.
    """
    with open(yaml_path, "r") as fh:
        data = yaml.safe_load(fh)

    tol = data.get("block_tolerances", {})
    plot = data.get("block_plotting", {})
    patterns = data.get("block_patterns", {})

    return Config(
        patterns=patterns,
        tolerances=ToleranceConfig(**tol),
        plotting=PlottingConfig(**plot)
    )
