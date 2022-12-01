"""Unit tests for fileutils.py."""
import shutil

from pathlib import Path

from askcell.utils.fileutils import load_toml, save_toml


testdata_dir = Path("tests/data/test_utils")


def test_load_toml() -> None:
    """Test load_toml."""
    result = load_toml(testdata_dir / "qc_metrics_thresholds.toml")
    expected = {
        "swift_sid": {
            "run": {
                "total_yield": {"threshold": 10},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 90},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 15},
                    "pct_target_bases_20x": {"threshold": 90},
                },
            },
        },
        "swift_lung": {
            "run": {
                "total_yield": {"threshold": 30},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 75},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 25},
                    "pct_target_bases_20x": {"threshold": 75},
                },
            },
        },
        "chp2": {
            "run": {
                "total_yield": {"threshold": 10},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 95},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 10},
                    "pct_target_bases_20x": {"threshold": 95},
                },
            },
        },
    }
    assert result == expected


def test_save_toml() -> None:
    """Test save_toml."""
    out_dir = testdata_dir / "test_dir_out"
    out_dir.mkdir(exist_ok=True)
    to_be_saved = {
        "swift_sid": {
            "run": {
                "total_yield": {"threshold": 10},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 90},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 15},
                    "pct_target_bases_20x": {"threshold": 90},
                },
            },
        },
        "swift_lung": {
            "run": {
                "total_yield": {"threshold": 30},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 75},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 25},
                    "pct_target_bases_20x": {"threshold": 75},
                },
            },
        },
        "chp2": {
            "run": {
                "total_yield": {"threshold": 10},
                "pct_bases_gt_q30": {"threshold": 80},
                "frac_clusters_pf": {"threshold": 0.8},
                "cluster_density": {"high": {"threshold": 270}, "low": {"threshold": 120}},
            },
            "sample": {
                "hard": {"n_frags": {"threshold": 300}},
                "soft": {
                    "median_target_cov": {"threshold": 500},
                    "pct_on_target_bases": {"threshold": 95},
                    "q90_q10_ratio": {"threshold": 50},
                    "fold_80_penalty": {"threshold": 10},
                    "pct_target_bases_20x": {"threshold": 95},
                },
            },
        },
    }
    save_toml(out_dir / "tmp.toml", data=to_be_saved)
    assert (testdata_dir / "test_dir_out" / "tmp.toml").exists()
    shutil.rmtree(out_dir)
