#!/usr/bin/env python3
"""Create the South Tyrol snow bulletin and SCD products."""

from __future__ import annotations

import argparse
import time
from datetime import datetime as dt
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_AOI = PROJECT_ROOT / "aux" / "SouthTyrol.geojson"
DEFAULT_WORK_DIR = Path("/mnt/CEPH_PRODUCTS/EURAC_SNOW/ELABORATIONS")
RASTER_PATTERN = "EURAC_SNOW_MERGE.alps.south-tyrol.20*.tif"

PRODUCT_CONFIG = {
    "modis": {
        "suffix": "modis",
        "csv_name": "ST_modis.csv",
        "input_dir": Path("/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/ST"),
    },
    "viirs": {
        "suffix": "vnp10a1f",
        "csv_name": "ST_viirs.csv",
        "input_dir": Path(
            "/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol"
        ),
    },
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "product",
        choices=PRODUCT_CONFIG,
        help="Snow product to process",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        help="Override the product's default raster directory",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=DEFAULT_WORK_DIR,
        help=f"Output root (default: {DEFAULT_WORK_DIR})",
    )
    parser.add_argument(
        "--aoi",
        type=Path,
        default=DEFAULT_AOI,
        help=f"Area of interest (default: {DEFAULT_AOI})",
    )
    return parser


def discover_rasters(input_dir: Path) -> list[Path]:
    """Find product rasters below an input directory in chronological order."""
    return sorted(input_dir.rglob(RASTER_PATTERN))


def run(args: argparse.Namespace) -> int:
    # Keep the geospatial dependencies out of argument parsing so --help remains
    # available even when the processing environment is not loaded.
    if __package__:
        from .utils import (
            get_scd_statistics,
            monthly_anomaly_scd,
            snow_bullettin,
            updateCSV,
        )
    else:
        from utils import (  # type: ignore[no-redef]
            get_scd_statistics,
            monthly_anomaly_scd,
            snow_bullettin,
            updateCSV,
        )

    config = PRODUCT_CONFIG[args.product]
    suffix = str(config["suffix"])
    input_dir = args.input_dir or config["input_dir"]

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if not args.aoi.is_file():
        raise FileNotFoundError(f"AOI file does not exist: {args.aoi}")

    snow_maps = discover_rasters(input_dir)
    if not snow_maps:
        raise FileNotFoundError(
            f"No rasters matching {RASTER_PATTERN!r} found below {input_dir}"
        )

    bulletin_dir = args.work_dir / "BULLETIN"
    bulletin_dir.mkdir(parents=True, exist_ok=True)
    csv_path = bulletin_dir / str(config["csv_name"])

    print(f"Processing {args.product}: {len(snow_maps)} raster(s) from {input_dir}")
    start_time = time.time()
    sca_data = updateCSV(
        str(csv_path),
        [str(path) for path in snow_maps],
        shp_fileName=args.aoi,
    )
    print(f"SCA statistics updated in {time.time() - start_time:.1f} seconds")

    if sca_data is None or sca_data.empty:
        raise RuntimeError(f"No SCA data are available in {csv_path}")

    date_end_value = sca_data.index[-1]
    if date_end_value.month >= 10:
        date_start_value = dt(date_end_value.year, 10, 1)
    else:
        date_start_value = dt(date_end_value.year - 1, 10, 1)

    date_start = date_start_value.strftime("%Y-%m-%d")
    date_end = date_end_value.strftime("%Y-%m-%d")
    print(f"Bulletin period: {date_start} through {date_end}")
    snow_bullettin(
        str(csv_path),
        date_start,
        date_end,
        str(bulletin_dir),
        suffix,
    )

    scd_root = args.work_dir / "SCD" / suffix
    snow_map_names = [str(path) for path in snow_maps]
    get_scd_statistics(
        snow_map_names,
        str(scd_root / "trimester"),
        max_missing_days=20,
        shp_fileName=args.aoi,
        window=2,
        mode="trimester",
    )
    get_scd_statistics(
        snow_map_names,
        str(scd_root / "yearly"),
        max_missing_days=71,
        shp_fileName=args.aoi,
        window=2,
        mode="yearly",
    )
    scd_monthly = get_scd_statistics(
        snow_map_names,
        str(scd_root / "monthly"),
        max_missing_days=15,
        shp_fileName=args.aoi,
        window=2,
        mode="monthly",
    )

    if scd_monthly:
        monthly_anomaly_scd(scd_monthly, str(bulletin_dir), suffix=suffix)
    else:
        print("No monthly SCD results are available; skipping the anomaly plot")

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
