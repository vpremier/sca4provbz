#!/usr/bin/env python3
"""Download VNP10A1F granules and create South Tyrol SCF rasters."""

from __future__ import annotations

import argparse
import logging
import re
import sys
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Callable, Iterable, Sequence

import earthaccess
import geopandas as gpd
import numpy as np

LOGGER = logging.getLogger(__name__)
PRODUCT = "VNP10A1F"
VERSION = "2"
DATE_PATTERN = re.compile(r"^VNP10A1F\.A(?P<date>\d{7})(?:\.|$)")
OUTPUT_DATE_PATTERN = re.compile(
    r"^EURAC_SNOW_MERGE\.alps\.south-tyrol\.(?P<date>\d{8})T120000\.vnp10a1f\.tif$"
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_AOI = PROJECT_ROOT / "aux" / "SouthTyrol.geojson"
DEFAULT_WATER_MASK = PROJECT_ROOT / "aux" / "Water_Mask_aligned.tif"
DEFAULT_DOWNLOAD_DIR = Path("/mnt/CEPH_PROJECTS/PROSNOW/raw_data/VIIRS/VNP10A1F")
DEFAULT_OUTPUT_DIR = Path(
    "/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol"
)
DEFAULT_EXTENT = (597972.0, 5117987.0, 767972.0, 5221987.0)
DEFAULT_ARCHIVE_START_DATE = date(2012, 1, 19)


def granule_date(path: Path) -> date:
    """Extract the acquisition date from a VNP10A1F filename."""
    match = DATE_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Not a recognized {PRODUCT} filename: {path.name}")
    return datetime.strptime(match.group("date"), "%Y%j").date()


def discover_granules(directory: Path) -> list[Path]:
    """Return valid local product granules ordered by date and filename."""
    granules = []
    for path in directory.glob(f"{PRODUCT}.A*.h5"):
        try:
            granule_date(path)
        except ValueError as error:
            LOGGER.warning("Ignoring %s", error)
        else:
            granules.append(path)
    return sorted(granules, key=lambda path: (granule_date(path), path.name))


def discover_processed_dates(directory: Path) -> set[date]:
    """Extract acquisition dates from existing SCF output filenames."""
    processed = set()
    for path in directory.glob("EURAC_SNOW_MERGE*.vnp10a1f.tif"):
        match = OUTPUT_DATE_PATTERN.match(path.name)
        if match is None:
            LOGGER.warning("Ignoring unrecognized output filename: %s", path.name)
            continue
        try:
            processed.add(datetime.strptime(match.group("date"), "%Y%m%d").date())
        except ValueError:
            LOGGER.warning("Ignoring output filename with an invalid date: %s", path.name)
    return processed


def dates_between(start: date, end: date) -> set[date]:
    """Return every calendar date in an inclusive interval."""
    return {start + timedelta(days=offset) for offset in range((end - start).days + 1)}


def contiguous_ranges(days: Iterable[date]) -> list[tuple[date, date]]:
    """Collapse dates into inclusive consecutive ranges."""
    ordered = sorted(set(days))
    if not ordered:
        return []

    ranges = []
    range_start = previous = ordered[0]
    for current in ordered[1:]:
        if current != previous + timedelta(days=1):
            ranges.append((range_start, previous))
            range_start = current
        previous = current
    ranges.append((range_start, previous))
    return ranges


def describe_dates(days: Iterable[date]) -> str:
    """Format date gaps compactly for log messages."""
    descriptions = []
    for start, end in contiguous_ranges(days):
        descriptions.append(str(start) if start == end else f"{start}..{end}")
    return ", ".join(descriptions)


def parse_iso_date(value: str) -> date:
    """Parse an ISO date for argparse."""
    try:
        return date.fromisoformat(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"Expected a date in YYYY-MM-DD format, got {value!r}"
        ) from error


def validate_extent(
    values: Sequence[float],
) -> tuple[float, float, float, float]:
    """Validate and normalize an extent as xmin, ymin, xmax, ymax."""
    xmin, ymin, xmax, ymax = (float(value) for value in values)
    if not np.isfinite((xmin, ymin, xmax, ymax)).all():
        raise ValueError("Target extent coordinates must be finite numbers")
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("Target extent must satisfy xmin < xmax and ymin < ymax")
    return xmin, ymin, xmax, ymax


def read_wgs84_bounds(aoi_path: Path) -> tuple[float, float, float, float]:
    """Read an AOI and return its WGS84 bounding box."""
    if not aoi_path.is_file():
        raise FileNotFoundError(f"AOI file does not exist: {aoi_path}")

    aoi = gpd.read_file(aoi_path)
    if aoi.empty:
        raise ValueError(f"AOI contains no features: {aoi_path}")
    if aoi.crs is None:
        raise ValueError(f"AOI has no coordinate reference system: {aoi_path}")

    aoi = aoi.to_crs("EPSG:4326")
    bounds = tuple(float(value) for value in aoi.total_bounds)
    if not np.isfinite(bounds).all():
        raise ValueError(f"AOI has invalid bounds: {bounds}")
    return bounds


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aoi", type=Path, default=DEFAULT_AOI)
    parser.add_argument("--water-mask", type=Path, default=DEFAULT_WATER_MASK)
    parser.add_argument("--download-dir", type=Path, default=DEFAULT_DOWNLOAD_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--archive-start-date",
        "--start-date",
        dest="archive_start_date",
        type=parse_iso_date,
        default=DEFAULT_ARCHIVE_START_DATE,
        help=(
            "First date to reconcile (YYYY-MM-DD; default: 2012-01-19); "
            "--start-date is an alias"
        ),
    )
    parser.add_argument(
        "--end-date", type=parse_iso_date, default=date.today(), help="Default: today"
    )
    parser.add_argument("--resolution", type=int, default=250)
    parser.add_argument("--target-epsg", default="32632")
    parser.add_argument(
        "--extent",
        type=float,
        nargs=4,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        default=DEFAULT_EXTENT,
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing SCF rasters"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=-1,
        help="Maximum search results; -1 requests all matching granules",
    )
    return parser


def run(args: argparse.Namespace, scf_processor: Callable[..., object]) -> int:
    if args.resolution <= 0:
        raise ValueError("Resolution must be greater than zero")
    if args.max_results == 0 or args.max_results < -1:
        raise ValueError("max-results must be -1 or a positive integer")

    extent = validate_extent(args.extent)
    args.download_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    start_date = args.archive_start_date
    if start_date > args.end_date:
        raise ValueError(f"Start date {start_date} is later than end date {args.end_date}")

    expected_dates = dates_between(start_date, args.end_date)
    processed_dates = discover_processed_dates(args.output_dir)
    dates_to_process = expected_dates if args.overwrite else expected_dates - processed_dates

    if not dates_to_process:
        LOGGER.info("All dates in the requested archive range are already processed")
        return 0

    LOGGER.info(
        "%d of %d date(s) need processing",
        len(dates_to_process),
        len(expected_dates),
    )

    existing = discover_granules(args.download_dir)
    local_dates = {granule_date(path) for path in existing}
    missing_raw_dates = dates_to_process - local_dates

    if missing_raw_dates:
        bounds = read_wgs84_bounds(args.aoi)
        LOGGER.info("Authenticating with Earthdata using .netrc")
        earthaccess.login(strategy="netrc")

        missing_ranges = contiguous_ranges(missing_raw_dates)
        LOGGER.info(
            "Searching %d missing raw date(s) in %d contiguous range(s)",
            len(missing_raw_dates),
            len(missing_ranges),
        )

        for range_start, range_end in missing_ranges:
            LOGGER.info("Searching %s through %s", range_start, range_end)
            results = earthaccess.search_data(
                count=args.max_results,
                short_name=PRODUCT,
                version=VERSION,
                bounding_box=bounds,
                temporal=(range_start.isoformat(), range_end.isoformat()),
            )
            LOGGER.info("Found %d matching granule(s)", len(results))
            if results:
                earthaccess.download(results, args.download_dir)

    granules = [
        path
        for path in discover_granules(args.download_dir)
        if granule_date(path) in dates_to_process
    ]

    if not granules:
        LOGGER.warning("No raw granules are available for the missing output dates")
        return 0

    available_dates = {granule_date(path) for path in granules}
    unavailable_dates = dates_to_process - available_dates
    if unavailable_dates:
        LOGGER.warning(
            "%d date(s) remain unavailable and will be retried next run: %s",
            len(unavailable_dates),
            describe_dates(unavailable_dates),
        )

    LOGGER.info("Processing %d local granule(s)", len(granules))
    scf_processor(
        [str(path) for path in granules],
        str(args.output_dir),
        res=args.resolution,
        extent_target=extent,
        epsg_target=str(args.target_epsg).removeprefix("EPSG:"),
        ow=args.overwrite,
        water_mask=args.water_mask,
    )

    return 0


def main(
    argv: Sequence[str] | None = None,
    scf_processor: Callable[..., object] | None = None,
) -> int:
    if scf_processor is None:
        if __package__:
            from .scf import get_scf_viirs
        else:
            from scf import get_scf_viirs

        scf_processor = get_scf_viirs

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        return run(args, scf_processor)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        LOGGER.error("%s", error)
        return 1


if __name__ == "__main__":
    sys.exit(main())
