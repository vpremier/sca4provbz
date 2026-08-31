#!/usr/bin/env python3
"""Convert VNP10A1F granules into consistently gridded SCF GeoTIFFs."""

from __future__ import annotations

import logging
import os
import re
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
import rasterio
import rioxarray  # noqa: F401 - registers the xarray ``rio`` accessor
import xarray as xr
from affine import Affine
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject


LOGGER = logging.getLogger(__name__)
PRODUCT = "VNP10A1F"
UNAVAILABLE_VALUE = 205.0
DEFAULT_WATER_MASK = (
    Path(__file__).resolve().parent.parent / "aux" / "Water_Mask_aligned.tif"
)
GRANULE_DATE_PATTERN = re.compile(r"^VNP10A1F\.A(?P<date>\d{7})(?:\.|$)")
SINUSOIDAL_RADIUS_PATTERN = re.compile(
    r"ProjParams=\(\s*(?P<radius>[0-9.]+)", re.IGNORECASE
)


def _decode_metadata(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _source_crs(hdf: h5py.File) -> CRS:
    """Build the sinusoidal CRS described by the HDF-EOS metadata."""
    try:
        metadata_value = hdf["HDFEOS INFORMATION"]["StructMetadata.0"][()]
    except KeyError as error:
        raise ValueError("HDF file has no HDF-EOS structural metadata") from error

    metadata = _decode_metadata(metadata_value)
    if "Projection=HE5_GCTP_SNSOID" not in metadata:
        raise ValueError("HDF grid does not declare the expected sinusoidal projection")

    radius_match = SINUSOIDAL_RADIUS_PATTERN.search(metadata)
    if radius_match is None:
        raise ValueError("HDF sinusoidal projection does not declare a sphere radius")
    radius = float(radius_match.group("radius"))
    return CRS.from_proj4(
        f"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R={radius} +units=m +no_defs"
    )


def _regular_spacing(values: np.ndarray, name: str) -> float:
    if values.ndim != 1 or values.size < 2:
        raise ValueError(f"{name} must be a one-dimensional coordinate array")
    differences = np.diff(values.astype(np.float64))
    spacing = float(differences[0])
    if spacing == 0 or not np.allclose(differences, spacing):
        raise ValueError(f"{name} coordinates are not regularly spaced")
    return spacing


def read_vnp10a1f(filename: str | Path) -> xr.Dataset:
    """Read the cloud-gap-filled NDSI layer and its georeferencing."""
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"Input granule does not exist: {path}")

    with h5py.File(path, "r") as hdf:
        try:
            grid = hdf["HDFEOS"]["GRIDS"]["VIIRS_Grid_IMG_2D"]
            field = grid["Data Fields"]["CGF_NDSI_Snow_Cover"]
            values = np.asarray(field)
            x_edges = np.asarray(grid["XDim"], dtype=np.float64)
            y_edges = np.asarray(grid["YDim"], dtype=np.float64)
            fill_value = int(np.asarray(field.attrs.get("_FillValue", 255)).item())
            crs = _source_crs(hdf)
        except KeyError as error:
            raise ValueError(f"Missing required VNP10A1F dataset in {path}") from error

    if values.shape != (y_edges.size, x_edges.size):
        raise ValueError(
            f"Data shape {values.shape} does not match coordinate dimensions "
            f"({y_edges.size}, {x_edges.size}) in {path}"
        )

    x_resolution = _regular_spacing(x_edges, "XDim")
    y_resolution = _regular_spacing(y_edges, "YDim")
    if x_resolution <= 0 or y_resolution >= 0:
        raise ValueError("Expected eastward XDim and north-to-south YDim coordinates")

    transform = Affine(
        x_resolution,
        0,
        float(x_edges[0]),
        0,
        y_resolution,
        float(y_edges[0]),
    )
    data = xr.DataArray(
        values,
        name="CGF_NDSI_Snow_Cover",
        dims=("y", "x"),
        coords={
            "x": x_edges + x_resolution / 2,
            "y": y_edges + y_resolution / 2,
        },
        attrs={"long_name": "Cloud-gap-filled NDSI snow cover"},
    )
    data.rio.write_crs(crs, inplace=True)
    data.rio.write_transform(transform, inplace=True)
    data.rio.write_nodata(fill_value, inplace=True)
    return data.to_dataset()


def _target_crs(epsg_target: str | int | CRS) -> CRS:
    if isinstance(epsg_target, CRS):
        return epsg_target
    value = str(epsg_target).strip()
    if value.isdigit():
        return CRS.from_epsg(int(value))
    return CRS.from_user_input(value)


def _reference_grid(image_path: str | Path) -> tuple[tuple[float, ...], CRS]:
    path = Path(image_path)
    if not path.is_file():
        raise FileNotFoundError(f"Reference image does not exist: {path}")
    with rasterio.open(path) as source:
        if source.crs is None:
            raise ValueError(f"Reference image has no CRS: {path}")
        bounds = source.bounds
        extent = (bounds.left, bounds.bottom, bounds.right, bounds.top)
        return extent, source.crs


def _validate_extent(
    extent: Sequence[float], resolution: float
) -> tuple[tuple[float, float, float, float], int, int]:
    if resolution <= 0 or not np.isfinite(resolution):
        raise ValueError("Output resolution must be a finite number greater than zero")
    if len(extent) != 4:
        raise ValueError("Extent must contain xmin, ymin, xmax, and ymax")

    xmin, ymin, xmax, ymax = (float(value) for value in extent)
    if not np.isfinite((xmin, ymin, xmax, ymax)).all():
        raise ValueError("Extent coordinates must be finite")
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("Extent must satisfy xmin < xmax and ymin < ymax")

    columns_exact = (xmax - xmin) / resolution
    rows_exact = (ymax - ymin) / resolution
    columns = round(columns_exact)
    rows = round(rows_exact)
    if not np.isclose(columns_exact, columns, rtol=0, atol=1e-7) or not np.isclose(
        rows_exact, rows, rtol=0, atol=1e-7
    ):
        raise ValueError(
            "Extent width and height must be exact multiples of the output resolution"
        )
    return (xmin, ymin, xmax, ymax), columns, rows


def _target_template(
    extent: Sequence[float], resolution: float, crs: CRS
) -> xr.DataArray:
    (xmin, _ymin, _xmax, ymax), columns, rows = _validate_extent(
        extent, resolution
    )
    x = xmin + (np.arange(columns, dtype=np.float64) + 0.5) * resolution
    y = ymax - (np.arange(rows, dtype=np.float64) + 0.5) * resolution
    template = xr.DataArray(
        np.zeros((rows, columns), dtype=np.uint8),
        coords={"y": y, "x": x},
        dims=("y", "x"),
    )
    template.rio.write_crs(crs, inplace=True)
    template.rio.write_transform(
        Affine(resolution, 0, xmin, 0, -resolution, ymax), inplace=True
    )
    return template


def _granule_output_date(filename: str | Path) -> str:
    match = GRANULE_DATE_PATTERN.match(Path(filename).name)
    if match is None:
        raise ValueError(f"Unrecognized {PRODUCT} filename: {Path(filename).name}")
    try:
        return datetime.strptime(match.group("date"), "%Y%j").strftime("%Y%m%d")
    except ValueError as error:
        raise ValueError(f"Invalid acquisition date in {Path(filename).name}") from error


def _conservative_invalid_mask(
    source: xr.DataArray, template: xr.DataArray
) -> xr.DataArray:
    """Mark targets touched by any invalid bilinear source contribution."""
    # Both 0 and 1 are real data here: GDAL must not exclude invalid pixels
    # from interpolation as it normally does for a source nodata value.
    invalid_source = np.asarray(
        (source.values < 0) | (source.values > 100), dtype=np.float32
    )
    # Initialize to invalid so pixels outside the source footprint remain masked.
    invalid_target = np.ones(template.shape, dtype=np.float32)
    reproject(
        source=invalid_source,
        destination=invalid_target,
        src_transform=source.rio.transform(),
        src_crs=source.rio.crs,
        src_nodata=None,
        dst_transform=template.rio.transform(),
        dst_crs=template.rio.crs,
        dst_nodata=1.0,
        resampling=Resampling.bilinear,
        init_dest_nodata=True,
    )
    # Bilinear weights are non-negative, so any positive result means at least
    # one invalid source pixel contributed with nonzero weight.
    return xr.DataArray(
        invalid_target > 0,
        coords=template.coords,
        dims=template.dims,
        name="invalid_contribution",
    )


def _read_water_mask(
    mask_path: str | Path, template: xr.DataArray
) -> xr.DataArray:
    """Read a water mask that is already aligned to the target grid."""
    path = Path(mask_path)
    if not path.is_file():
        raise FileNotFoundError(f"Water mask does not exist: {path}")

    with rasterio.open(path) as source:
        if source.count < 1:
            raise ValueError(f"Water mask has no raster bands: {path}")
        if source.crs is None:
            raise ValueError(f"Water mask has no coordinate reference system: {path}")
        if source.shape != template.shape:
            raise ValueError(
                f"Water mask shape {source.shape} does not match target shape "
                f"{template.shape}: {path}"
            )
        if source.crs != template.rio.crs:
            raise ValueError(f"Water mask CRS does not match the target grid: {path}")
        if source.transform != template.rio.transform():
            raise ValueError(
                f"Water mask transform does not match the target grid: {path}"
            )
        water = source.read(1) == 1

    return xr.DataArray(
        water,
        coords=template.coords,
        dims=template.dims,
        name="water_mask",
    )


def _write_raster_atomically(data: xr.DataArray, destination: Path) -> None:
    """Write and validate a TIFF before moving it into its final location."""
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.stem}.",
            suffix=".tmp.tif",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)

        data.rio.to_raster(
            temporary_path,
            driver="GTiff",
            dtype="float32",
            compress="DEFLATE",
            predictor=3,
            tiled=True,
        )
        with rasterio.open(temporary_path) as output:
            if output.width != data.sizes["x"] or output.height != data.sizes["y"]:
                raise RuntimeError("Written raster has unexpected dimensions")
            if output.crs != data.rio.crs:
                raise RuntimeError("Written raster has an unexpected CRS")
            if output.nodata is not None:
                raise RuntimeError("Written raster unexpectedly declares a no-data value")

        os.replace(temporary_path, destination)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def get_scf_viirs(
    fileList: Sequence[str | Path],
    outdir: str | Path,
    res: float = 500,
    img4ext: str | Path | None = None,
    extent_target: Sequence[float] | None = None,
    epsg_target: str | int | CRS | None = None,
    ow: bool = False,
    water_mask: str | Path | None = DEFAULT_WATER_MASK,
) -> list[Path]:
    """Reproject VNP10A1F granules and save snow-cover-fraction TIFFs."""
    output_directory = Path(outdir)
    output_directory.mkdir(parents=True, exist_ok=True)

    if img4ext is not None:
        extent, target_crs = _reference_grid(img4ext)
    else:
        if extent_target is None or epsg_target is None:
            raise ValueError(
                "Provide img4ext, or provide both extent_target and epsg_target"
            )
        extent = tuple(extent_target)
        target_crs = _target_crs(epsg_target)

    template = _target_template(extent, res, target_crs)
    water_target = (
        _read_water_mask(water_mask, template)
        if water_mask is not None
        else xr.zeros_like(template, dtype=bool).rename("water_mask")
    )
    LOGGER.info(
        "Water mask marks %d target pixel(s) as unavailable",
        int(water_target.sum().item()),
    )
    outputs: list[Path] = []
    for filename in fileList:
        input_path = Path(filename)
        output_date = _granule_output_date(input_path)
        output_path = output_directory / (
            "EURAC_SNOW_MERGE.alps.south-tyrol."
            f"{output_date}T120000.vnp10a1f.tif"
        )
        if output_path.exists() and not ow:
            LOGGER.info("Skipping existing output %s", output_path)
            continue

        LOGGER.info("Processing %s", input_path)
        source = read_vnp10a1f(input_path)["CGF_NDSI_Snow_Cover"]

        invalid_target = _conservative_invalid_mask(source, template)
        valid_target = ~invalid_target & ~water_target
        has_valid_data = bool(valid_target.any().item())

        ndsi = source.where((source >= 0) & (source <= 100)).astype(np.float32) / 100
        ndsi.rio.write_nodata(np.nan, inplace=True)
        ndsi_target = ndsi.rio.reproject_match(
            template, resampling=Resampling.bilinear, nodata=np.nan
        )

        scf = ((-0.01 + 1.45 * ndsi_target) * 100).clip(min=0, max=100)
        scf = scf.where(valid_target, other=UNAVAILABLE_VALUE).fillna(
            UNAVAILABLE_VALUE
        )
        scf = scf.astype(np.float32).rename("snow_cover_fraction")
        scf.attrs.update(
            {
                "long_name": "Snow cover fraction",
                "units": "percent",
                "unavailable_value": UNAVAILABLE_VALUE,
                "water_mask": str(water_mask) if water_mask is not None else "none",
            }
        )
        scf.rio.write_crs(target_crs, inplace=True)
        scf.rio.write_nodata(None, encoded=True, inplace=True)

        if not has_valid_data:
            LOGGER.warning(
                "%s has no valid pixels in the target area; writing an all-205 raster",
                input_path.name,
            )

        _write_raster_atomically(scf, output_path)
        outputs.append(output_path)

    return outputs
