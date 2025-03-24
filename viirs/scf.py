#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:13:00 2023

@author: vpremier
"""
import os
import sys
import numpy as np
from datetime import datetime as dt
from osgeo import osr, gdal
import xarray as xr
import rioxarray
from affine import Affine
from rasterio.enums import Resampling
import h5py



def open_image(image_path):
    """
    Opens a raster image using GDAL and extracts geospatial metadata.

    Parameters:
    -----------
    image_path : str
        The file path to the raster image.

    Returns:
    --------
    tuple:
        - image : gdal.Dataset
            The opened raster image as a GDAL dataset.
        - information : dict
            A dictionary containing geospatial metadata:
            - 'geotransform' : tuple
                A six-element tuple describing the affine transformation:
                (top-left x, pixel width, rotation, top-left y, rotation, pixel height)
            - 'extent' : list
                A list defining the spatial extent of the raster: [minx, miny, maxx, maxy]
            - 'X_Y_raster_size' : list
                A list with the number of columns (width) and rows (height) of the raster: [cols, rows]
            - 'projection' : str
                The projection information in WKT (Well-Known Text) format.
    """
    image = gdal.Open(image_path)

    if image is None:
        print('Could not open ' + image_path)
        sys.exit(1)

    cols = image.RasterXSize
    rows = image.RasterYSize
    geotransform = image.GetGeoTransform()
    proj = image.GetProjection()
    
    minx = geotransform[0]
    maxy = geotransform[3]
    maxx = minx + geotransform[1] * cols
    miny = maxy + geotransform[5] * rows

    X_Y_raster_size = [cols, rows]
    extent = [minx, miny, maxx, maxy]

    information = {
        'geotransform': geotransform,
        'extent': extent,
        'X_Y_raster_size': X_Y_raster_size,
        'projection': proj
    }

    return image, information



def read_vnp10a1f(filename):
    """
    Reads a VNP10A1F HDF5 file and converts it into an xarray.Dataset
    containing NDSI-related variables.

    Parameters:
    -----------
    filename : str
        Path to the HDF5 file.

    Returns:
    --------
    ds : xarray.Dataset
        An xarray dataset containing:
        - 'CGF_NDSI_Snow_Cover': DataArray representing the CGF_NDSI_Snow_Cover field.
        - 'Daily_NDSI_Snow_Cover': DataArray representing the Daily_NDSI_Snow_Cover field.
    
    Dataset Attributes:
    -------------------
    - The dataset includes spatial coordinates (x, y).
    - A sinusoidal projection is defined for geospatial reference.
    - The transform and CRS (Coordinate Reference System) are embedded in the dataset using rasterio.

    Notes:
    ------
    - The function reads two variables: 'CGF_NDSI_Snow_Cover' and 'Daily_NDSI_Snow_Cover'.
    - The projection information is based on a sinusoidal projection.
    - Uses the rasterio and affine transformations for proper geospatial referencing.

    For details about the dataset, check: https://nsidc.org/data/vnp10a1f/versions/2
    """
    # Open the HDF5 file
    f = h5py.File(filename, 'r')

    # Read the variables of interest
    CGF_NDSI_Snow_Cover = np.array(f['HDFEOS']['GRIDS']['VIIRS_Grid_IMG_2D']['Data Fields']['CGF_NDSI_Snow_Cover'])
    Daily_NDSI_Snow_Cover = np.array(f['HDFEOS']['GRIDS']['VIIRS_Grid_IMG_2D']['Data Fields']['Daily_NDSI_Snow_Cover'])

    # Define the projection (sinusoidal)
    projInfo = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid", DATUM["Not specified (based on custom spheroid)", SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]',\
               'GEOGCS["Unknown datum based upon the Clarke 1866 ellipsoid", DATUM["Not specified (based on Clarke 1866 spheroid)", SPHEROID["Clarke 1866",6378206.4,294.9786982139006]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]]'

    # Read x and y coordinates
    XDim = np.array(f['HDFEOS']['GRIDS']['VIIRS_Grid_IMG_2D']['XDim'])
    YDim = np.array(f['HDFEOS']['GRIDS']['VIIRS_Grid_IMG_2D']['YDim'])

    # Compute geotransform
    geotransform = (XDim[0], XDim[1] - XDim[0], 0, YDim[0], 0, YDim[1] - YDim[0])

    # Function to prepare a DataArray
    def prepare_dataarray(array, varname):
        da = xr.DataArray(
            name=varname,
            data=array,
            dims=["y", "x"],
            coords=dict(
                y=(["y"], YDim + (YDim[1] - YDim[0]) / 2),
                x=(["x"], XDim + (XDim[1] - XDim[0]) / 2)
            ),
            attrs=dict(
                transform=Affine.from_gdal(*geotransform),
                crs=projInfo
            ),
        )
        return da

    # Create xarray Dataset
    ds = xr.Dataset({
        'CGF_NDSI_Snow_Cover': prepare_dataarray(CGF_NDSI_Snow_Cover, 'CGF_NDSI_Snow_Cover'),
        'Daily_NDSI_Snow_Cover': prepare_dataarray(Daily_NDSI_Snow_Cover, 'Daily_NDSI_Snow_Cover')
    })

    # Set CRS and transform using rasterio
    ds.rio.write_crs(projInfo, inplace=True) \
        .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True) \
        .rio.write_coordinate_system(inplace=True) \
        .rio.write_transform(Affine.from_gdal(*geotransform))

    return ds



def get_scf_viirs(fileList, outdir, res=500, img4ext=None, extent_target=None, 
                  epsg_target=None, ow=False):
    """
    Reads VIIRS snow cover fraction (SCF) data from HDF files, performs 
    resampling and reprojection, and saves the output as GeoTIFF.

    Parameters:
    -----------
    fileList : list
        List of paths to the input HDF files.
    outdir : str
        Directory where the output GeoTIFF files will be saved.
    res : int, optional (default=500)
        Spatial resolution (in meters) of the output image.
    img4ext : str, optional
        Path to an image file for reading the target extent and CRS.
    extent_target : tuple, optional
        Manually specified extent in the format (xMin, yMin, xMax, yMax).
    epsg_target : str, optional
        EPSG code of the target coordinate system.
    ow : bool, optional (default=False)
        Whether to overwrite existing output files.

    Returns:
    --------
    None : NoneType
        Outputs are saved as GeoTIFF files in the specified directory.

    Notes:
    ------
    - If `img4ext` is provided, it determines the extent and EPSG from the given image.
    - If `img4ext` is not provided, `extent_target` and `epsg_target` must be specified.
    - The function reprojects and resamples the NDSI snow cover data.
    - Cloud mask values are retained as 205 in the output.
    - The snow cover fraction is computed using: `SCF = (-0.01 + 1.45 * NDSI) * 100`
      with values capped between 0 and 100.
    """
    
    if img4ext:
        print('Reading extent and EPSG from reference image...')
        ds, info = open_image(img4ext)
        
        # Extract target EPSG
        srOut = osr.SpatialReference(str(info['projection']))
        epsg = 'EPSG:' + srOut.GetAttrValue("AUTHORITY", 1)   
        
    else:
        if epsg_target is None:
            print('Please specify the target EPSG or provide a reference image.')
            return
        if extent_target is None:
            print('Please specify the target extent or provide a reference image.')
            return
        epsg = 'EPSG:' + epsg_target
        info = {'extent': list(extent_target)}

    # Compute grid dimensions
    nx = int((info['extent'][2] - info['extent'][0]) / res)
    ny = int((info['extent'][3] - info['extent'][1]) / res)

    # Create target grid
    x = np.linspace(info['extent'][0] + res / 2, info['extent'][2] - res / 2, nx)
    y = np.flip(np.linspace(info['extent'][1] + res / 2, info['extent'][3] - res / 2, ny))
    
    # Initialize an empty dataset
    ones = np.ones((ny, nx))
    output_da = xr.DataArray(ones, coords=[y, x], dims=["y", "x"])
    output_da.rio.write_crs(epsg, inplace=True)

    for filename in fileList:
        print(f'Processing {filename}')
        
        # Generate output filename
        date = os.path.basename(filename).split('.')[1][1:] 
        new_date = dt.strptime(date, '%Y%j').strftime('%Y%m%d')
        fileName_output = os.path.join(outdir, f'VNP10A1F_{new_date}.tif')

        # Skip existing files if overwrite is False
        if os.path.exists(fileName_output) and not ow:
            print(f'File {fileName_output} already exists. Set `ow=True` to overwrite.')
            continue  
        
        # Read VIIRS dataset
        viirs_ds = read_vnp10a1f(filename)
        
        # Reproject NDSI Snow Cover
        ndsi = viirs_ds['CGF_NDSI_Snow_Cover'].rio.reproject(epsg, resampling=Resampling.bilinear) / 100
        
        # Reproject cloud mask
        cloud = viirs_ds['CGF_NDSI_Snow_Cover'].rio.reproject(epsg, resampling=Resampling.nearest)
        cloud_rsmp = cloud.rio.reproject_match(output_da, resampling=Resampling.nearest)
        
        # Skip files with only invalid values
        if (cloud_rsmp.values[cloud_rsmp.values > 100] == 0).all():
            continue
        
        # Resample NDSI
        ndsi_rsmp = ndsi.rio.reproject_match(output_da, resampling=Resampling.cubic)

        # Compute Snow Cover Fraction (SCF)
        scf = (-0.01 + 1.45 * ndsi_rsmp) * 100
        scf = scf.where(scf < 100, other=100).where(scf > 0, other=0)
        scf = scf.where(cloud_rsmp <= 100, other=205)

        # Save output as GeoTIFF
        scf.rio.to_raster(fileName_output)


