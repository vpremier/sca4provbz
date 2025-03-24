#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 13:16:44 2025

@author: vpremier
"""

import os
import numpy as np
import xarray as xr
import rioxarray  # Ensure rioxarray is imported to enable raster operations
import geopandas as gpd
import pandas as pd
from osgeo import gdal, osr
from shapely.geometry import mapping
from shapely.ops import unary_union  
import re
import datetime
import matplotlib.pyplot as plt
import progressbar



def dateFromFileName(string, date_format = '%Y%m%dT%H%M%S', 
                     regex_pattern = r'\d{8}T\d{6}' ):
    """
    Extracts a date from a filename string based on a specific date-time format.
    
    Parameters:
    -----------
    string : str
        The filename (or string) from which to extract the date. It should contain a date-time substring
        in the format 'YYYYMMDDTHHMMSS'.
    date_format : str, default is %Y%m%dT%H%M%S
        The date format 
    regex_pattern : str, default is r'\d{8}T\d{6}'
        The pattern that will match the general format of 'YYYYMMDDTHHMMSS or 
        another pattern defined by the user.
        
    Returns:
    --------
    date : datetime.date
        A `datetime.date` object representing the extracted date (YYYY-MM-DD).
    
    Raises:
    -------
    ValueError : If the date format is not found in the string or is malformed.
    """
    # Search for the date-time string
    match = re.search(regex_pattern, string)
    
    if not match:
        raise ValueError(f"Date-time not found in the string: {string}")
    
    # Extract the date-time string
    date_string = match.group()
    
    # Convert the extracted string to a datetime object and return the date part
    date = datetime.datetime.strptime(date_string, date_format).date()
    return date



def get_sca_metrics(array, snow, nosnow, cloud, nodata):
    """
    Computes Snow Cover Area (SCA) metrics from a classified raster array.

    Parameters:
    -----------
    array : np.ndarray
        A NumPy array representing the classified raster dataset, where each pixel 
        has a classification value (e.g., snow, no snow, cloud, or no data).
    snow : int or float
        The pixel value representing snow in the classification.
    nosnow : int or float
        The pixel value representing no snow (bare ground or other land cover).
    cloud : int or float
        The pixel value representing clouds in the dataset.
    nodata : int or float
        The pixel value representing missing or invalid data.

    Returns:
    --------
    dict
        A dictionary containing the following computed metrics:
        - `cloudPercent`: Percentage of the image covered by clouds.
        - `SCA_cloudFree`: Snow Cover Area (SCA) percentage considering only cloud-free pixels.
        - `SCA`: Snow Cover Area percentage including clouds.
        - `SCA_error`: Estimated error in the snow percentage due to cloud coverage.
    """
    # Compute pixel counts
    cloudPixels = np.sum(np.logical_or(array == cloud, array == nodata))
    snowPixels = np.sum(array == snow)
    snowPixels_max = snowPixels + cloudPixels
    snowPixels_min = snowPixels
    nosnowPixels = np.sum(array == nosnow)
    nodataPixels = np.sum(array == nodata)

    # Compute total number of pixels
    N = float(snowPixels + cloudPixels + nosnowPixels + nodataPixels)
    if N == 0:
        raise ValueError("The input array contains no valid pixels for computation.")

    # Compute metrics
    cloudPercent = (cloudPixels / N) * 100
    snowPixels_error = float(snowPixels_max - snowPixels_min) / 2.0
    snowPercent_error = (snowPixels_error / N) * 100
    snowPercent = ((snowPixels_min + snowPixels_error) / N) * 100

    # Compute cloud-free snow percentage
    noCloudPixels = snowPixels + nosnowPixels
    snowPercent_cloudFree = (snowPixels / float(noCloudPixels) * 100) if noCloudPixels > 0 else -1

    # Return results as a dictionary
    return {
        'cloudPercent': cloudPercent,
        'SCA_cloudFree': snowPercent_cloudFree,
        'SCA': snowPercent,
        'SCA_error': snowPercent_error
    }



def apply_mask(snowMap, shp_fileName=None):
    """
    Applies a mask to a raster dataset using a shapefile. 
    Masked values are set to 255.

    Parameters:
    -----------
    snowMap : str
        Path to the raster file (GeoTIFF).
    shp_fileName : str, optional
        Path to the shapefile used for masking (default is None).

    Returns:
    --------
    xarray.DataArray
        Masked raster dataset clipped to the shapefile extent.

    """
    # Open the raster to get CRS
    img = gdal.Open(snowMap)
    if img is None:
        raise ValueError(f"Could not open raster file: {snowMap}")

    prj_img = img.GetProjection()
    crs_img = osr.SpatialReference(wkt=prj_img)
    img = None  # Close file to free memory

    # Open the raster using rioxarray
    ds = rioxarray.open_rasterio(snowMap)

    # Open and process the shapefile if provided
    if shp_fileName is not None:
        shp = gpd.read_file(shp_fileName)

        if shp.empty:
            raise ValueError(f"Shapefile {shp_fileName} is empty.")

        # Convert to the same CRS as the raster
        shp_rpj = shp.to_crs(crs=crs_img.ExportToProj4())

        # Merge multiple polygons into one if there are multiple layers or features
        merged_geometry = unary_union(shp_rpj["geometry"])

        # Clip the raster using the merged polygon
        ds = ds.rio.write_nodata(255)
        ds = ds.rio.clip([mapping(merged_geometry)], shp_rpj.crs, drop=True)

        
    return ds



def updateCSV(csv_path, snowMap_fileNameList, shp_fileName=None):
    """
    This function updates or creates a CSV with new data based on files in a given folder.

    Parameters:
    -----------
    csv_path : str
        Path to the existing CSV file to be updated or created.
    snowMap_fileNameList : list of str
        List of snow map file names to process.
    work_folder : str
        Path to the working folder for intermediate results.
    shp_fileName : str, optional
        Path to the shapefile, if needed for spatial processing.
    """
    
    # Initialize the DataFrame and date range
    df = None
    new_data = []

    # Check if the CSV file exists and load it
    if os.path.exists(csv_path):
        print('The CSV file has already been created and will be updated.')
        
        try:
            df = pd.read_csv(csv_path)
            df.index = pd.to_datetime(df['Unnamed: 0'], format='%Y-%m-%d') if 'Unnamed: 0' in df.columns else None
            df.drop(columns='Unnamed: 0', errors='ignore', inplace=True)  # Remove the 'Unnamed: 0' column if it exists
            dateStart = df.index[-1].date() if df is not None else datetime.datetime.min.date()
            print(f'The old CSV contains data till {dateStart}')
            
            # Filter files based on the date in the filename
            snowMap_fileNameList = [f for f in snowMap_fileNameList if dateFromFileName(f) > dateStart]
        
        except Exception as e:
            print(f"Error reading or processing CSV: {e}")
            return
    
    else:
        print('The CSV file does not exist. A new one will be created.')
    
    # Process the new files
    for snowMap in snowMap_fileNameList:
        date = dateFromFileName(snowMap)
        print(f'Processing date: {date}')
        
        try:
            # Apply mask to the snow map using the shapefile
            ds = apply_mask(snowMap, shp_fileName=shp_fileName)
            array = ds.sel(band=1).values  # Assuming the first band is of interest
            
            # Get the metrics from the array
            d = get_sca_metrics(array, 1, 2, 3, 0)
            new_data.append((date, d))
        
        except Exception as e:
            print(f"Error processing {snowMap}: {e}")
            continue  # Skip the current file and continue with the next
    
    if not new_data:
        print("No new data to add.")
        return df

    # Create a DataFrame from the collected new data
    newdf = pd.DataFrame(data=[item[1] for item in new_data], index=[item[0] for item in new_data])
    newdf.index = pd.to_datetime(newdf.index, format='%Y-%m-%d')
    newdf.sort_index(inplace=True)
    
    # Combine the old and new DataFrames
    if df is not None:
        updated_df = pd.concat([df, newdf], axis=0)
    else:
        updated_df = newdf

    # Write the updated data to the CSV file
    updated_df.to_csv(csv_path)
    print(f"CSV file {csv_path} has been updated")
    
    return updated_df
    
    

def snow_bullettin(pathToCSV, date_start, date_end, work_folder, CCA_threshold=30):
    """
    Plots the daily Snow Cover Area (SCA) based on a given dataset, filtering out 
    days with excessive cloud coverage. The function also computes historical mean, 
    minimum, and maximum SCA values for comparison.

    Parameters:
    -----------
    pathToCSV : str
        Path to the CSV file containing the SCA data.
    date_start : str
        Start date for analysis in 'YYYY-MM-DD' format.
    date_end : str
        End date for analysis in 'YYYY-MM-DD' format.
    work_folder : str
        Directory where the plot will be saved.
    CCA_threshold : int, optional
        Maximum allowed cloud percentage for a day to be considered valid (default is 30%).

    Returns:
    --------
    newdf : pandas.DataFrame
        A DataFrame containing the SCA values for the observation period, 
        along with historical mean, min, and max values.
        
    Output:
    -------
    - Saves a plot of the daily SCA trend compared to historical statistics.
    """

    # Read CSV
    df = pd.read_csv(pathToCSV)
    df.index = pd.to_datetime(df['Unnamed: 0'], format='%Y-%m-%d')
    
    # filter up to the end date
    df = df[:date_end]
    
    # Ensure all dates are present
    df = df.reindex(pd.date_range(df.index[0], df.index[-1]))
    
    # Calculate number of days in observation period
    date_start_dt = pd.to_datetime(date_start)
    date_end_dt = pd.to_datetime(date_end)
    nr_days = (date_end_dt - date_start_dt).days
    
    # Extract day, month, year
    df['day'] = df.index.day
    df['month'] = df.index.month
    df['year'] = df.index.year

    # Apply cloud threshold filter
    df['correctedSCA'] = np.nan
    df.loc[df['cloudPercent'] <= CCA_threshold, 'correctedSCA'] = df['SCA_cloudFree']
    df['correctedSCA'] = df['correctedSCA'].interpolate()
    
    # Select data for given period
    sca = df.loc[date_start:date_end, 'correctedSCA']
    
    # Compute historical statistics (excluding the observation period)
    historical_data = df[df.index < date_start]  # Exclude current period
    # daily_stats = historical_data.groupby(['month', 'day'])['correctedSCA'].agg(['mean', 'min', 'max'])
    daily_stats = historical_data.groupby(['month', 'day'])['correctedSCA'].agg(
    mean='mean',
    min=lambda x: np.percentile(x.dropna(), 10),
    max=lambda x: np.percentile(x.dropna(), 90)
    )

    # Reorder index to start from October 1st
    year_start = datetime.date(int(date_start[:4]), 1, 1)
    index_shift = (year_start - date_start_dt.date()).days + 365
    daily_stats.index = list(range(index_shift, 367)) + list(range(1, index_shift))
    
    # Remove Feb 29 if exists
    daily_stats = daily_stats[~((daily_stats.index == 60) & (~df.index.is_leap_year.any()))]
    daily_stats.sort_index(inplace=True)
    
    # Create DataFrame for plotting
    newdf = pd.DataFrame(index=sca.index)
    newdf['current SCA'] = sca.values
    newdf['min SCA'] = daily_stats['min'].iloc[:nr_days+1].values
    newdf['max SCA'] = daily_stats['max'].iloc[:nr_days+1].values
    newdf['mean SCA'] = daily_stats['mean'].iloc[:nr_days+1].values
    
    newdf.sort_index(inplace=True)


    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(sca.index, sca, label='Current SCA')
    plt.plot(sca.index, newdf['mean SCA'], label='Mean SCA', linestyle='dashed')
    plt.fill_between(sca.index, newdf['min SCA'], newdf['max SCA'], color='b', alpha=0.2)
    plt.ylim([0, 100])
    plt.xticks(rotation=30)
    plt.ylabel('% SCA')
    plt.legend()
    plt.grid()
    
    # Save plot
    plt.savefig(os.path.join(work_folder, f"snow_bulletin_{date_start}_{date_end}.png"), bbox_inches="tight")
    
    return newdf



def cloud_temporalFilter(fileList, dateList, cloud, nodata, output_root, water=None, waterMask_fileName=None, confidence=False,tmax=3):

    # This function generate a cloud filtered map considering a time window of +- 2 days. Only pixel having snow
    # (or no snow) in one image before and after inside this window will be cleared from cloud

    # INPUTS:
    # fileList: list of file names complete of path
    # dateList: list of dates corresponding to the fileList (same order) in datetime
    # cloud, nodata: cloud and nodata value
    # output_root: output root name complete of path (e.g. /.../snowMap_). Date and '_cloudTempFil.tif' will be added to
    #              the fileName of each map generated
    # confidence: True if the confidence layer is present in band 2, False if confidence layer is not present
    #
    # OPTIONAL INPUTS:
    # water: water value of the input and output
    # waterMask_fileName: water mask filename complete of path (boolean raster)


    # Read info from the first image
    img_ref = gdal.Open(fileList[0])
    Ncols = img_ref.RasterXSize
    Nrows = img_ref.RasterYSize
    geoTransform = img_ref.GetGeoTransform()
    projection = img_ref.GetProjection()
    colorTable = img_ref.GetRasterBand(1).GetColorTable()
    #img = None

    # Open the water mask
    if waterMask_fileName is not None:
        img = gdal.Open(waterMask_fileName)
        water_mask = img.GetRasterBand(1).ReadAsArray()
        img = None


    # Run the temporal filter:
    bar = progressbar.ProgressBar(len(dateList))
    bar.start()
    
    
    def fuse_adj_snowmap(file_list, dim, cloud, nodata):

        # File list must be ordered by date
        fused = np.zeros(dim,dtype=np.uint8) + cloud
        for n in file_list:
            img = gdal.Open(n)
            snowMap = img.GetRasterBand(1).ReadAsArray()
            cloudMask = np.logical_or(fused==cloud, fused==nodata)
            fused[cloudMask] = snowMap[cloudMask]

        return fused

    for i,d in enumerate(dateList):

        print( "Elaborating date: " + d.strftime('%Y%m%d%H%M%S'))

        # Fuse the images acquired max tmax days before the current date d

        # Find the list of pre-images to be fused
        indexes = [n for n,x in enumerate(dateList) if x<d and x>=d-datetime.timedelta(days=tmax)]
        if not indexes:
            print( "Empty list!")

        file_list = list(reversed(sorted([fileList[ix] for ix in (indexes)])))                
        fused_pre = fuse_adj_snowmap(file_list, (Nrows, Ncols), cloud, nodata)

        indexes = [n for n,x in enumerate(dateList) if x>d and x<=d+timedelta(days=tmax)]
        if not indexes:
            print ("Empty list!")

        file_list = (sorted([fileList[ix] for ix in indexes]))
        fused_post = fuse_adj_snowmap(file_list, (Nrows, Ncols), cloud, nodata)
        



        # Load the snow map of the time 0 and apply the temporal filter
        snowMap_curr, conf_curr = load_snowMap_curr(d, dateList, fileList, (Nrows,Ncols), cloud, confidence)
        if snowMap_curr is None:
            bar.update(i + 1)
            continue
        cloudMask = np.logical_or(snowMap_curr==cloud, snowMap_curr==nodata)
        if np.sum(cloudMask):
            snowMap_fill = np.logical_and(cloudMask, fused_pre==fused_post)
            snowMap_curr[snowMap_fill] = fused_pre[snowMap_fill]

        # Mask and write the filtered map
        if np.sum(snowMap_curr == cloud) < Nrows*Ncols:

            # Mask the snow map with the water bodies
            if water is not None:
                if waterMask_fileName is not None:
                    snowMap_curr[water_mask] = water
                else:
                    water_mask_curr = np.logical_or(snowMap_pre==water,snowMap_post==water)
                    snowMap_curr[water_mask_curr] = water

            # Write the filtered image to file
            output_fileName = os.path.join(output_root + d.strftime('%Y%m%d%H%M%S') + '_sca_fused.tif')
            img = gdal.GetDriverByName('GTiff').Create(output_fileName, Ncols, Nrows, 1, gdal.GDT_Byte)
            img.SetGeoTransform(geoTransform)
            img.SetProjection(projection)
            img.GetRasterBand(1).WriteArray(snowMap_curr)
            if colorTable is not None:
                img.GetRasterBand(1).SetColorTable(colorTable.Clone())
            img = None

            if confidence:
                conf_curr[snowMap_fill] = conf_pre[snowMap_fill]
                changeMask = np.logical_and(snowMap_fill, conf_post<conf_pre)
                conf_curr[changeMask] = conf_post[changeMask]

                # Write the filtered image to file
                outputConf_fileName = output_fileName[:-4] + '_conf.tif'
                img = gdal.GetDriverByName('GTiff').Create(outputConf_fileName, Ncols, Nrows, 1, gdal.GDT_Float32)
                img.SetGeoTransform(geoTransform)
                img.SetProjection(projection)
                img.GetRasterBand(1).WriteArray(conf_curr)
                img = None

                # Create the vrt
                cmd = "gdalbuildvrt -separate " + output_fileName[:-4] + '_stack.vrt ' + output_fileName + " " + outputConf_fileName
                os.system(cmd)

        bar.update(i + 1)
    bar.finish()
    img_ref = None
 
    

def load_stack_xarray(snowMap_fileNameList, shp_fileName):
    """
    Load and stack masked snow maps into an xarray DataArray with a time dimension.

    Parameters:
    - snowMap_fileNameList (list of str): List of file paths for snow maps.
    - shp_fileName (str): Path to shapefile used for masking.

    Returns:
    - xarray.DataArray: A 3D DataArray (time, y, x) containing the stacked masked maps.
    """
    
    # Ensure the list is sorted alphabetically
    snowMap_fileNameList = sorted(snowMap_fileNameList)

    data_list = []
    time_list = []

    for snowMap in snowMap_fileNameList:
        date = dateFromFileName(snowMap)  # Extract date from filename
        print(f"Elaborating date: {date.strftime('%Y-%m-%d')}")
        
        time_list.append(pd.to_datetime(date))  # Convert date to pandas datetime format

        # Apply mask to the snow map using the shapefile
        ds = apply_mask(snowMap, shp_fileName=shp_fileName)

        # Select the first band and convert to NumPy array
        try:
            array = ds.sel(band=1, drop=True).values  # Drop band dimension for cleaner output
        except KeyError:
            print(f"Warning: Band 1 not found in {snowMap}. Skipping...")
            continue

        data_list.append(array)

    if not data_list:
        return None  # Return None if no valid data

    # Stack arrays along a new time dimension
    data_array = xr.DataArray(
        np.stack(data_list),
        dims=["time", "y", "x"],
        coords={"time": time_list, "y": ds.y, "x": ds.x},
        name="snow_map"
    )

    return data_array








def multitemporal_filter(data_array, window=2):
    """
    Applies a multitemporal filter to an xarray DataArray, checking a window of ±2 days.
    
    Parameters:
    - data_array (xarray.DataArray): Input time-series data.
    - window (int): Number of days for forward and backward filling (default: 2).
    
    Returns:
    - xarray.DataArray: Filtered DataArray with NaNs where conditions are met.
    """
    # Define values to mask
    mask_values = [0, 3, 255]
    
    # Create mask
    mask = ~data_array.isin(mask_values)
    masked_data = data_array.where(mask, np.nan)
    
    # Forward and backward fill
    data_ffill = masked_data.ffill(dim='time', limit=window)
    data_bfill = masked_data.bfill(dim='time', limit=window)

    # Copy original data to avoid in-place modifications
    filtered_data = data_array.copy()

    for date in masked_data.time:
        # Condition 1: forward and backward filled values must be the same
        mask_mt = data_bfill.sel(time=date) == data_ffill.sel(time=date)
        
        # Condition 2: Original data was masked (was in mask_values)
        mask_nan = ~mask.sel(time=date)
        
        # Combine all conditions
        mask_new = np.logical_and(mask_mt, mask_nan)

        # Apply the filtered values
        filtered_data.loc[dict(time=date)] = xr.where(mask_new, data_bfill.sel(time=date), data_array.sel(time=date))

    return filtered_data



def multitemporal_filter(data_array, window=2):
    
    # Define values to mask
    mask_values = [0, 3, 255]
    
    # Apply mask
    mask = ~data_array.isin(mask_values)
    masked_data = data_array.where(mask, np.nan)
    
    data_ffill = masked_data.ffill(dim='time', limit = window)
    data_bfill = masked_data.bfill(dim='time', limit = window)



    for i,date in enumerate(masked_data.time):
        
        # condizione 1 il gap deve essere inferiore o uguale a 5 giorni
        delta = np.array(data_bfill.sel(time=date)) - np.array(data_ffill.sel(time=date)) 
        mask_delta = delta <= 30
        
        # il forward fill e backward fill devono essere uguali
        mask_mt = np.array(data_bfill.sel(time=date)) == np.array(data_ffill.sel(time=date))
        mask_nan = np.array(~mask.sel(time=date))
        
        # deve essere nan
        mask_new = np.logical_and.reduce((mask_mt,mask_nan))
 
        
        data_array.sel(time=date).values = np.array(data_bfill.sel(time=date))[mask_new] 
        

    


        