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
from datetime import datetime as dt
import matplotlib.pyplot as plt
import progressbar



def dateFromFileName(string, 
                     date_formats=('%Y%m%dT%H%M%S', '%Y%m%d'), 
                     regex_patterns=(r'\d{8}T\d{6}', r'\d{8}')):
    """
    Extracts a date from a filename string based on specific date-time formats.

    Parameters:
    -----------
    string : str
        The filename (or string) from which to extract the date. It should contain 
        a date-time substring in one of the expected formats.
    date_formats : tuple, default is ('%Y%m%dT%H%M%S', '%Y%m%d')
        The possible date formats.
    regex_patterns : tuple, default is (r'\d{8}T\d{6}', r'\d{8}')
        Patterns matching the date formats.

    Returns:
    --------
    date : datetime.date
        A `datetime.date` object representing the extracted date (YYYY-MM-DD).

    Raises:
    -------
    ValueError : If no valid date is found in the string.
    """

    for regex_pattern, date_format in zip(regex_patterns, date_formats):
        match = re.search(regex_pattern, string)
        if match:
            date_string = match.group()
            return datetime.datetime.strptime(date_string, date_format).date()

    raise ValueError(f"No valid date found in the string: {string}")


    
def get_sca_metrics(array, snow, nosnow, cloud, nodata):
    """
    Computes Snow Cover Area (SCA) metrics from a classified raster array.

    Parameters:
    -----------
    array : np.ndarray
        A NumPy array representing the classified raster dataset, where each pixel 
        has a classification value (e.g., snow, no snow, cloud, or no data).
    snow : int or float or list
        The pixel value representing snow in the classification. If a list, it 
        is interpreted as a fractional snow cover map.
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
    cloudPixels = np.sum(array == cloud)
    nodataPixels = np.sum(array == nodata)

    if type(snow) is list:
        validPixels = np.sum(array <= snow[1]) 
        snowPixels = np.mean(array[array <= snow[1]]) * validPixels/100
        nosnowPixels = (100 - np.mean(array[array <= snow[1]])) * validPixels/100

    else:
        snowPixels = np.sum(array == snow)
        nosnowPixels = np.sum(array == nosnow)

    snowPixels_max = snowPixels + cloudPixels
    snowPixels_min = snowPixels


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
            if os.path.basename(snowMap).startswith('VNP10A1F'):
                # 'snow': 0-100, nosnow': 0, 'cloud': 205, 'nodata': 255
                d = get_sca_metrics(array, [0,100], 0, 205, 255)
            elif os.path.basename(snowMap).startswith('EURAC_SNOW'):
                # 'snow': 1, nosnow': 2, 'cloud': 3, 'nodata': 0
                d = get_sca_metrics(array, 1, 2, 3, 0)
            else: 
                raise ValueError(f"Insert classification values for: {snowMap}")


                
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



def load_stack_xarray(snowMap_fileNameList, shp_fileName=None):
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
    mask_values = [0, 3, 5, 255]
    
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



def get_scd(snowMap_fileNameList, date_start, date_end, shp_fileName=None, window=2):
    """
    Computes the Snow Cover Duration (SCD) from a list of snow cover maps.
    
    Parameters:
    -----------
    snowMap_fileNameList : list of str
        List of file paths for snow cover maps.
    date_start : str
        Start date in the format 'DD-MM-YYYY'.
    date_end : str
        End date in the format 'DD-MM-YYYY'.
    shp_fileName : str, optional
        Path to a shapefile for masking the data (default is None).
    window : int, optional
        Time window (in days) for extending the start and end dates (default is 2).
    
    Returns:
    --------
    xarray.DataArray
        A dataset representing the computed Snow Cover Duration (SCD) [0-1].
        Missing values are replaced with -999.
    """
    
    # Convert start and end dates to datetime.date objects
    date_start = dt.strptime(date_start, '%d-%m-%Y').date() - datetime.timedelta(days=window)
    date_end = dt.strptime(date_end, '%d-%m-%Y').date() + datetime.timedelta(days=window)
    
    # Extract dates and filter files
    filtered_files = [
        f for f in snowMap_fileNameList if date_start <= dateFromFileName(f) <= date_end
    ]
    
    # Sort files by extracted date
    sorted_files = sorted(filtered_files, key=dateFromFileName)
    
    snowMap_stack = load_stack_xarray(sorted_files, shp_fileName)
    
    # Create a complete time range since some dates might be missing
    full_time = pd.date_range(start=date_start, end=date_end, freq="D")

    # Reindex with full time and fill missing values with 255
    snowMap_stack_filled = snowMap_stack.reindex(time=full_time, 
                                                 fill_value=255)
    
    # MODIS
    if os.path.basename(snowMap_fileNameList[0]).startswith('EURAC_SNOW'):   
        # apply multi-temporal filter
        snowMap_stack_fltd = multitemporal_filter(snowMap_stack_filled, window=window)
    
        # keep only snow and no snow values
        # snow:1, nosnow:2
        snowMap_stack_fltd = xr.where(snowMap_stack_fltd == 1, 1, 
                                      xr.where(snowMap_stack_fltd == 2, 0, np.nan))
        
        # interpolate over time
        snowMap_stack_interp = snowMap_stack_fltd.interpolate_na(dim="time", method="linear")
        
        # compute snow cover duration (SCD)
        scd = snowMap_stack_interp.mean(dim='time')

    
    # VIIRS
    elif os.path.basename(snowMap_fileNameList[0]).startswith('VNP10A1F'):   
        # keep only snow and no snow values (0-100)
        # snowT = 30
        # snowMap_stack_fltd = xr.where((snowMap_stack_filled >= snowT) & (snowMap_stack_filled <= 100), 1, 
        #                       xr.where((snowMap_stack_filled >= 0) & (snowMap_stack_filled < snowT), 0, np.nan))
        
        snowMap_stack_fltd = xr.where((snowMap_stack_filled >= 0) & (snowMap_stack_filled <= 100), 
                                      snowMap_stack_filled, np.nan)
        
        # interpolate over time
        snowMap_stack_interp = snowMap_stack_fltd.interpolate_na(dim="time", method="linear")
        
        # compute snow cover duration (SCD)
        scd = snowMap_stack_interp.mean(dim='time')/100


    else: 
        raise ValueError("Product not recognized") 
               
    return scd 



def get_scd_statistics(snowMap_fileNameList, outdir, max_missing_days=30, 
                       shp_fileName=None, window=2, mode="yearly"):
    """
    Computes Snow Cover Duration (SCD) statistics for either full seasons (yearly) or individual months.

    Parameters:
    -----------
    snowMap_fileNameList : list of str
        List of file paths for snow cover maps.
    max_missing_days : int, optional
        Maximum allowed missing days in a time period (default is 30).
    shp_fileName : str, optional
        Path to a shapefile for masking the data (default is None).
    window : int, optional
        Time window (in days) for extending the start and end dates (default is 2).
    mode : str, optional
        - `"yearly"` (default) → Compute SCD for full snow seasons (1st Oct – 30th Sept).
        - `"monthly"` → Compute SCD for each individual month.

    Returns:
    --------
    dict
        Dictionary where keys are either **seasons (e.g., "2020-2021")** or **months ("YYYY-MM")** 
        and values are SCD datasets.
    """
    
    # Extract dates from filenames
    dates = sorted([dateFromFileName(f) for f in snowMap_fileNameList])
    
    results = {}

    # --- Yearly (Seasonal) Processing ---
    if mode == "yearly":
        for year in range(dates[0].year, dates[-1].year):  
            date_start = pd.Timestamp(year=year, month=10, day=1).date()  # 1st October
            date_end = pd.Timestamp(year=year + 1, month=9, day=30).date()  # 30th September
            
            # Filter files for this season
            season_files = [f for f in snowMap_fileNameList if date_start <= dateFromFileName(f) <= date_end]
            
            # Check missing days
            full_time_range = pd.date_range(start=date_start, end=date_end, freq="D")
            available_dates = set(dateFromFileName(f) for f in season_files)
            missing_days = len(full_time_range) - len(available_dates)
            
            if missing_days <= max_missing_days:  
                print(f"Processing season {year}-{year+1} (missing days: {missing_days})")
                results[f"{year}-{year+1}"] = get_scd(season_files, date_start.strftime('%d-%m-%Y'), 
                                                      date_end.strftime('%d-%m-%Y'), shp_fileName, window)
            else:
                print(f"Skipping season {year}-{year+1} (too many missing days: {missing_days})")
    
    # --- Monthly Processing ---
    elif mode == "monthly":
        for date in pd.date_range(start=dates[0], end=dates[-1], freq="MS"):  # MS = Month Start
            date_start = pd.Timestamp(year=date.year, month=date.month, day=1).date()
            date_end = pd.Timestamp(year=date.year, month=date.month, 
                                    day=pd.Period(date, freq='D').days_in_month).date()
            
            # Filter files for this month
            month_files = [f for f in snowMap_fileNameList if date_start <= dateFromFileName(f) <= date_end]
            
            # Check missing days
            full_time_range = pd.date_range(start=date_start, end=date_end, freq="D")
            available_dates = set(dateFromFileName(f) for f in month_files)
            missing_days = len(full_time_range) - len(available_dates)
            
            if missing_days <= max_missing_days:  
                print(f"Processing month {date.strftime('%Y-%m')} (missing days: {missing_days})")
                results[f"{date.strftime('%Y-%m')}"] = get_scd(month_files, date_start.strftime('%d-%m-%Y'), 
                                                               date_end.strftime('%d-%m-%Y'), shp_fileName, window)
            else:
                print(f"Skipping month {date.strftime('%Y-%m')} (too many missing days: {missing_days})")
    
    else:
        raise ValueError("Invalid mode! Use 'yearly' or 'monthly'.")

    return results


        