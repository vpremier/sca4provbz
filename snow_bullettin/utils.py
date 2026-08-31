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
from collections import defaultdict
import pymannkendall as mk



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
    cloudPixels = np.sum(array == cloud) + np.sum(array == nodata)

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
    N = float(snowPixels + cloudPixels + nosnowPixels)
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



def apply_mask(snowMap, shp_fileName=None, mask_raster_fileName=None):
    """
    Applies a mask to a raster dataset using a shapefile and/or another raster mask.
    Masked values are set to 255.

    Parameters:
    -----------
    snowMap : str
        Path to the raster file (GeoTIFF).
    shp_fileName : str, optional
        Path to the shapefile for masking.
    mask_raster_fileName : str, optional
        Path to a raster file used as an additional mask.

    Returns:
    --------
    xarray.DataArray
        Masked raster dataset.
    """
   

    # Open the raster to get CRS
    img = gdal.Open(snowMap)
    if img is None:
        raise ValueError(f"Could not open raster file: {snowMap}")

    prj_img = img.GetProjection()
    crs_img = osr.SpatialReference(wkt=prj_img)
    img = None  # Close GDAL dataset

    # Open the raster with rioxarray
    ds = rioxarray.open_rasterio(snowMap)

    # Ensure nodata is defined
    ds = ds.rio.write_nodata(255)

    # --- Apply vector mask ---
    if shp_fileName is not None:
        shp = gpd.read_file(shp_fileName)

        if shp.empty:
            raise ValueError(f"Shapefile {shp_fileName} is empty.")

        # Reproject shapefile to raster CRS
        shp_rpj = shp.to_crs(crs=crs_img.ExportToProj4())

        # Merge all geometries
        merged_geometry = unary_union(shp_rpj["geometry"])

        # Clip by shapefile geometry
        ds = ds.rio.clip([mapping(merged_geometry)], shp_rpj.crs, drop=True)
        ds = ds.rio.write_crs(prj_img, inplace=True)

    # --- Apply raster mask ---
    if mask_raster_fileName is not None:
        mask_ds = rioxarray.open_rasterio(mask_raster_fileName)

        # Resample mask raster to match the main raster shape/resolution if needed
        if ds.shape[1:] != mask_ds.shape[1:]:
            mask_ds = mask_ds.rio.reproject_match(ds)

        # Build mask: here we assume mask=0 means keep, mask>0 means mask out.
        # Adapt this logic to your mask raster values!
        mask_array = mask_ds.sel(band=1).values

        # Apply mask to main raster: where mask is 1 (or >0), set raster to nodata (255)
        ds_array = ds.sel(band=1).values
        ds_array[mask_array== 0] = 255

        # Replace the band values with the masked array
        ds[0] = ds_array

    return ds


                
def updateCSV(csv_path, snowMap_fileNameList, shp_fileName=None, mask_raster_fileName=None):
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
    mask_raster_fileName : str, optional
        Path to a raster file to use as an additional mask.
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
            ds = apply_mask(snowMap, shp_fileName=shp_fileName, mask_raster_fileName=mask_raster_fileName)
            array = ds.sel(band=1).values  # Assuming the first band is of interest
            
            # Get the metrics from the array
            if 'vnp10a1f' in os.path.basename(snowMap):
                # 'snow': 0-100, nosnow': 0, 'cloud': 205, 'nodata': 205
                d = get_sca_metrics(array, [0,100], 0, 205, 205)
            elif 'modis' in os.path.basename(snowMap):
                # 'snow': 1, nosnow': 2, 'cloud': 3, 'nodata': 0
                d = get_sca_metrics(array, 1, 2, 3, 0)
            elif 'complete' in os.path.basename(snowMap).lower():
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
    
    

def snow_bullettin(pathToCSV, 
                   date_start, 
                   date_end, 
                   work_folder, 
                   suffix, 
                   CCA_threshold=30):
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
    suffix : str
        Suffix for the image path.
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
    
    # Remove Feb 29 from historical data
    historical_data = historical_data[~((historical_data['month'] == 2) & (historical_data['day'] == 29))]


    # daily_stats = historical_data.groupby(['month', 'day'])['correctedSCA'].agg(['mean', 'min', 'max'])
    daily_stats = historical_data.groupby(['month', 'day'])['correctedSCA'].agg(
                        mean='mean',
                        min='min',
                        max='max',
                        perc5=lambda x: np.percentile(x.dropna(), 5),
                        perc25=lambda x: np.percentile(x.dropna(), 25),
                        perc50=lambda x: np.percentile(x.dropna(), 50),
                        perc75=lambda x: np.percentile(x.dropna(), 75),
                        perc95=lambda x: np.percentile(x.dropna(), 95)
                        )


    # Convert (month, day) to datetime using a dummy non-leap year for consistency  
    dummy_year = int(date_start.split('-')[0])
    new_index = []
    
    for month, day in daily_stats.index:
        year = dummy_year if month >= 10 else dummy_year + 1
        try:
            new_index.append(pd.Timestamp(year=year, month=month, day=day))
        except ValueError:
            # Skip invalid dates like Feb 29 if it slipped through
            continue
    
    daily_stats.index = pd.DatetimeIndex(new_index)


    # Filter from October 1 to date_end
    oct_start = pd.to_datetime(f"{dummy_year}-10-01")
    end_dt = pd.to_datetime(date_end)
    
    # If end date is earlier in the year (e.g. April), it means it belongs to the *next* calendar year in water year logic
    if end_dt.month < 10:
        end_dt = end_dt.replace(year=dummy_year + 1)
    else:
        end_dt = end_dt.replace(year=dummy_year)
    
    # Keep only the range October 1 to desired end date
    daily_stats = daily_stats[(daily_stats.index >= oct_start) & (daily_stats.index <= end_dt)]

    daily_stats.sort_index(inplace=True)

    
    # Remove Feb 29 from observation period if present
    sca = sca[~((sca.index.month == 2) & (sca.index.day == 29))]

    # Create DataFrame for plotting
    newdf = pd.DataFrame(index=sca.index)
    newdf['current SCA'] = sca.values
    newdf['mean SCA'] = daily_stats['mean'].iloc[:nr_days+1].values
    newdf['min SCA'] = daily_stats['min'].iloc[:nr_days+1].values
    newdf['max SCA'] = daily_stats['max'].iloc[:nr_days+1].values

    newdf['perc5 SCA'] = daily_stats['perc5'].iloc[:nr_days+1].values
    newdf['perc25 SCA'] = daily_stats['perc25'].iloc[:nr_days+1].values
    newdf['perc50 SCA'] = daily_stats['perc50'].iloc[:nr_days+1].values
    newdf['perc75 SCA'] = daily_stats['perc75'].iloc[:nr_days+1].values
    newdf['perc95 SCA'] = daily_stats['perc95'].iloc[:nr_days+1].values
    
    newdf.sort_index(inplace=True)


    # Plot
    plt.figure(figsize=(10, 5))
    plt.fill_between(sca.index, newdf['perc25 SCA'], newdf['perc75 SCA'], 
                     color='gray', alpha=0.2, label='25-75% range')


    # Fill lower tail: 5th to 25th percentile
    plt.fill_between(sca.index, newdf['perc5 SCA'], newdf['perc25 SCA'],
                     color='orangered', alpha=0.2, label='5-95% range')
    
    # Fill upper tail: 75th to 95th percentile
    plt.fill_between(sca.index, newdf['perc75 SCA'], newdf['perc95 SCA'],
                     color='orangered', alpha=0.2)

    plt.plot(sca.index, newdf['perc50 SCA'], label='50th percentile', color='gray', linestyle='-')
    plt.plot(sca.index, sca, color='orangered', linewidth=2, label='Current SCA')

    # plt.plot(sca.index, newdf['mean SCA'], label='Mean SCA', linestyle='dashed')

    plt.ylim([0, 100])
    plt.xticks(rotation=30)
    plt.ylabel('SCA [%]')
    plt.legend()
    plt.grid()
    
    # Save plot
    plt.savefig(os.path.join(work_folder, f"snow_bulletin_{date_start}_{date_end}_{suffix}.png"), bbox_inches="tight")

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
    mask_values = [0, 3, 5, 205, 255]
    
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
    bn = os.path.basename(snowMap_fileNameList[0])
    if 'modis' in bn or "complete" in bn:   
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
    elif 'vnp10a1f' in bn:   
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
                       shp_fileName=None, window=2, mode="yearly", 
                       save=True, ow=False):
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
    save: bool, optional
        Whether to save as GeoTIFF or not the snow cover duration maps.
    ow: bool, optional
        Whether to overwrite or not the output GeoTiff maps.

    Returns:
    --------
    dict
        Dictionary where keys are either **seasons (e.g., "2020-2021")** or **months ("YYYY-MM")** 
        and values are SCD datasets.
    """
    
    today = dt.today()
    
    # Extract dates from filenames
    dates = sorted([dateFromFileName(f) for f in snowMap_fileNameList])
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    # Open the raster to get CRS
    img = gdal.Open(snowMap_fileNameList[0])
    if img is None:
        raise ValueError(f"Could not open raster file: {snowMap_fileNameList[0]}")

    prj_img = img.GetProjection()
    img = None  # Close file to free memory
        
    results = {}

    # --- Yearly (Seasonal) Processing ---
    if mode == "yearly":
        for year in range(dates[0].year, dates[-1].year):  
            
            if year == today.year-1:
                print("Skipping current year")
                continue
            
            date_start = pd.Timestamp(year=year, month=10, day=1).date()  # 1st October
            date_end = pd.Timestamp(year=year + 1, month=9, day=30).date()  # 30th September
            
            outname = os.path.join(outdir, f'scd_{year}_{year+1}.tif')
            
            if not os.path.exists(outname) or ow:
                
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
                    
                    if save:
                        results[f"{year}-{year+1}"] = results[f"{year}-{year+1}"].rio.write_crs(prj_img, inplace=True)
                        results[f"{year}-{year+1}"].rio.to_raster(outname)
                else:
                    print(f"Skipping season {year}-{year+1} (too many missing days: {missing_days})")
            
            else:
                results[f"{year}-{year+1}"] = xr.open_dataset(outname)['band_data']
                
    
    # --- Monthly Processing ---
    elif mode == "monthly":
        for date in pd.date_range(start=dates[0], end=dates[-1], freq="MS"):  # MS = Month Start
        
            if date.year == today.year and date.month == today.month:
                print("Skipping current month")
                continue
        
            date_start = pd.Timestamp(year=date.year, month=date.month, day=1).date()
            date_end = pd.Timestamp(year=date.year, month=date.month, 
                                    day=pd.Period(date, freq='D').days_in_month).date()
            
            outname = os.path.join(outdir, f"{date.strftime('%Y-%m')}.tif")
            
            if not os.path.exists(outname) or ow:
    
    
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
                    
                    if save:
                        results[f"{date.strftime('%Y-%m')}"] = results[f"{date.strftime('%Y-%m')}"].rio.write_crs(prj_img, inplace=True)
                        results[f"{date.strftime('%Y-%m')}"].rio.to_raster(outname)
                        
                else:
                    print(f"Skipping month {date.strftime('%Y-%m')} (too many missing days: {missing_days})")
                    
            else:
                results[f"{date.strftime('%Y-%m')}"] = xr.open_dataset(outname)['band_data']
    
    # --- Trimester Processing ---
    elif mode == "trimester":
        for date in pd.date_range(start=dates[0], end=dates[-1], freq="QS"):  # quarter starts: Jan, Apr, Jul, Oct
            today_ts = pd.Timestamp.today().normalize()
            q_start = pd.Timestamp(date.year, date.month, 1).normalize()
            q_end = (q_start + pd.offsets.QuarterEnd()).normalize()
    
            # Skip only the trimester that includes 'today'
            if q_start <= today_ts <= q_end:
                print("Skipping current trimester")
                continue
    
            tnum = ((q_start.month - 1) // 3) + 1
            outname = os.path.join(outdir, f"{q_start.year}-T{tnum}.tif")
    
            if not os.path.exists(outname) or ow:
                trimester_files = [f for f in snowMap_fileNameList
                                   if q_start.date() <= dateFromFileName(f) <= q_end.date()]
    
                full_time_range = pd.date_range(start=q_start.date(), end=q_end.date(), freq="D")
                available_dates = set(dateFromFileName(f) for f in trimester_files)
                missing_days = len(full_time_range) - len(available_dates)
    
                if missing_days <= max_missing_days:
                    print(f"Processing trimester {q_start.year}-T{tnum} (missing days: {missing_days})")
                    results[f"{q_start.year}-T{tnum}"] = get_scd(
                        trimester_files, q_start.strftime('%d-%m-%Y'), q_end.strftime('%d-%m-%Y'),
                        shp_fileName, window
                    )
                    if save:
                        results[f"{q_start.year}-T{tnum}"] = results[f"{q_start.year}-T{tnum}"].rio.write_crs(prj_img, inplace=True)
                        results[f"{q_start.year}-T{tnum}"].rio.to_raster(outname)
                else:
                    print(f"Skipping trimester {q_start.year}-T{tnum} (too many missing days: {missing_days})")
            else:
                results[f"{q_start.year}-T{tnum}"] = xr.open_dataset(outname)['band_data']
    
    
    
    else:
        raise ValueError("Invalid mode! Use 'yearly', 'monthly' or 'trimester'.")
            
        
    return results



def compute_scd_anomaly(results, target_year):
    # Extract the target map
    target_scd = results[target_year]
    
    # Compute the mean of all other years
    other_years = [year for year in results if year != target_year]
    all_others = xr.concat([results[year] for year in other_years], dim='year')
    mean_scd = all_others.mean(dim='year')
    
    # Compute anomaly
    anomaly = target_scd - mean_scd
    
    plt.figure()
    plt.imshow(np.squeeze(anomaly.values), vmin=-1, vmax=1)
    plt.colorbar()
    return anomaly



def monthly_anomaly_scd(results, outdir):
    # Organize maps by month (e.g., "10" → [Oct_2021, Oct_2022, ...])
    monthly_data = defaultdict(list)
    for key, da in results.items():
        year_month = pd.to_datetime(key)
        month_str = f"{year_month.month:02d}"
        monthly_data[month_str].append((year_month.year, da))
    
    monthly_anomalies = {}
    
    # Assuming you want anomalies for the latest year in each month group
    for month_str, values in monthly_data.items():
        # Sort by year
        values = sorted(values, key=lambda x: x[0])
        
        # Separate latest from historical
        *historical, latest = values
        historical_maps = xr.concat([da for _, da in historical], dim="year")
        mean_map = historical_maps.mean(dim="year")
        
        latest_year, latest_map = latest
        anomaly = latest_map - mean_map
        monthly_anomalies[month_str] = (f"{latest_year}-{month_str}", anomaly)
    
    
    # Sort by actual datetime for proper chronological order
    ordered_items = sorted(monthly_anomalies.items(), key=lambda kv: pd.to_datetime(kv[1][0]))
    
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(16, 10))
    axes = axes.flatten()
    
    for i, (month_str, (label, da)) in enumerate(ordered_items):
        da.plot(ax=axes[i], cmap="RdBu", vmin=-1, vmax=1, center=0, cbar_kwargs={'label': 'SCD Anomaly'})
        axes[i].set_title(pd.to_datetime(label).strftime("%B %Y"))
        axes[i].tick_params(axis='x', labelrotation=20)
    
    # Hide any unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    
    plt.suptitle("Monthly Snow Cover Duration Anomalies", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(outdir, 'SCD_anomalies.png'))
    
    

def calculate_trends(csv_name):
    """
    Calculate seasonal and hydrological year snow cover area (SCA) trends
    and the Mann-Kendall trend test.


    Parameters
    ----------
    csv_name : str
        Path to the CSV file containing daily SCA data. The CSV must have a date column 
        as index (or in column `Unnamed: 0`) and a column named 'SCA'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
            - 'Nov-Dec' : mean SCA for November–December
            - 'Jan-Feb' : mean SCA for January–February (shifted to previous year)
            - 'Mar-Apr' : mean SCA for March–April (shifted to previous year)
            - 'Mean(Nov-Apr)' : mean SCA for the full hydrological year

        The index shows hydrological years in `YYYY/YYYY` format.
    """
    # Load and prepare
    df = pd.read_csv(csv_name)
    df.index = pd.to_datetime(df['Unnamed: 0'], format='%Y-%m-%d')
    df.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Seasonal means
    df_novdec = df[df.index.month.isin([11, 12])].copy()
    df_novdec['year'] = df_novdec.index.year
    df_novdec_mean = df_novdec.groupby('year').mean()
    
    df_janfeb = df[df.index.month.isin([1, 2])].copy()
    df_janfeb['year'] = df_janfeb.index.year
    df_janfeb_mean = df_janfeb.groupby('year').mean()
    
    df_marapr = df[df.index.month.isin([3, 4])].copy()
    df_marapr['year'] = df_marapr.index.year
    df_marapr_mean = df_marapr.groupby('year').mean()
    
    # Shift indices
    df_janfeb_mean.index = df_janfeb_mean.index - 1
    df_marapr_mean.index = df_marapr_mean.index - 1

    # Combine to final DataFrame
    newdf = pd.DataFrame({
        'Nov-Dec': df_novdec_mean['SCA'],
        'Jan-Feb': df_janfeb_mean['SCA'],
        'Mar-Apr': df_marapr_mean['SCA']
    })

    # Add hydrological mean
    newdf['Mean(Nov-Apr)'] = newdf[['Nov-Dec', 'Jan-Feb', 'Mar-Apr']].mean(axis=1)

    # Add hydrological year label YYYY/YYYY
    newdf.index = [f"{y}/{y+1}" for y in newdf.index]

    # Trend test
    mk.original_test(newdf['Mean(Nov-Apr)'].dropna())

    return newdf

