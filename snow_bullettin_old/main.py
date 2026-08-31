#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 11:30:59 2025

@author: vpremier
"""

import glob
import os
import numpy as np
import time
from datetime import datetime as dt

from utils import *

# shape with the area of interest
shp_fileName = r'/mnt/CEPH_PROJECTS/PROSNOW/SCA4provbz/shapefile/ST_shape/SouthTyrol.shp'
# shp_fileName = r'/mnt/CEPH_PROJECTS/PROSNOW/1.data/shapefiles/Cortina/Cortina_buffer.geojson'
# shp_fileName = r'/mnt/CEPH_PROJECTS/PROSNOW/1.data/Casera/area_Casera.shp'

# rilanciare forse per il problema no data??

suffix = 'modis'

if suffix == 'modis':
    # for EURAC SNOW dataset
    csv_path = r'/mnt/CEPH_PROJECTS/PROSNOW/SCA4provbz/results/snow_bullettin/ST_modis.csv'
    pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/ST/*'
    # csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin/TAA_modis.csv'
    # csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin/Cortina_modis.csv'
    # pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/alps/*'
    

elif suffix == 'vnp10a1f':
    # for VNP10A1F dataset
    csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin/ST_viirs.csv'
    # csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin/Cortina_viirs.csv'
    # pathToDataFolder = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_Cortina'
    pathToDataFolder = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol'
    # csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin/TAA_viirs.csv'
    # pathToDataFolder = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_TAA'

snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'EURAC_SNOW_MERGE.alps.south-tyrol.20*.tif')
# snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'EURAC_SNOW.alps.complete.20*.tif')


outdir = r'/home/vpremier/Documents/git/sca4provbz/results/snow_bullettin'

start_time = time.time()

# update a csv file containing the statistics
sca_data = updateCSV(csv_path, snowMap_fileNameList, shp_fileName=shp_fileName)

print("--- %s seconds ---" % (time.time() - start_time))

ddd
"""
SNOW BULLETTIN: Create daily plot with current SCA, mean SCA, min and max 
                (or percentiles: check the function!)
"""


# Get last available date
date_end = sca_data.index[-1]

# Determine date_start as the last October 1st
if date_end.month >= 10:
    date_start = dt(date_end.year, 10, 1)  # If date_end is Oct or later, use this year's October 1st
else:
    date_start = dt(date_end.year - 1, 10, 1)  # If before October, use last year's October 1st

# Convert to string format
date_start = date_start.strftime('%Y-%m-%d')
date_end = date_end.strftime('%Y-%m-%d')

# date_start = '2024-10-01'
# date_end = '2025-09-30'

print(f"date_start: {date_start}, date_end: {date_end}")

statistics = snow_bullettin(csv_path, date_start, date_end, outdir, suffix)


            
            
          

"""
SNOW COVER DURATION (SCD)
"""
outdir_y = f'/home/vpremier/Documents/git/sca4provbz/results/scd/{suffix}/yearly'
outdir_m = f'/home/vpremier/Documents/git/sca4provbz/results/scd/{suffix}/monthly'
outdir_t = f'/home/vpremier/Documents/git/sca4provbz/results/scd/{suffix}/trimester'


# outdir_y = f'/home/vpremier/Documents/git/sca4provbz/results/scd_Casera'


# scd_trimester = get_scd_statistics(snowMap_fileNameList, outdir_t, max_missing_days=20, 
#                                 shp_fileName=shp_fileName, window=2, mode="trimester")


scd_yearly = get_scd_statistics(snowMap_fileNameList, outdir_y, max_missing_days=71, 
                                shp_fileName=shp_fileName, window=2, mode="yearly")

scd_monthly = get_scd_statistics(snowMap_fileNameList, outdir_m, max_missing_days=15, 
                                shp_fileName=shp_fileName, window=2, mode="monthly")


# # SCD anomalies
monthly_anomaly_scd(scd_monthly, outdir)

aa


# Extract date from filename using regex
dates = []
for f in snowMap_fileNameList:
    match = re.search(r"\.(\d{8})T\d{6}\.", f)
    if match:
        dates.append(dt.strptime(match.group(1), "%Y%m%d"))

# Create DataFrame
df = pd.DataFrame(dates, columns=["date"]).sort_values("date").reset_index(drop=True)

# Generate complete date range (from min to max date)
full_range = pd.date_range(df["date"].min(), df["date"].max(), freq="D")

# Find missing dates
missing = full_range.difference(df["date"])

print("DataFrame of available dates:")
print(df)

print("\nMissing dates:")
print(missing)