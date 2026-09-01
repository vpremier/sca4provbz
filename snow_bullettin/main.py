#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 11:30:59 2025

@author: vpremier
"""

import glob
import os
import time
from datetime import datetime as dt

from utils import (
    updateCSV, 
    snow_bullettin,
    get_scd_statistics,
    monthly_anomaly_scd,
    )

# shape with the area of interest
shp_fileName = r'/mnt/CEPH_PROJECTS/PROSNOW/SCA4provbz/shapefile/ST_shape/SouthTyrol.shp'
wd = "/mnt/CEPH_PROJECTS/PROSNOW/SCA4provbz/results"

snow_bullettin_dir = os.path.join(wd, "snow_bullettin")

suffix = 'modis'

if suffix == 'modis':
    
    # for EURAC SNOW dataset
    csv_path = os.path.join(snow_bullettin_dir, "ST_modis.csv")
    pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/ST/*'
    # pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/alps/*'
    

elif suffix == 'vnp10a1f':
    
    # for VNP10A1F dataset
    csv_path = os.path.join(snow_bullettin_dir, "ST_viirs.csv")
    pathToDataFolder = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol'



snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'EURAC_SNOW_MERGE.alps.south-tyrol.20*.tif')
# snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'EURAC_SNOW.alps.complete.20*.tif')




start_time = time.time()

# update a csv file containing the statistics
sca_data = updateCSV(csv_path, snowMap_fileNameList, shp_fileName=shp_fileName)

print("--- %s seconds ---" % (time.time() - start_time))


"""
SNOW BULLETTIN: Create daily plot with current SCA, mean SCA, min and max 
                (or percentiles: check the function!)
                
    To be run every 15 days
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


print(f"date_start: {date_start}, date_end: {date_end}")

statistics = snow_bullettin(csv_path, date_start, date_end, snow_bullettin_dir, suffix)


                     
          

"""
SNOW COVER DURATION (SCD)
    To be run every month 
"""
outdir_y = os.path.join(wd, f"scd/{suffix}/yearly")
outdir_m = os.path.join(wd, f"scd/{suffix}/monthly")
outdir_t = os.path.join(wd, f"scd/{suffix}/trimester") 



scd_trimester = get_scd_statistics(snowMap_fileNameList, outdir_t, max_missing_days=20, 
                                shp_fileName=shp_fileName, window=2, mode="trimester")


scd_yearly = get_scd_statistics(snowMap_fileNameList, outdir_y, max_missing_days=71, 
                                shp_fileName=shp_fileName, window=2, mode="yearly")

scd_monthly = get_scd_statistics(snowMap_fileNameList, outdir_m, max_missing_days=15, 
                                shp_fileName=shp_fileName, window=2, mode="monthly")


# SCD anomalies
monthly_anomaly_scd(scd_monthly, snow_bullettin_dir)

