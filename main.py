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
shp_fileName = r'/home/vpremier/Documents/git/sca4provbz/input/ST_shape/SouthTyrol.shp'

satellite = 'MODIS'

if satellite == 'MODIS':
    # for EURAC SNOW dataset
    csv_path = r'/home/vpremier/Documents/git/sca4provbz/input/ST_shape_alps.csv'
    pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/ST'
    snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + '*/EURAC_SNOW*.tif')

elif satellite == 'VIIRS':
    # for VNP10A1F dataset
    csv_path = r'/home/vpremier/Documents/git/sca4provbz/input/ST_viirs_shape.csv'
    pathToDataFolder = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol'
    snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'VNP10A1F*.tif')


work_folder = os.getcwd()

start_time = time.time()

# update a csv file containing the statistics
sca_data = updateCSV(csv_path, snowMap_fileNameList, shp_fileName=shp_fileName)

print("--- %s seconds ---" % (time.time() - start_time))


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

# date_start = '2023-10-01'
# date_end = '2024-06-29'

print(f"date_start: {date_start}, date_end: {date_end}")

statistics = snow_bullettin(csv_path, date_start, date_end, work_folder)

sorted_files = sorted(snowMap_fileNameList, key=dateFromFileName)
# Open a text file to save the results
with open(r"/home/vpremier/Documents/git/sca4provbz/ST_to_reprocess.txt", "a") as outfile:
    ref = None  # Reference bounds

    for i, f in enumerate(sorted_files):
        print(f)
        ds = rioxarray.open_rasterio(f)
        
        # Get bounds
        xmin, ymin, xmax, ymax = ds.rio.bounds()
        bounds = (xmin, ymin, xmax, ymax)
        
        # Set the first file as the reference
        if i == 0:
            ref = bounds
        
        # If bounds differ, print and write to file
        if bounds != ref:
            print(f"{f} has different bounds: {bounds}")
            outfile.write(f"{f}")
            
            
sss            

"""
SNOW COVER DURATION (SCD)
"""
outdir_y = r'/home/vpremier/Documents/git/sca4provbz/scd/VIIRS/yearly'
outdir_m = r'/home/vpremier/Documents/git/sca4provbz/scd/VIIRS/monthly'


scd_yearly = get_scd_statistics(snowMap_fileNameList, outdir_y, max_missing_days=30, 
                                shp_fileName=shp_fileName, window=2, mode="yearly")

scd_monthly = get_scd_statistics(snowMap_fileNameList, outdir_m, max_missing_days=2, 
                                shp_fileName=shp_fileName, window=2, mode="monthly")


