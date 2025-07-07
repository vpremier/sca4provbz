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
import pandas as pd

from utils import *

"""
Updated script to produce the results for the Eurac climate monitoring indicators page
see https://www.eurac.edu/en/data-in-action/climate-change-monitoring/snow-coverage

Used inputs:
    - SCA derived from MODIS/ST considering the South Tyrol box (no shapefile applied)
    
    
Needed information:
    - csv file with the statistics current, min, max and mean SCA
    - csv file with the SCA for elevation belts (<1000, 1000-2000, 2000-3000, >3000)

"""

suffix = 'modis'

csv_path = r'/home/vpremier/Documents/git/sca4provbz/results/indicators/ST_box_modis.csv'
pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/ST/*'
    


snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + 'EURAC_SNOW*.tif')


outdir = r'/home/vpremier/Documents/git/sca4provbz/results/indicators'

start_time = time.time()

# update a csv file containing the statistics
sca_data = updateCSV(csv_path, snowMap_fileNameList)

print("--- %s seconds ---" % (time.time() - start_time))


"""
SNOW BULLETTIN: Create daily plot with current SCA, mean SCA, min and max 
                (or percentiles: check the function!)
"""


date_start = '2024-10-01'
date_end = '2025-06-25'

# Parse years 
year_start = pd.to_datetime(date_start).year
year_end = pd.to_datetime(date_end).year

statistics = snow_bullettin(csv_path, date_start, date_end, outdir, suffix)
stats_name = f'/home/vpremier/Documents/git/sca4provbz/results/indicators/statistics_{year_start}_{year_end}.csv'

statistics[['current SCA', 'min SCA', 'max SCA', 'mean SCA']].to_csv(stats_name)

            
            

# List of (mask raster path, output CSV path) tuples
mask_csv_pairs = [
    (r'/home/vpremier/Documents/git/sca4provbz/shapefile/dem_1000.tif',
     os.path.join(outdir,'SouthTyrol_1000.csv')),

    (r'/home/vpremier/Documents/git/sca4provbz/shapefile/dem_1000_2000.tif',
     os.path.join(outdir,'SouthTyrol_1000_2000.csv')),

    (r'/home/vpremier/Documents/git/sca4provbz/shapefile/dem_2000_3000.tif',
     os.path.join(outdir,'SouthTyrol_2000_3000.csv')),
    
    (r'/home/vpremier/Documents/git/sca4provbz/shapefile/dem_3000.tif',
     os.path.join(outdir,'SouthTyrol_3000.csv'))
]

# Loop through each pair and run updateCSV
for mask_raster_fileName, csv_output_path in mask_csv_pairs:
    print(f"Updating: {csv_output_path} with mask: {mask_raster_fileName}")
    start_time = time.time()

    updated_sca = updateCSV(
        csv_output_path,
        snowMap_fileNameList,
        mask_raster_fileName=mask_raster_fileName
    )

    elapsed = time.time() - start_time
    print(f"✅ Finished: {csv_output_path} --- {elapsed:.2f} seconds ---\n")



# compute trends
trend_1000 = calculate_trends(os.path.join(outdir,'SouthTyrol_1000.csv'))
trend_1000.to_csv(os.path.join(outdir,'trend_1000.csv'))

trend_1000_2000 = calculate_trends(os.path.join(outdir,'SouthTyrol_1000_2000.csv'))
trend_1000_2000.to_csv(os.path.join(outdir,'trend_1000_2000.csv'))

trend_2000_3000 = calculate_trends(os.path.join(outdir,'SouthTyrol_2000_3000.csv'))
trend_2000_3000.to_csv(os.path.join(outdir,'trend_2000_3000.csv'))

trend_3000 = calculate_trends(os.path.join(outdir,'SouthTyrol_3000.csv'))
trend_3000.to_csv(os.path.join(outdir,'trend_3000.csv'))


