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

shp_fileName = r'/home/vpremier/Documents/git/sca4provbz/input/ST_box/SouthTyrol.shp'
csv_path = r'/home/vpremier/Documents/git/sca4provbz/input/ST_box_alps.csv'
pathToDataFolder = r'/mnt/CEPH_PRODUCTS/EURAC_SNOW/MODIS/alps'

snowMap_fileNameList = glob.glob(pathToDataFolder + os.sep + '*/EURAC_SNOW**.tif')

work_folder = os.getcwd()


start_time = time.time()

# funziona con la codifica dei dati Eurac, modificare per le VIIRS!
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
    date_start = dt(today.year, 10, 1)  # If date_end is Oct or later, use this year's October 1st
else:
    date_start = dt(today.year - 1, 10, 1)  # If before October, use last year's October 1st

# Convert to string format
date_start = date_start.strftime('%Y-%m-%d')
date_end = date_end.strftime('%Y-%m-%d')

# date_start = '2023-10-01'
# date_end = '2024-06-29'

print(f"date_start: {date_start}, date_end: {date_end}")

statistics = snow_bullettin(csv_path, date_start, date_end, work_folder)



"""
SNOW COVER DURATION (SCD)
"""


dateStart = '20020101T120000'
    
fileList = [f for f in snowMap_fileNameList if dateFromFileName(f) > dt.strptime(dateStart, '%Y%m%dT%H%S%M').date()]
date = [dateFromFileName(f) for f in fileList if dateFromFileName(f) > dt.strptime(dateStart, '%Y%m%dT%H%S%M').date()]


fileList = sorted(fileList)[120:]
date = sorted(date)[120:]


filt.cloud_temporalFilter(fileList, date, 3, 0, output_root, water=None, waterMask_fileName=None, confidence=False,tmax=3)

