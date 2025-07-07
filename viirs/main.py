#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 14:09:31 2023

@author: vpremier
"""
import os
import glob
import geopandas as gpd
from datetime import datetime as dt
import earthaccess
from scf import get_scf_viirs


""" 
This scipt is used to download VNP10A1F snow cover data. The product is daily
and gap-filled.

DOWNLOAD
The data are downloaded from https://www.earthdata.nasa.gov/
(see script data_download.py)

RESAMPLING/REPROJECTION and SCF retrieval
The data are reprojected to the AOI. You can specify either a tif file
with the target extent and crs or they can be set manually. 
Please, specify the target resolution (defult is 500 m).
SCF is calculated from the NDSI.
"""

# credentials 
# username = 'vpremier'
# password = 'Eurac_2022'
auth = earthaccess.login(strategy='netrc')
# set your credentials: https://earthaccess.readthedocs.io/en/latest/howto/authenticate/

# shapefile with the AOI
shp = r'/home/vpremier/Documents/git/sca4provbz/shapefile/ST_shape/SouthTyrol.shp'

#read the shapefile
gdf = gpd.read_file(shp)

# directory where you want to store the raw data
download_dir = r'/mnt/CEPH_PROJECTS/PROSNOW/raw_data/VIIRS/VNP10A1F'

# list of downloaded files
fileList = glob.glob(download_dir + os.sep + 'VNP10A1F.A20*.h5')
fileList.sort()
       
# period for the download: from the last downloaded date until today
last_date = os.path.basename(fileList[-1]).split('.')[1][1:] 
date_start = dt.strptime(last_date, '%Y%j').strftime('%Y-%m-%d')
# date_start = '2025-01-01'
# date_end = '2024-10-23'
date_end = dt.today().strftime('%Y-%m-%d')


# convert crs (otherwise may result in an error)
if not gdf.crs == 'EPSG:4326':
    gdf = gdf.to_crs('EPSG:4326')  
    
# Extract the Bounding Box Coordinates
bounds = tuple(gdf.total_bounds)

    
results = earthaccess.search_data(
    short_name='VNP10A1F',
    version='2',
    bounding_box=bounds,
    temporal=(date_start, date_end),
    count=1000
)

# download the data
files = earthaccess.download(results, download_dir)



# list of downloaded files
fileList = glob.glob(download_dir + os.sep + 'VNP10A1F.A20*.h5')

# output directory where the SCF maps are saved
outdir = r'/mnt/CEPH_PROJECTS/PROSNOW/4.results/VNP10A1F_SouthTyrol'

epsg_target = '32632'

# get bounding box
# shapefile = gpd.read_file(shp)
# shp_rpj = shapefile.to_crs(crs=epsg_target)
# bbox = list(shp_rpj.bounds.iloc[0])

# keep the extent asked by the province
extent_target = [597972, 5117987, 767972, 5221987]

# create the SCF maps
get_scf_viirs(fileList, outdir, res = 250, extent_target = extent_target, 
                epsg_target = epsg_target, ow = False)
