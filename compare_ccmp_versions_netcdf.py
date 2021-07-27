'''
This is a short script to compare V2.1 and  V2.0 NRT CCMP to make sure winds in the file are the  same.
It is really a signle use piece of code, and  could probably be deleted soon.
I put a raise statement in the beginning to alert anyone that launches it
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/python_packages_cam/rss_stats/')
sys.path.append('C:/job_CCMP/python/CCMP/')
sys.path.append('C:/job_CCMP/python/eval_era5/')
sys.path.append('C:/job_CCMP/python/era5/')

from global_map import global_map
import logging
import datetime
from calendar import monthrange

raise ValueError("Why are you running this -- its purpose has been  served")


log_file = 'P:/CCMP_RT/MEaSUREs/flk/v2.0_vs_V2.1_compare/compare_v21_v20_log_.'+datetime.datetime.now().strftime("%Y_%m_%d.%H%M")+'.log'
logger = logging.getLogger('ccmp_compare')
file_hdlr = logging.FileHandler(log_file)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
file_hdlr.setFormatter(formatter)
file_hdlr.setLevel(logging.DEBUG)

logger.addHandler(file_hdlr )

screen_hdlr = logging.StreamHandler()
screen_hdlr.setLevel(logging.INFO)
screen_hdlr.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
logger.addHandler(screen_hdlr)

logger.setLevel(logging.DEBUG)
logger.info('Starting CCMP V2.1 vs CCMP v2.0 netcdf comparison')

year = 2020
for year in range(2015,2020):
    last_month = 12
    if year == 2020:
        last_month = 8
    for month in range(1,last_month+1):
        logger.info(f"Checking {month:02d}/{year:04d}")
        junk,num_days_in_month = monthrange(year,month)
        for day in range(1,num_days_in_month+1):
            print(f"{month:02d}/{day:02d}/{year:04d}")
            # read in CCMP with ERA5 background
            version = '2.0'
            version_path = f"P:/CCMP_RT/MEaSUREs/flk/v{version}/netcdf/"
            filename  = f"{version_path}Y{year:04d}/M{month:02d}/CCMP_RT_Wind_Analysis_{year:04d}{month:02d}{day:02d}_V0{version}_L3.0_RSS.nc"
            if year <= 2016:
               filename  = f"{version_path}Y{year:04d}/M{month:02d}/CCMP_RT_Wind_Analysis_{year:04d}{month:02d}{day:02d}_V02.1_L3.0_RSS.nc" 
            files_missing = False
            try:
                ccmp_v20 = xr.open_dataset(filename)
            except FileNotFoundError:
                logger.info(f"File Missing: {filename}")
                files_missing = True

            version = '2.1'
            version_path = f"P:/CCMP_RT/MEaSUREs/flk/v{version}/netcdf/"
            filename  = f"{version_path}Y{year:04d}/M{month:02d}/CCMP_RT_Wind_Analysis_{year:04d}{month:02d}{day:02d}_V0{version}_L3.0_RSS.nc"
            try:
                ccmp_v21 = xr.open_dataset(filename)
            except FileNotFoundError:
                logger.info(f"File Missing: {filename}")
                files_missing = True
            if files_missing:
                continue

            for timestep in [0,1,2,3]:
                # u_ccmp_v20 = np.zeros((720,1440))*np.nan
                # u_ccmp_v20[46:46+628,:] = ccmp_v20['uwnd'].values[timestep,:,:]

                # u_ccmp_v21 = np.zeros((720,1440))*np.nan
                # u_ccmp_v21[46:46+628,:] = ccmp_v21['uwnd'].values[timestep,:,:]

                # v_ccmp_v20 = np.zeros((720,1440))*np.nan
                # v_ccmp_v20[46:46+628,:] = ccmp_v20['vwnd'].values[timestep,:,:]

                # v_ccmp_v21 = np.zeros((720,1440))*np.nan
                # v_ccmp_v21[46:46+628,:] = ccmp_v21['vwnd'].values[timestep,:,:]

                #global_map(u_ccmp_v20,vmin=-20.0,vmax=20.0,title='CCMP V2.0')
                #global_map(u_ccmp_v21,vmin=-20.0,vmax=20.0,title='CCMP V2.1')

                #global_map(u_ccmp_v21-u_ccmp_v20,vmin=-2.0,vmax=2.0,title='CCMP V2.1 minus CCMP V2.0')
                date_string = f"{month:02d}/{day:02d}/{year:04d} {timestep*6:02d}Z"
                max_udiff = np.nanmax(np.abs(ccmp_v21['uwnd'].values[timestep,:,:]-
                                             ccmp_v20['uwnd'].values[timestep,:,:]))
                max_vdiff = np.nanmax(np.abs(ccmp_v21['vwnd'].values[timestep,:,:]-
                                             ccmp_v20['vwnd'].values[timestep,:,:]))

                if (max_vdiff > 0.0001) or (max_udiff > 0.0001):
                    logger.info(f"Compare Failed: {date_string} max udiff: {max_udiff}  max vdiff: {max_vdiff}")
print()
