'''
This compares to CCMP versions.

The output is all in  map  form.

Not as well constructed as the  CCMP satllite comparisons, but could be
used as a basis for a more complete comparison

'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import calendar
import copy
import os
import cProfile

import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/job_CCMP/python/binned_stats/')
sys.path.append('C:/job_CCMP/python/map_stats/')
sys.path.append('C:/job_CCMP/python/CCMP/')
sys.path.append('C:/job_CCMP/python/eval_era5/')
sys.path.append('C:/job_CCMP/python/era5/')
sys.path.append('C:/job_CCMP/python/buoys/')

sys.path.append('C:/python_packages_cam/')
sys.path.append('../')

from ccmp import read_ccmp_native,load_ccmp_daily_arrays
from era5 import convert_to_rss_grid
from time_interpolate_synoptic_maps import time_interpolate_synoptic_maps

#plotting
from global_map import global_map,global_map_w_zonal_mean
from plt_zonal_mean_array import plt_zonal_mean_array
from plotting import plot_scat

#RSS statistic accumulators
from rss_stats.binned_stats import BinnedStat
from rss_stats.hist_1d import Hist1D
from rss_stats.hist_2d import Hist2D
from rss_stats.map_stats import MapStat
from rss_maps import zonal_mean, global_mean

# one reader to read them all
from rss_sat_readers import read_rss_satellite_daily_xr
from rss_sat_readers import apply_scat_flags

def compare_ccmp_versions_month(year=2013,month = 1, version = 'era5_current_corrected_adj2'):
    u_ccmp_era5_adj_minus_ccmp_era5_stats =  MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = -179.875,max_lon=179.875)

    v_ccmp_era5_adj_minus_ccmp_era5_stats =  MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = -179.875,max_lon=179.875)

    w_ccmp_era5_adj_minus_ccmp_era5_stats =  MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = -179.875,max_lon=179.875)

    ccmp_nobs_stats =                     MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = -179.875,max_lon=179.875)


    plt_path = 'C:/job_CCMP/python/eval_ccmp_versions/plots/'+version+'/'
    os.makedirs(plt_path,exist_ok = True)

    num_days = calendar.monthrange(year,month)[1]
    first_time_step = True
    for day in range(1,num_days+1):
        for time_step in [0,1,2,3]:
            map_index = time_step + (day-1)*4

            date_str = f"{year:04}_{month:02}_{day:02}_{6*time_step:02}"
            print(date_str)
            
            # read in CCMP with ERA5 background
            version = 'era5'
            version_path = f"P:/CCMP/MEaSUREs/{version}/v2.0/native/"
            ccmp_era5 = read_ccmp_native(year=year,month=month,day=day,time_step = time_step,path = version_path,return_as_xr  = True)

            # read in CCMP with ERA5 background, corrected for OSCAR currents, scaled
            version = 'era5_current_corrected_adj2'
            version_path = f"P:/CCMP/MEaSUREs/{version}/v2.0/native/"
            ccmp_era5_adj = read_ccmp_native(year=year,month=month,day=day,time_step = time_step,path = version_path,return_as_xr  = True)

            #calculate wind speeds for components
            ccmp_era5_w    = np.sqrt(np.square(ccmp_era5['u10']) + np.square(ccmp_era5['v10'])).values
            ccmp_era5_adj_w = np.sqrt(np.square(ccmp_era5_adj['u10']) + np.square(ccmp_era5_adj['v10'])).values

            #calculate difference between ccmp versions
            u_ccmp_era5_adj_minus_ccmp_era5 = ccmp_era5_adj['u10'].values - ccmp_era5['u10'].values
            v_ccmp_era5_adj_minus_ccmp_era5 = ccmp_era5_adj['v10'].values - ccmp_era5['v10'].values
            w_ccmp_era5_adj_minus_ccmp_era5 = ccmp_era5_adj_w - ccmp_era5_w

            u_ccmp_era5_adj_minus_ccmp_era5_stats.add_map(u_ccmp_era5_adj_minus_ccmp_era5)
            v_ccmp_era5_adj_minus_ccmp_era5_stats.add_map(v_ccmp_era5_adj_minus_ccmp_era5)
            w_ccmp_era5_adj_minus_ccmp_era5_stats.add_map(w_ccmp_era5_adj_minus_ccmp_era5)

            #accumulate nobs used in CCMP
            ccmp_num_obs = ccmp_era5_adj['nobs'].values
            ccmp_nobs_stats.add_map(ccmp_num_obs)

            if first_time_step:

                u_ccmp_era5_adj_minus_ccmp_era5_mean = u_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
                ccmp_diff_map,ax = global_map_w_zonal_mean(u_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                                    vmin=-0.5, zmin=-0.5, 
                                                    vmax=0.5, zmax=0.5, 
                                                    plt_colorbar=True,
                                                    cmap = 'BrBG',
                                                    title='U10, CCMP ERA5_CC - CCMP_ERA5')
                png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.u10.png'
                ccmp_diff_map.savefig(png_file)

                v_ccmp_era5_adj_minus_ccmp_era5_mean = v_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
                ccmp_diff_map,ax =global_map_w_zonal_mean(v_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                                    vmin=-0.5, zmin=-0.5, 
                                                    vmax=0.5, zmax=0.5, 
                                                    plt_colorbar=True,
                                                    cmap = 'BrBG',
                                                    title='V10, CCMP ERA5_CC - CCMP_ERA5')
                png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.v10.png'
                ccmp_diff_map.savefig(png_file)

                w_ccmp_era5_adj_minus_ccmp_era5_mean = w_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
                ccmp_diff_map,ax = global_map_w_zonal_mean(w_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                                    vmin=-0.5, zmin=-0.5, 
                                                    vmax=0.5, zmax=0.5, 
                                                    plt_colorbar=True,
                                                    cmap = 'BrBG',
                                                    title='W10, CCMP ERA5_CC - CCMP_ERA5')
                png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.w10.png'
                ccmp_diff_map.savefig(png_file)

                num_map_to_plot = ccmp_nobs_stats.tot
                max_obs = np.nanmax(num_map_to_plot)
                temp = np.power(10,np.floor(np.log10(max_obs)))
                vmax = round(max_obs/temp)*temp

                num_map,ax = global_map(num_map_to_plot, 
                                                    vmin=0,
                                                    vmax=vmax, 
                                                    plt_colorbar=True,
                                                    title='Number of Obs, CCMP_ERA5_CC')
                png_file = plt_path+'num_obs_ccmp_era5_adj_'+date_str+'.png'
                num_map.savefig(png_file)

                first_time_step = False


    date_str = f"{year:04}_{month:02}"



    u_ccmp_era5_adj_minus_ccmp_era5_mean = u_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
    ccmp_diff_map,ax = global_map_w_zonal_mean(u_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                        vmin=-0.5, zmin=-0.5, 
                                        vmax=0.5, zmax=0.5,
                                        plt_colorbar=True,
                                        cmap = 'BrBG',
                                        title='U10, CCMP ERA5_CC - CCMP_ERA5')
    png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.u10.png'
    ccmp_diff_map.savefig(png_file)

    v_ccmp_era5_adj_minus_ccmp_era5_mean = v_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
    ccmp_diff_map,ax = global_map_w_zonal_mean(v_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                        vmin=-0.5, zmin=-0.5, 
                                        vmax=0.5, zmax=0.5, 
                                        plt_colorbar=True,
                                        cmap = 'BrBG',
                                        title='V10, CCMP ERA5_CC - CCMP_ERA5')
    png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.v10.png'
    ccmp_diff_map.savefig(png_file)

    w_ccmp_era5_adj_minus_ccmp_era5_mean = w_ccmp_era5_adj_minus_ccmp_era5_stats.mean()
    ccmp_diff_map,ax = global_map_w_zonal_mean(w_ccmp_era5_adj_minus_ccmp_era5_mean,  
                                        vmin=-0.5, zmin=-0.5, 
                                        vmax=0.5, zmax=0.5, 
                                        plt_colorbar=True,
                                        cmap = 'BrBG',
                                        title='W10, CCMP ERA5_CC - CCMP_ERA5')
    png_file = plt_path+'ccmp_era5_adj_minus_ccmp_era5_'+date_str+'.w10.png'
    ccmp_diff_map.savefig(png_file)

    num_map_to_plot = ccmp_nobs_stats.tot
    max_obs = np.nanmax(num_map_to_plot)
    temp = np.power(10,np.floor(np.log10(max_obs)))
    vmax = round(max_obs/temp)*temp

    num_map,ax = global_map_w_zonal_mean(num_map_to_plot, 
                                        vmin=0,
                                        vmax=vmax, 
                                        plt_colorbar=True,
                                        title='Number of Obs, CCMP_ERA5_CC')
    png_file = plt_path+'num_obs_ccmp_era5_adj_'+date_str+'.png'
    num_map.savefig(png_file)

if __name__ == '__main__':
    version = 'era5_current_corrected_adj2'
    cProfile.run('compare_ccmp_versions_month(year=2013,month = 1, version = version)')