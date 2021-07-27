import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import calendar
import copy
import os
import cmocean

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

#plotting
from global_map import global_map
from land_masks.land_fraction import read_land_fraction_1440_720

from plt_zonal_mean_array import plt_zonal_mean_array
from plotting import plot_scat

#RSS statistic accumulators

from rss_stats.hist_1d import Hist1D
from rss_stats.map_stats import MapStat
from rss_maps import zonal_mean, global_mean, lat_wt_map

#Time series analysis
from time_series_fits import harmonic_fit



start_year = 2000
end_year = 2018
num_years = 1+end_year-start_year

version = 'era5'
nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'

plot_trend_maps = False
plot_time_series = True

if plot_time_series:
    subversion = 'v2.0'
    nc_path = f"{nc_root}{version}_{subversion}/"
    nc_file_time_series = f"{nc_path}CCMP_time_series_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"
    time_series_v20 = xr.open_dataset(nc_file_time_series)

    subversion = 'current_corrected_adj2_v2.0'
    nc_path = f"{nc_root}{version}_{subversion}/"
    nc_file_time_series = f"{nc_path}CCMP_time_series_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"
    time_series_v20cca = xr.open_dataset(nc_file_time_series)

    print(time_series_v20.keys())
    yrs = time_series_v20['Years']
    diff_u = time_series_v20cca['u_time_series'].values - time_series_v20['u_time_series'].values
    diff_v = time_series_v20cca['v_time_series'].values - time_series_v20['v_time_series'].values
    diff_w = time_series_v20cca['w_time_series'].values - time_series_v20['w_time_series'].values

    fig,axs = plt.subplots(nrows=3,ncols=1,sharex=True)


    axs[0].plot(yrs,diff_u)
    axs[0].set_title('Global U Wind Difference, V2.0 CCA - V2.0')
    axs[1].plot(yrs,diff_v)
    axs[1].set_title('Global V Wind Difference, V2.0 CCA - V2.0')
    
    axs[2].plot(yrs,diff_w) 
    axs[2].set_title('Global W Wind Difference, V2.0 CCA - V2.0')
 
    plt.show()          
    print




if plot_trend_maps:
    subversion = 'v2.0'
    nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'
    nc_path = f"{nc_root}{version}_{subversion}/"   
    nc_file_trend_map = f"{nc_path}CCMP_trend_maps_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"

    trend_maps_v20 = xr.open_dataset(nc_file_trend_map)

    subversion = 'current_corrected_v2.0'
    nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'
    nc_path = f"{nc_root}{version}_{subversion}/"
    nc_file_trend_map = f"{nc_path}CCMP_trend_maps_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"

    trend_maps_v20cc = xr.open_dataset(nc_file_trend_map)

    subversion = 'current_corrected_adj2_v2.0'
    nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'
    nc_path = f"{nc_root}{version}_{subversion}/"
    nc_file_trend_map = f"{nc_path}CCMP_trend_maps_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"

    trend_maps_v20cca = xr.open_dataset(nc_file_trend_map)

    fig,ax = global_map(10.0/12.0*trend_maps_v20['u_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='U Trend, CCMP V2.0')
    fig,ax = global_map(10.0/12.0*trend_maps_v20cca['u_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='U TrendmCCMP V2.0 CC Adj2')

    diff_cca_20 = trend_maps_v20cca['u_trend'].values - trend_maps_v20['u_trend'].values
    fig,ax = global_map(10.0/12.0*diff_cca_20, vmin=-0.02, vmax=0.02, cmap='BrBG', plt_colorbar=True,title='U Trend, CCMP V2.0 CCA - CCMP V2.0')

    fig,ax = global_map(10.0*trend_maps_v20['v_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='V Trend, CCMP V2.0')
    fig,ax = global_map(10.0*trend_maps_v20cca['v_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='V TrendmCCMP V2.0 CC Adj2')

    diff_cca_20 = trend_maps_v20cca['v_trend'].values - trend_maps_v20['v_trend'].values
    fig,ax = global_map(10.0*diff_cca_20, vmin=-0.02, vmax=0.02, cmap='BrBG', plt_colorbar=True,title='V Trend, CCMP V2.0 CCA - CCMP V2.0')

    fig,ax = global_map(10.0*trend_maps_v20['w_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='W Trend, CCMP V2.0')
    fig,ax = global_map(10.0*trend_maps_v20cca['w_trend'].values, vmin=-0.1, vmax=0.1, cmap='BrBG', plt_colorbar=True,title='W TrendmCCMP V2.0 CC Adj2')

    diff_cca_20 = trend_maps_v20cca['w_trend'].values - trend_maps_v20['w_trend'].values
    fig,ax = global_map(10.0*diff_cca_20, vmin=-0.02, vmax=0.02, cmap='BrBG', plt_colorbar=True,title='W Trend, CCMP V2.0 CCA - CCMP V2.0')

    plt.show()
    print