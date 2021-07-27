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

def calc_trend_map(map_arr):

    sz = map_arr.shape
    num_lats = sz[0]
    num_lons = sz[1]
    num_months = sz[2]
    trend_map = np.zeros((num_lats,num_lons))*np.nan

    yrs = 2000.0 + np.arange(0.0,num_months)/12.0
    freq = 1.0
    num_harmonics=2 
    for ilat in range(0,num_lats):
        for ilon in range(0,num_lons):
            ts = map_arr[ilat,ilon,:]
            if np.sum(np.isfinite(ts)) == num_months:
                ft = harmonic_fit(yrs,ts,freq,num_harmonics=num_harmonics,include_trend=True)
                #fig, ax = plt.subplots(1,1)
                #ax.plot(yrs,ts,marker = '+')
                #ax.plot(yrs,ft['fit'],linewidth = 0.3)
                #plt.show()
                trend_map[ilat,ilon] = ft['fit_parameters'][1]
    return trend_map

def calc_ts(map_arr,lat_range=[-80.80],ocean_only = True):
    
    sz = map_arr.shape
 
    num_lats = sz[0]
    num_lons = sz[1]
    num_months = sz[2]

    lat_wt_mp = lat_wt_map(num_lats=num_lats,num_lons=num_lons,lat_range=[-80.0,80.0],ocean_only=True)
    
    ts = np.zeros((num_months))*np.nan
    for m in range(0,num_months):
        ts[m] = global_mean(map_arr[:,:,m],lat_wt_mp = lat_wt_mp)[0]

        
    return ts


start_year = 2000
end_year = 2018
num_years = 1+end_year-start_year

version = 'era5'
subversion = 'current_corrected_adj2_v2.0'
nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'
nc_path = f"{nc_root}{version}_{subversion}/"

u_mean_arr = np.zeros((720,1440,num_years*12))
v_mean_arr = np.zeros((720,1440,num_years*12))
w_mean_arr = np.zeros((720,1440,num_years*12))

for year in range(start_year,end_year+1):
    for month in range(1,13):
        date_str = f"{month:02}_{year:04}"
        nc_file = f"{nc_path}CCMP_maps_{version}_{subversion}_{date_str}.nc"
        ds = xr.open_dataset(nc_file)
        print(nc_file)

        mean_u = ds['u_tot']/ds['num']
        mean_v = ds['v_tot']/ds['num']
        mean_w = ds['w_tot']/ds['num']
    

        month_index = 12*(year-start_year)+month-1

        u_mean_arr[:,:,month_index] = mean_u
        v_mean_arr[:,:,month_index] = mean_v
        w_mean_arr[:,:,month_index] = mean_w

trend_u = 12.0*calc_trend_map(u_mean_arr)
print(f'U trend: Min = {np.nanmin(trend_u)}, Max = {np.nanmax(trend_u)}')
trend_v = calc_trend_map(v_mean_arr)
print(f'V trend: Min = {np.nanmin(trend_v)}, Max = {np.nanmax(trend_v)}')
trend_w = calc_trend_map(w_mean_arr)
print(f'W trend: Min = {np.nanmin(trend_w)}, Max = {np.nanmax(trend_w)}')

ds_ccmp_trend_maps = xr.Dataset(
    data_vars={ 'u_trend'    : (('Latitude', 'Longitude'), trend_u),
                'v_trend'    : (('Latitude', 'Longitude'), trend_v),
                'w_trend'    : (('Latitude', 'Longitude'), trend_w),
                },
    coords={'Latitude': ds['Latitude'],
            'Longitude': ds['Longitude']},
    attrs={ 'CCMP_version'  : f"{version}_{subversion}"}
)

nc_file_out = f"{nc_path}CCMP_trend_maps_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"
print(nc_file_out)
ds_ccmp_trend_maps.to_netcdf(nc_file_out)
print

ocean_only = True
if ocean_only:
    region_str = 'ocean only'
else:
    region_str = 'land and ocean'
glob_u_ts = calc_ts(u_mean_arr,lat_range=[-80,80],ocean_only=True)
glob_v_ts = calc_ts(v_mean_arr,lat_range=[-80,80],ocean_only=True)
glob_w_ts = calc_ts(w_mean_arr,lat_range=[-80,80],ocean_only=True)

yrs = 2000.0 + np.arange(0.0,len(glob_u_ts))/12.0

ds_ccmp_time_series = xr.Dataset(
    data_vars={ 'u_time_series'    : (('Years'),glob_u_ts),
                'v_time_series'    : (('Years'),glob_v_ts),
                'w_time_series'    : (('Years'),glob_w_ts)
                },
    coords={'Years': yrs},
    attrs={ 'CCMP_version'  : f"{version}_{subversion}",
            'extent' : region_str}
)
nc_file_out = f"{nc_path}CCMP_time_series_{version}_{subversion}_{start_year:04d}_{end_year:04d}.nc"
print(nc_file_out)
ds_ccmp_time_series.to_netcdf(nc_file_out)

fig,axs = plt.subplots(nrows=3,ncols=1,sharex=True)


axs[0].plot(yrs,glob_u_ts)
axs[1].plot(yrs,glob_v_ts)
axs[2].plot(yrs,glob_w_ts)  
plt.show()  
print
