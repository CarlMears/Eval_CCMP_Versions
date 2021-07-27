from calendar import monthrange
import numpy as np 
import xarray as xr 
import os
import matplotlib.pyplot as plt

import sys
sys.path.append('C:/python_packages_cam/')
sys.path.append('C:/job_CCMP/python/eval_model_wind_vs_buoys/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/job_CCMP/python/ncep/')
sys.path.append('C:/job_CCMP/python/era5/')
sys.path.append('C:/job_CCMP/python/buoys/')
#sys.path.append('../')

from rss_maps import zonal_mean, global_mean
from plt_zonal_mean_array import plt_zonal_mean_array
from plotting import plot_scat
from era5 import convert_to_rss_grid
from era5 import convert_to_rss_grid
from global_map import global_map
from regional_map import regional_map

def read_era5_025deg_data_for_vam(filename):

    u10 = np.zeros((721,1440),dtype=np.float)
    v10 = np.zeros((721,1440),dtype=np.float)
    
    with open(filename,'rb') as f:
        temp = np.fromfile(filename, dtype='>f')
        temp = np.reshape(temp,(2,721,1440))

        u10 = temp[0,:,:]
        v10 = temp[1,:,:]

    return u10,v10

if __name__ == '__main__':

    CCMP_oceanwinds_path = 'P:/CCMP/oceanWinds/vam/'
    
    year = 2013
    month = 1
    day = 2
    time_step = 0

    short_year = year%100

    date_str = f'{month:02d}/{day:02d}/{year:04d} {6*time_step:02d}z'
    date_str_file = f'{year:04d}_{month:02d}_{day:02d}_{6*time_step:02d}z'

    file1 = f"{CCMP_oceanwinds_path}Y{year}/M{month:02d}/D{day:02d}/era5_oscar_corr.025deg.{short_year}{month:02d}{day:02d}.{6*time_step:02d}z.dat"
    file2 = f"{CCMP_oceanwinds_path}Y{year}/M{month:02d}/D{day:02d}/era5_oscar_corr_adj2.{short_year}{month:02d}{day:02d}.{6*time_step:02d}z.dat"

    u10_1,v10_1 = read_era5_025deg_data_for_vam(file1)
    u10_2,v10_2 = read_era5_025deg_data_for_vam(file2)

    w10_1 = np.sqrt(np.square(u10_1)+np.square(v10_1))
    w10_2 = np.sqrt(np.square(u10_2)+np.square(v10_2))
    

    fig_u_ADJ,ax = global_map(u10_2, vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='U, ERA5 ADJ2 '+date_str,units='U Wind (m/s)')
    fig_u_baseline,ax = global_map(u10_1, vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='U, CCMP ERA5 '+date_str,units='U Wind (m/s)')
    fig_u_diff,ax = global_map(u10_2 - u10_1,vmin=-2.0,vmax=2.0,cmap='BrBG',plt_colorbar=True,title='U, ERA5 ADJ2 - ERA5 '+date_str,units=' U Wind Difference (m/s)')

    fig_W_ADJ,ax = global_map(w10_2, vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='W, ERA5 ADJ2 '+date_str,units='Wind Speed (m/s)')
    fig_W_baseline,ax = global_map(w10_1, vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='W, CCMP ERA5 '+date_str,units='Wind Speed (m/s)')
    fig_W_diff,ax = global_map(w10_2 - w10_1,vmin=-2.0,vmax=2.0,cmap='BrBG',plt_colorbar=True,title='W, ERA5 ADJ2 - ERA5 '+date_str,units='Wind Speed Difference (m/s)')

    plt.show()

    '''  Regional Map Stuff
    center = [-29.0,16.0]
    #center = [-10.0,0.0]
    half_size = 8.0
    extent = [center[0] - half_size,
              center[0] + half_size,
              center[1] - half_size,
              center[1] + half_size]
    #extent = [-80.0,-65.,25.0,40.0]

    fig_dudy_ADJ,ax = regional_map(ccmp_ADJ_dudy,vmin=-8.0,vmax=8.0,plt_colorbar=True,cmap='BrBG',title='dudy '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_dvdx_ADJ,ax = regional_map(ccmp_ADJ_dvdx,vmin=-8.0,vmax=8.0,plt_colorbar=True,cmap='BrBG',title='dvdx '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_curl_ADJ,ax = regional_map(ccmp_ADJ_curl,vmin=-10.0,vmax=10.0,plt_colorbar=True,cmap='BrBG',title='cygnss curl '+date_str,extent=extent,units='Wind Speed (m/s)')
    
    fig_dudy_baseline,ax = regional_map(ccmp_baseline_dudy,vmin=-8.0,vmax=8.0,plt_colorbar=True,cmap='BrBG',title='baseline dudy '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_dvdx_baseline,ax = regional_map(ccmp_baseline_dvdx,vmin=-8.0,vmax=8.0,plt_colorbar=True,cmap='BrBG',title='baseline dvdx '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_curl_baseline,ax = regional_map(ccmp_baseline_curl,vmin=-10.0,vmax=10.0,plt_colorbar=True,cmap='BrBG',title='baseline curl '+date_str,extent=extent,units='Wind Speed (m/s)')
  
    fig_curl_diff,ax = regional_map(ccmp_ADJ_curl-ccmp_baseline_curl,vmin=-5.0,vmax=5.0,plt_colorbar=True,cmap='BrBG',title='baseline curl '+date_str,extent=extent,units='Wind Speed (m/s)')

    plt.show()

    fig_w_ADJ,ax = regional_map(ccmp_ADJ_w10,vmin=0.0,vmax=20.0,plt_colorbar=True,title='W, CCMP ERA5 ADJ2 '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_w_baseline,ax = regional_map(ccmp_baseline_w10,vmin=0.0,vmax=20.0,plt_colorbar=True,title='W, CCMP ERA5 '+date_str,extent=extent,units='Wind Speed (m/s)')
    fig_w_diff,ax = regional_map(ccmp_ADJ_w10-ccmp_baseline_w10,vmin=-2.0,vmax=2.0,cmap='BrBG',
    plt_colorbar=True,title='W, CCMP ADJ2 - CCMP ERA5 '+date_str,extent=extent,units='Wind Speed Difference (m/s)')


    fig_u_ADJ,ax = regional_map(ccmp_ADJ_u10,vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='U, CCMP ERA5 ADJ2 '+date_str,extent=extent,units='Zonal Wind (m/s)')
    fig_u_baseline,ax = regional_map(ccmp_baseline_u10,vmin=-20.0,vmax=20.0,cmap='BrBG',plt_colorbar=True,title='U, CCMP ERA5 '+date_str,extent=extent,units='Zonal Wind (m/s)')
    fig_u_diff,ax = regional_map(ccmp_ADJ_u10-ccmp_baseline_u10,vmin=-2.0,vmax=2.0,cmap='BrBG',plt_colorbar=True,title='U, CCMP ADJ2 - CCMP ERA5 '+date_str,extent=extent,units='Zonal Wind Difference(m/s)')

    
    #fig_n_ADJ,ax = regional_map(ccmp_ADJ_nobs,vmin=0,vmax=4.0,plt_colorbar=True,title='N, CCMP ERA5 ADJ2 '+date_str,extent=extent,units='Number of Obs')
    temp = np.copy(ccmp_baseline_nobs)
    temp[temp>1.0] = 1.0

    temp2 = np.copy(temp)
    temp2[(ccmp_ADJ_nobs>ccmp_baseline_nobs)] = 2.0

    fig_n_ADJ,ax = regional_map(temp2,vmin=0,vmax=2.0,plt_colorbar=True,title='N, CCMP ERA5 ADJ2 '+date_str,extent=extent,units='Number of Obs')
    fig_n_baseline,ax = regional_map(temp,vmin=0.0,vmax=2.0,plt_colorbar=True,title='N, CCMP ERA5 '+date_str,extent=extent,units='Number of Obs')

    new_obs = np.zeros((720,1440),dtype = np.int32)
    new_obs[np.all([(ccmp_baseline_nobs == 0),(ccmp_ADJ_nobs > 0)],axis=0)] = 1
    fig_new,ax = regional_map(new_obs,vmin=0.0,vmax=1.0,plt_colorbar=True,title='Obs in gaps  '+date_str,extent=extent)

    plt_path = f'C:/job_CCMP/eval_cygnss/plots/y{year:04d}/m{month:02d}/d{day:02d}/'
    os.makedirs(plt_path,exist_ok = True)
    png_file_u_ADJ = f"{plt_path}u_ADJ_{date_str_file}.png"
    fig_u_ADJ.savefig(png_file_u_ADJ)

    png_file_u_baseline = f"{plt_path}u_baseline_{date_str_file}.png"
    fig_u_baseline.savefig(png_file_u_baseline)

    png_file_u_diff = f"{plt_path}u_diff_{date_str_file}.png"
    fig_u_diff.savefig(png_file_u_diff)

    png_file_w_ADJ = f"{plt_path}w_ADJ_{date_str_file}.png"
    fig_w_ADJ.savefig(png_file_w_ADJ)

    png_file_w_baseline = f"{plt_path}w_baseline_{date_str_file}.png"
    fig_w_baseline.savefig(png_file_w_baseline)

    png_file_w_diff = f"{plt_path}w_diff_{date_str_file}.png"
    fig_w_diff.savefig(png_file_w_diff)

    png_file_n_ADJ = f"{plt_path}n_w_ADJ_{date_str_file}.png"
    fig_n_ADJ.savefig(png_file_n_ADJ)

    png_file_n_baseline = f"{plt_path}n_baseline_{date_str_file}.png"
    fig_n_baseline.savefig(png_file_n_baseline)

    png_file_new = f"{plt_path}new_obs_baseline_{date_str_file}.png"
    fig_new.savefig(png_file_new)
    '''




    