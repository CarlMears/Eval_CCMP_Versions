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

def read_ccmp_native(year=2016,month=1,day=1,time_step = 0,path = 'P:/CCMP/MEaSUREs/scl_eraint/v2.0/native/',return_as_xr = False,check_exists=False,verbose=False):

    assert year >= 1988,'year out of range, year must be larger than 1988'
    assert month >= 1,'month out of range'
    assert month <= 12,'month out of range'
    num_days_in_month = monthrange(year, month)[1]
    assert day >= 1,'day out of range'
    assert day <= num_days_in_month,'day out of range'
    time_step = int(time_step)
    assert time_step >= 0,'time step out of range'
    assert time_step <= 3,'time step out of range'

    native_file = path + ('y' + str(year).zfill(4) + '/M'+str(month).zfill(2)+'/D'+str(day).zfill(2) + '/vam_v2.0.analysis.t' +
                         str(year).zfill(4)+str(month).zfill(2)+ str(day).zfill(2) +'.z'+ str(6*time_step).zfill(2)+'.dat')
    if verbose:
        print(native_file)
    if check_exists:
        #just check for existance, return  true if file exists, false otherwise
        if os.path.isfile(native_file):
            return True
        else:
            return False
    print(native_file)
    w = np.fromfile(native_file, dtype='float32').reshape(3, 628, 1440).byteswap()
    
        
    u10 = np.zeros((720,1440),dtype = np.float32)*np.nan
    v10 = np.zeros((720,1440),dtype = np.float32)*np.nan
    nobs = np.zeros((720,1440),dtype = np.float32)*np.nan

    u10[46:674,:]  = w[0, :, :]
    v10[46:674,:]  = w[1, :, :]
    nobs[46:674,:] = w[2, :, :]
    lats = -90.0 + 0.125 + np.arange(0,720)*0.25
    lons = 0.125 + np.arange(0,1440)*0.25

    if return_as_xr:
        ccmp_wind = xr.Dataset(
                        data_vars={'u10':    (('latitude','longitude'),u10),
                                   'v10':    (('latitude','longitude'),v10),
                                   'nobs':   (('latitude','longitude'),nobs)},
                        coords={'latitude'  : lats,
                                'longitude' : lons}
                               )
    else:
        ccmp_wind = dict(lats=lats,
                        lons=lons,
                        u10=u10,
                        v10=v10,
                        nobs=nobs)
    return ccmp_wind

if __name__ == '__main__':

    CCMP_baseline_path = 'P:/CCMP/MEaSUREs/era5_current_corrected/v2.0/native/'
    CCMP_adj_path = 'P:/CCMP/MEaSUREs/era5_current_corrected_adj2/v2.0/native/'

    year = 2016
    month = 6
    day = 2
    time_step = 0
    date_str = f'{month:02d}/{day:02d}/{year:04d} {6*time_step:02d}z'
    date_str_file = f'{year:04d}_{month:02d}_{day:02d}_{6*time_step:02d}z'
    
    ccmp_ADJ = read_ccmp_native(year=year,month=month,day=day,time_step = time_step,path = CCMP_adj_path,return_as_xr = False,check_exists=False,verbose=True)
    ccmp_baseline = read_ccmp_native(year=year,month=month,day=day,time_step = time_step,path = CCMP_baseline_path,return_as_xr = False,check_exists=False,verbose=True)

    ccmp_ADJ_u10 = ccmp_ADJ['u10']
    ccmp_ADJ_v10 = ccmp_ADJ['v10']
    ccmp_ADJ_nobs = ccmp_ADJ['nobs']
    ccmp_baseline_u10 = ccmp_baseline['u10']
    ccmp_baseline_v10 = ccmp_baseline['v10']
    ccmp_baseline_nobs = ccmp_baseline['nobs']

    ccmp_ADJ_dudx,ccmp_ADJ_dudy = np.gradient(ccmp_ADJ_u10,0.25,-0.25)
    ccmp_ADJ_dvdx,ccmp_ADJ_dvdy = np.gradient(ccmp_ADJ_v10,0.25,-0.25)
    ccmp_ADJ_curl = ccmp_ADJ_dvdx - ccmp_ADJ_dudy

    ccmp_baseline_dudx,ccmp_baseline_dudy = np.gradient(ccmp_baseline_u10,0.25,-0.25)
    ccmp_baseline_dvdx,ccmp_baseline_dvdy = np.gradient(ccmp_baseline_v10,0.25,-0.25)
    ccmp_baseline_curl = ccmp_baseline_dvdx - ccmp_baseline_dudy
    

    ccmp_ADJ_w10 = np.sqrt(np.square(ccmp_ADJ['u10']) + np.square(ccmp_ADJ['v10']))
    ccmp_baseline_w10 = np.sqrt(np.square(ccmp_baseline['u10']) + np.square(ccmp_baseline['v10']))

    fig_w_ADJ,ax = global_map(ccmp_ADJ_w10,vmin=0.0,vmax=20.0,plt_colorbar=True,title='W, CCMP ERA5 ADJ2 '+date_str,units='Wind Speed (m/s)')
    fig_w_baseline,ax = global_map(ccmp_baseline_w10,vmin=0.0,vmax=20.0,plt_colorbar=True,title='W, CCMP ERA5 '+date_str,units='Wind Speed (m/s)')
    fig_w_diff,ax = global_map(ccmp_ADJ_w10-ccmp_baseline_w10,vmin=-1.0,vmax=1.0,cmap='BrBG',
    plt_colorbar=True,title='W, CCMP ADJ2 - CCMP ERA5 '+date_str,units='Wind Speed Difference (m/s)')

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




    