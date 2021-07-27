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
from land_fraction.land_fraction import read_land_fraction_1440_720

from plt_zonal_mean_array import plt_zonal_mean_array
from plotting import plot_scat

#RSS statistic accumulators

from rss_stats.hist_1d import Hist1D
from rss_stats.map_stats import MapStat
from rss_maps import zonal_mean, global_mean




def compute_monthly_stats(*,month,year,version,subversion,nc_root,plt_root,save_plots,save_nc_files):

    land_fraction = read_land_fraction_1440_720()

    w_ccmp_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = 0.125,max_lon=360.0-0.125)
    u_ccmp_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                                min_lat = -89.875,max_lat = 89.875,
                                                min_lon = 0.125,max_lon=360.0-0.125)
    v_ccmp_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                                min_lat = -89.875,max_lat = 89.875,
                                                min_lon = 0.125,max_lon=360.0-0.125)

    w_ccmp_hist_1d = Hist1D(num_xbins = 250,min_xval = 0.0,max_xval = 50.0, name = f'W CCMP',units='W (m/s)')
    u_ccmp_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'U CCMP',units='W (m/s)')
    v_ccmp_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'V CCMP',units='W (m/s)')
    
    num_days_in_month = calendar.monthrange(year,month)[1]

    nc_path = f"{nc_root}{version}_{subversion}/"
    plt_path = f"{plt_root}{version}_{subversion}/"
    plot_path_date = f"{plt_path}{year:04}/{month:02}/"
    for day in range(1,num_days_in_month+1):
        version_path = f"P:/CCMP/MEaSUREs/{version}/{subversion}/native/"
        if version == 'nrt':
            version_path = f"P:/CCMP_RT/MEaSUREs/{subversion}/v2.0/native/"
        try:
            ccmp = load_ccmp_daily_arrays(year=year, 
                            month=month,
                            day=day,
                            verbose = True,
                            ccmp_path = version_path,
                            calc_wind_speed = True)
        except FileNotFoundError:
            continue 

        for step in range(0,4):

            u_ccmp_maps.add_map(ccmp['u10'][step,:,:])
            v_ccmp_maps.add_map(ccmp['v10'][step,:,:])
            w_ccmp_maps.add_map(ccmp['w10'][step,:,:])

            ok = land_fraction < 0.01

            u_ccmp_hist_1d.add_data(ccmp['u10'][step,:,:][ok])
            v_ccmp_hist_1d.add_data(ccmp['v10'][step,:,:][ok])
            w_ccmp_hist_1d.add_data(ccmp['w10'][step,:,:][ok])

    u_ccmp_mean_map    = u_ccmp_maps.mean()
    u_ccmp_stddev_map  = u_ccmp_maps.stddev()
    
    v_ccmp_mean_map    = v_ccmp_maps.mean()
    v_ccmp_stddev_map  = v_ccmp_maps.stddev()
 
    w_ccmp_mean_map    = w_ccmp_maps.mean()
    w_ccmp_stddev_map  = w_ccmp_maps.stddev()


    if save_nc_files:
        # write results:
        # accumulated maps:
        ds_ccmp_maps = xr.Dataset(
                data_vars={ 'num'    : (('Latitude', 'Longitude'), u_ccmp_maps.num),
                            'u_tot'    : (('Latitude', 'Longitude'), u_ccmp_maps.tot),
                            'u_tot_sqr'    : (('Latitude', 'Longitude'),u_ccmp_maps.totsqr),
                            'v_tot'    : (('Latitude', 'Longitude'), v_ccmp_maps.tot),
                            'v_tot_sqr'    : (('Latitude', 'Longitude'),v_ccmp_maps.totsqr),                
                            'w_tot'    : (('Latitude', 'Longitude'), w_ccmp_maps.tot),
                            'w_tot_sqr'    : (('Latitude', 'Longitude'),w_ccmp_maps.totsqr),
                            },
                coords={'Latitude': u_ccmp_maps.lats,
                        'Longitude': u_ccmp_maps.lons},
                attrs={ 'CCMP_version'  : f"{version}_{subversion}",
                        'date'          : f"{year:04}_{month:02}"}
            )
        date_str = f"{month:02}_{year:04}"
        ccmp_str = f"CCMP {version} {subversion}"
        os.makedirs(nc_path,exist_ok=True)
        nc_file = f"{nc_path}CCMP_maps_{version}_{subversion}_{date_str}.nc"
        print(nc_file)
        ds_ccmp_maps.to_netcdf(nc_file)

        # 1-D histograms
        ds = xr.Dataset(data_vars = {f'n_u_ccmp' : (('xcenters_vec'),u_ccmp_hist_1d.data['n']),
                                    f'n_v_ccmp' : (('xcenters_vec'),v_ccmp_hist_1d.data['n']),
                                    f'n_w_ccmp'  : (('xcenters_spd'),w_ccmp_hist_1d.data['n']),
                                            },
                                coords = {'xedges_vec'     : u_ccmp_hist_1d.data['xedges'].values,
                                        'xcenters_vec'   : u_ccmp_hist_1d.data['xcenters'].values,
                                        'xedges_spd'     : w_ccmp_hist_1d.data['xedges'].values,
                                        'xcenters_spd'   : w_ccmp_hist_1d.data['xcenters'].values})

        nc_file = f"{nc_path}CCMP_1D_hists_{version}_{subversion}_{date_str}.nc"
        print(nc_file)
        ds.to_netcdf(nc_file)



    fig_umap,ax = global_map(u_ccmp_mean_map, cmap='BrBG', vmin=-15.0, vmax=15.0, plt_colorbar=True,title=f"{ccmp_str}, Mean U (m/s), {date_str}")
    fig_vmap,ax = global_map(v_ccmp_mean_map, cmap='BrBG', vmin=-15.0, vmax=15.0, plt_colorbar=True,title=f"{ccmp_str}, Mean V (m/s), {date_str}")
    fig_wmap,ax = global_map(w_ccmp_mean_map, cmap='viridis', vmin=0.0, vmax=15.0, plt_colorbar=True,title=f"{ccmp_str}, Mean W (m/s), {date_str}")
    
    figu,ax_u = u_ccmp_hist_1d.plot()
    figv,ax_v = v_ccmp_hist_1d.plot()
    figw,ax_w = w_ccmp_hist_1d.plot()
    
    os.makedirs(plot_path_date,exist_ok=True)
    if save_plots:
        png_file = f"{plot_path_date}ccmp_u_mean_map_{date_str}.png"
        fig_umap.savefig(png_file)
        plt.close(fig_umap)

        png_file = f"{plot_path_date}ccmp_v_mean_map_{date_str}.png"
        fig_vmap.savefig(png_file)
        plt.close(fig_vmap)

        png_file = f"{plot_path_date}ccmp_w_mean_map_{date_str}.png"
        fig_wmap.savefig(png_file)
        plt.close(fig_wmap)

        

        png_file = f"{plot_path_date}ccmp_u_1d_hists_{date_str}.png"
        figu.savefig(png_file)
        png_file = f"{plot_path_date}ccmp_v_1d_hists_{date_str}.png"
        figv.savefig(png_file)
        png_file = f"{plot_path_date}ccmp_w_1d_hists_{date_str}.png"
        figw.savefig(png_file)






if __name__ == '__main__':

    for version in ['era5','era5_current_corrected','era5_current_corrected_adj2']:
        subversion = 'v2.0'
        for year in range(2000,2019):
            for month in range(1,13):
                nc_root = f'B:/job_CCMP/monthly_ccmp_stats/nc_files/'
                plt_root = f'B:/job_CCMP/monthly_ccmp_stats/plots/'
                save_plots = True
                save_nc_files = True
                compute_monthly_stats(month=month,
                                    year=year,
                                    version=version,
                                    subversion=subversion,
                                    nc_root=nc_root,
                                    plt_root=plt_root,
                                    save_plots=save_plots,
                                    save_nc_files=save_nc_files)