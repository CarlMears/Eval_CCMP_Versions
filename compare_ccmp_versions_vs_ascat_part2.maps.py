import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import calendar
import copy
import os

import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/job_CCMP/python/binned_stats/')
sys.path.append('C:/job_CCMP/python/map_stats/')
sys.path.append('C:/job_CCMP/python/CCMP/')
sys.path.append('C:/job_CCMP/python/eval_era5/')
sys.path.append('C:/job_CCMP/python/era5/')
sys.path.append('C:/job_CCMP/python/buoys/')
sys.path.append('B:/job_CCMP/python/bytemaps/')
sys.path.append('B:/job_CCMP/python/quikscat/')
sys.path.append('B:/job_CCMP/python/ascat/')
sys.path.append('B:/job_CCMP/python/ssmi/')
sys.path.append('B:/job_CCMP/python/ssmis/')
sys.path.append('B:/job_CCMP/python/amsre/')
sys.path.append('B:/job_CCMP/python/amsr2/')
sys.path.append('B:/job_CCMP/python/windsat/')
sys.path.append('C:/python_packages_cam/')
sys.path.append('../')

#from ccmp import read_ccmp_native,load_ccmp_daily_arrays
#from era5 import convert_to_rss_grid
#from time_interpolate_synoptic_maps import time_interpolate_synoptic_maps

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

import copy
import os


compare_sats = ['ASCAT-B']
ccmp_versions = ['era5_current_corrected']
#ccmp_versions = ['era5']
ccmp_subversions = ['v2.0']
exclude_ccmp_with_sat = True
exclude_ccmp_no_sat = False

start_year=2013
end_year = 2019
num_months = 12*(1 + end_year - start_year)
spatial_mean_ts = np.zeros((2,2,3,num_months))
spatial_stddev_ts = np.zeros((2,2,3,num_months))

zonal_means = np.zeros((2,2,3,720))  #ccmp_version,compare_sat,u_v_w,latitude
zonal_mean_stddevs = np.zeros((2,2,3,720))  #ccmp_version,compare_sat,u_v_w,latitude

year_list = np.arange(start_year,end_year+1)

for j,ccmp_version in enumerate(ccmp_versions):
    ccmp_subversion = ccmp_subversions[j] 

    year_list = np.arange(start_year,end_year+1)
    
    for k,compare_sat in enumerate(compare_sats):

        nc_path = f"C:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/{ccmp_version}_{ccmp_subversion}/" 
        zonal_mean_u_time_lat = np.zeros((num_months,720))
        zonal_mean_v_time_lat = np.zeros((num_months,720))
        zonal_mean_w_time_lat = np.zeros((num_months,720))

        zonal_stddev_u_time_lat = np.zeros((num_months,720))
        zonal_stddev_v_time_lat = np.zeros((num_months,720))
        zonal_stddev_w_time_lat = np.zeros((num_months,720))

        first_month = True
        for  year in year_list:
            for month in np.arange(1,13,dtype=np.int):
                date_str = f"{month:02}_{year:04}"
                if exclude_ccmp_no_sat:
                    date_str = f"{month:02}_{year:04}_sat_only"
                if exclude_ccmp_with_sat:
                    date_str = f"{month:02}_{year:04}_no_sat"

                month_index = (month -1) + 12*(year - start_year)
                nc_file = f"{nc_path}CCMP_vs_{compare_sat}_{ccmp_version}_{date_str}.nc"
                print(nc_file)
                ccmp_vs_ascat_u_maps,ccmp_vs_ascat_v_maps,ccmp_vs_ascat_w_maps = MapStat.from_netcdf_triple(nc_file = nc_file)
                
                if first_month:
                    ccmp_vs_ascat_u_maps_tot = copy.copy(ccmp_vs_ascat_u_maps)
                    ccmp_vs_ascat_v_maps_tot = copy.copy(ccmp_vs_ascat_v_maps)
                    ccmp_vs_ascat_w_maps_tot = copy.copy(ccmp_vs_ascat_w_maps)
                    first_month = False
                else:
                    ccmp_vs_ascat_u_maps_tot.combine(ccmp_vs_ascat_u_maps)
                    ccmp_vs_ascat_v_maps_tot.combine(ccmp_vs_ascat_v_maps)
                    ccmp_vs_ascat_w_maps_tot.combine(ccmp_vs_ascat_w_maps)

                #compute spatial mean and spatial stddev
                for n,map in enumerate([ccmp_vs_ascat_u_maps,ccmp_vs_ascat_v_maps,ccmp_vs_ascat_w_maps]):
                    mn_diff_map = map.mean()
                    mn_diff_map[map.num_obs() <= 5] = np.nan
                    spatial_mean = np.nanmean(mn_diff_map)
                    spatial_stddev = np.nanstd(mn_diff_map)
                    print(j,k,n,spatial_mean,spatial_stddev)
                    spatial_mean_ts[j,k,n,month_index] = spatial_mean
                    spatial_stddev_ts[j,k,n,month_index] = spatial_stddev
                

                zm_u  = ccmp_vs_ascat_u_maps.zonal_mean(min_good_pix = 200)
                zm_v  = ccmp_vs_ascat_v_maps.zonal_mean(min_good_pix = 200)
                zm_w  = ccmp_vs_ascat_w_maps.zonal_mean(min_good_pix = 200)


                zs_u  = ccmp_vs_ascat_u_maps.zonal_stddev(min_good_pix = 200)
                zs_v  = ccmp_vs_ascat_v_maps.zonal_stddev(min_good_pix = 200)
                zs_w  = ccmp_vs_ascat_w_maps.zonal_stddev(min_good_pix = 200)

                
                zonal_mean_u_time_lat[month_index,:] = zm_u
                zonal_mean_v_time_lat[month_index,:] = zm_v
                zonal_mean_w_time_lat[month_index,:] = zm_w

                zonal_stddev_u_time_lat[month_index,:] = zs_u
                zonal_stddev_v_time_lat[month_index,:] = zs_v
                zonal_stddev_w_time_lat[month_index,:] = zs_w

                '''fig,ax = plt.subplots()
                ax.plot(zs_u)
                ax.plot(zs_v,color='red')
                ax.plot(zs_w,color='green')
                
                plt.show()
                print()'''

        num_obs          = ccmp_vs_ascat_u_maps_tot.num_obs()
        mx_num = np.nanmax(num_obs)
        mean_ccmp_u_diff = ccmp_vs_ascat_u_maps_tot.mean()
        mean_ccmp_v_diff = ccmp_vs_ascat_v_maps_tot.mean()
        mean_ccmp_w_diff = ccmp_vs_ascat_w_maps_tot.mean()

        std_ccmp_u_diff = ccmp_vs_ascat_u_maps_tot.stddev()
        std_ccmp_v_diff = ccmp_vs_ascat_v_maps_tot.stddev()
        std_ccmp_w_diff = ccmp_vs_ascat_w_maps_tot.stddev()

        zonal_mean_ccmp_u_diff = ccmp_vs_ascat_u_maps_tot.zonal_mean(min_good_pix = 144)
        zonal_mean_ccmp_v_diff = ccmp_vs_ascat_v_maps_tot.zonal_mean(min_good_pix = 144)
        zonal_mean_ccmp_w_diff = ccmp_vs_ascat_w_maps_tot.zonal_mean(min_good_pix = 144)

        zonal_means[j,k,0,:] =  zonal_mean_ccmp_u_diff  #ccmp_version,compare_sat,u_v_w,latitude
        zonal_means[j,k,1,:] =  zonal_mean_ccmp_v_diff
        zonal_means[j,k,2,:] =  zonal_mean_ccmp_w_diff

        zonal_mean_ccmp_u_stddev = ccmp_vs_ascat_u_maps_tot.zonal_stddev(min_good_pix = 144)
        zonal_mean_ccmp_v_stddev = ccmp_vs_ascat_v_maps_tot.zonal_stddev(min_good_pix = 144)
        zonal_mean_ccmp_w_stddev = ccmp_vs_ascat_w_maps_tot.zonal_stddev(min_good_pix = 144)

        zonal_mean_stddevs[j,k,0,:] =  zonal_mean_ccmp_u_stddev  #ccmp_version,compare_sat,u_v_w,latitude
        zonal_mean_stddevs[j,k,1,:] =  zonal_mean_ccmp_v_stddev
        zonal_mean_stddevs[j,k,2,:] =  zonal_mean_ccmp_w_stddev



        mean_ccmp_u_diff[num_obs < 100.0] = np.nan
        mean_ccmp_v_diff[num_obs < 100.0] = np.nan
        mean_ccmp_w_diff[num_obs < 100.0] = np.nan

        map_bias_vmin = -1.0
        map_bias_vmax = 1.0

        title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, Number of Collocations"
        if exclude_ccmp_with_sat:
            title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, Number of Collocations, No Sat"
        num_fig,ax       = global_map_w_zonal_mean(num_obs,         cmap='magma', vmin= 0.0, vmax=1000.0, plt_colorbar=True,title=title_str)
        
        title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, U Bias (m/s)"
        if exclude_ccmp_with_sat:
            title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, U Bias (m/s), No Sat"
        umap_diff_fig,ax = global_map_w_zonal_mean(mean_ccmp_u_diff, cmap='BrBG', vmin=map_bias_vmin, vmax=map_bias_vmax,zmin=map_bias_vmin, zmax=map_bias_vmax, plt_colorbar=True,title=title_str)
        
        title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, V Bias (m/s)"
        if exclude_ccmp_with_sat:
            title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, V Bias (m/s), No Sat"
        vmap_diff_fig,ax = global_map_w_zonal_mean(mean_ccmp_v_diff, cmap='BrBG', vmin=map_bias_vmin, vmax=map_bias_vmax, zmin=map_bias_vmin, zmax=map_bias_vmax,plt_colorbar=True,title=title_str)
        
        title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, W Bias (m/s)"
        if exclude_ccmp_with_sat:
            title_str = f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, W Bias (m/s), No Sat"
        wmap_diff_fig,ax = global_map_w_zonal_mean(mean_ccmp_w_diff, cmap='BrBG', vmin=map_bias_vmin, vmax=map_bias_vmax, zmin=map_bias_vmin, zmax=map_bias_vmax,plt_colorbar=True,title=title_str)
       
        #umap_diff_fig3,ax = global_map(mean_ccmp_u_diff, cmap='BrBG', vmin=-1.0, vmax=1.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year},U Bias (m/s)",extent=[90,270,-30,30],central_longitude=180.0)
        #vmap_diff_fig3,ax = global_map(mean_ccmp_v_diff, cmap='BrBG', vmin=-1.0, vmax=1.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, V Bias (m/s)",extent=[90,270,-30,30],central_longitude=180.0)
        #wmap_diff_fig3,ax = global_map(mean_ccmp_w_diff, cmap='BrBG', vmin=-1.0, vmax=1.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, W Bias (m/s)",extent=[90,270,-30,30],central_longitude=180.0)

        umap_diff_fig4,ax = global_map_w_zonal_mean(std_ccmp_u_diff, cmap='magma', vmin=0.0, vmax=2.0, zmin=0.0, zmax=2.0,plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year},U Std. Dev. (m/s)")
        vmap_diff_fig4,ax = global_map_w_zonal_mean(std_ccmp_v_diff, cmap='magma', vmin=0.0, vmax=2.0, zmin=0.0, zmax=2.0,plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, V Std. Dev. (m/s)")
        wmap_diff_fig4,ax = global_map_w_zonal_mean(std_ccmp_w_diff, cmap='magma', vmin=0.0, vmax=2.0, zmin=0.0, zmax=2.0,plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, W Std. Dev. (m/s)")
      
        #umap_diff_fig4,ax = global_map(std_ccmp_u_diff, cmap='magma', vmin=0.0, vmax=2.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year},U Std. Dev. (m/s)",extent=[90,270,-30,30],central_longitude=180.0)
        #vmap_diff_fig4,ax = global_map(std_ccmp_v_diff, cmap='magma', vmin=0.0, vmax=2.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, V Std. Dev. (m/s)",extent=[90,270,-30,30],central_longitude=180.0)
        #wmap_diff_fig4,ax = global_map(std_ccmp_w_diff, cmap='magma', vmin=0.0, vmax=2.0, plt_colorbar=True,title=f"{ccmp_version} - {compare_sat}, {start_year}-{end_year}, W Std. Dev. (m/s)",extent=[90,270,-30,30],central_longitude=180.0)
      
        plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/map_summaries/'

        os.makedirs(plt_path,exist_ok=True)
        png_ending = '.png'
        title_ending = ', All'
        if exclude_ccmp_with_sat:
            png_ending = '.nosat.png'
            title_ending = ', No Sat'
        
        png_file = f"{plt_path}nobs_CCMP_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        num_fig.savefig(png_file)
        plt.close(num_fig)
        png_file = f"{plt_path}mean_CCMP_U_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        umap_diff_fig.savefig(png_file)
        plt.close(umap_diff_fig)
        png_file = f"{plt_path}mean_CCMP_V_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        vmap_diff_fig.savefig(png_file)
        plt.close(vmap_diff_fig)
        png_file = f"{plt_path}mean_CCMP_W_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        wmap_diff_fig.savefig(png_file)
        plt.close(wmap_diff_fig)

        '''png_file = f"{plt_path}mean_CCMP_U_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}.trop_pac{png_ending}"
        umap_diff_fig3.savefig(png_file)
        plt.close(umap_diff_fig3)
        png_file = f"{plt_path}mean_CCMP_V_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}.trop_pac{png_ending}"
        vmap_diff_fig3.savefig(png_file)
        plt.close(vmap_diff_fig3)
        png_file = f"{plt_path}mean_CCMP_W_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}.trop_pac{png_ending}"
        wmap_diff_fig3.savefig(png_file)
        plt.close(wmap_diff_fig3)
        '''

        png_file = f"{plt_path}std_CCMP_U_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        umap_diff_fig4.savefig(png_file)
        plt.close(umap_diff_fig4)
        png_file = f"{plt_path}std_CCMP_V_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        vmap_diff_fig4.savefig(png_file)
        plt.close(vmap_diff_fig4)
        png_file = f"{plt_path}std_CCMP_W_{ccmp_version}_minus_{compare_sat}_{start_year}-{end_year}{png_ending}"
        wmap_diff_fig4.savefig(png_file)
        plt.close(wmap_diff_fig4)


        plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/time_lat_summaries/'
        os.makedirs(plt_path,exist_ok = True)
       
        yrs = start_year + np.arange(0,num_months)/12.0 + 1.0/24.0
        X,Y = np.meshgrid(yrs,ccmp_vs_ascat_w_maps_tot.lats)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_mean_u_time_lat), cmap='BrBG', vmin  = -0.6,vmax=0.6)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Mean U Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}{title_ending}")
        png_file = f"{plt_path}mean_CCMP_U_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_mean_v_time_lat), cmap='BrBG', vmin  = -0.6,vmax=0.6)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Mean V Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}{title_ending}")
        png_file = f"{plt_path}mean_CCMP_V_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)
        
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_mean_w_time_lat), cmap='BrBG', vmin  = -0.6,vmax=0.6)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Mean W Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}{title_ending}")
        png_file = f"{plt_path}mean_CCMP_W_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)
        
        

        

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_stddev_u_time_lat), cmap='magma', vmin  = 0.0,vmax=2.0)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Stddev U Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}")
        png_file = f"{plt_path}stddev_CCMP_U_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)
                

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_stddev_v_time_lat), cmap='magma', vmin  = 0.0,vmax=2.0)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Stddev V Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}")
        png_file = f"{plt_path}stddev_CCMP_V_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X,Y, np.transpose(zonal_stddev_w_time_lat), cmap='magma', vmin  = 0.0,vmax=2.0)
        cbar = fig.colorbar(im)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel('Stddev W Difference', fontsize=16)
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Year')
        ax.set_title(f"{ccmp_subversion} minus {compare_sat}")
        png_file = f"{plt_path}stddev_CCMP_W_{ccmp_version}_minus_{compare_sat}_time_lat{png_ending}"
        fig.savefig(png_file)
        plt.close(fig)

plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/zonal_mean_summaries/'
os.makedirs(plt_path,exist_ok = True)

lats = -90.0 + 0.125 + 0.25*np.arange(0,720)

for ivar,var_name in enumerate(['U','V','W']):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    for j,ccmp_version in enumerate(ccmp_versions):
        ccmp_subversion = ccmp_subversions[j] 

        for k,compare_sat in enumerate(['ASCAT-A','ASCAT-B']):
            ax.plot(zonal_means[j,k,ivar,:],lats,label = f'{ccmp_subversion}-{compare_sat}')
    ax.plot([0.0,0.0],[-90.0,90.0],color='grey')
    ax.legend()
    ax.set_xlim((-0.4,0.4))
    ax.set_ylim((-90.0,90.0))
    ax.legend(loc='upper left')
    ax.set_xlabel = f"{var_name}, CCMP - ASCAT"
    ax.set_ylabel = 'Latitude'
    png_file = f"{plt_path}CCMP_minus_ASCAT_zonal_mean_{var_name}{png_ending}"
    fig.savefig(png_file)
print()

plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/time_series_summaries/'
os.makedirs(plt_path,exist_ok = True)

for ivar,var_name in enumerate(['U','V','W']):
    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    for j,ccmp_version in enumerate(ccmp_versions):
        ccmp_subversion = ccmp_subversions[j] 

        for k,compare_sat in enumerate(compare_sats):
            ax.plot(yrs,spatial_mean_ts[j,k,ivar],label = f'{ccmp_version}-{compare_sat}')
    ax.plot(yrs,np.zeros(len(yrs)),color='grey')
    ax.set_title(f"Spatial Mean, {var_name}")
    ax.set_ylim((-1.0,1.0))
    ax.legend()
    png_file = f"{plt_path}{var_name}_CCMP_minus_ASCAT_spatial_mean_ts{png_ending}"
    fig.savefig(png_file)


for ivar,var_name in enumerate(['U','V','W']):
    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    for j,ccmp_version in enumerate(ccmp_versions):
        ccmp_subversion = ccmp_subversions[j] 

        for k,compare_sat in enumerate(compare_sats):
            ax.plot(yrs,spatial_stddev_ts[j,k,ivar],label = f'{ccmp_version}-{compare_sat}')

    ax.set_title(f"Spatial Std. Dev., {var_name}")
    ax.set_ylim((0.0,1.0))
    ax.legend()
    png_file = f"{plt_path}{var_name}_CCMP_minus_ASCAT_spatial_stddev_ts{png_ending}"
    fig.savefig(png_file)

plt.show()



print()