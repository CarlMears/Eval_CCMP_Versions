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

def compare_ccmp_to_sat(  year = 2015,month=1,
                            version = 'flk',
                            version_sh = 'V2.0',
                            spatial_subset='globe',
                            compare_sat = 'ASCAT-B',
                            plt_root = None,
                            nc_root = None,
                            save_maps = True, save_binned_plots = True, 
                            plot_maps = True, plot_binned_stats = True, 
                            save_binned_stats = True,
                            plot_1d_hists = True, save_1d_hists = True,save_1d_plots = True,
                            plot_hist_2d = True, save_2d_hist = True,save_2d_plots=True,
                            exclude_ccmp_with_sat = False, exclude_ccmp_no_sat = False):

    if (exclude_ccmp_with_sat and exclude_ccmp_no_sat):
        raise ValueError('Both exclude_ccmp_with_sat and exclude_ccmp_no_sat can not be True') 

    plt_path = f"{plt_root}{version}_{version_sh}/"
    nc_path = f"{nc_root}{version}_{version_sh}/"

    os.makedirs(plt_path,exist_ok=True)
    os.makedirs(nc_path,exist_ok=True)

    if compare_sat in ['ASCAT-A','ASCAT-B','QSCAT']:
        has_vectors = True
    else:
        has_vectors = False

    w_ccmp_minus_sat_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                            min_lat = -89.875,max_lat = 89.875,
                                            min_lon = 0.125,max_lon=360.0-0.125)
    w_ccmp_minus_sat_binned_stats = BinnedStat(num_bins = 100,x_rng=[0.0,50.0])
    w_ccmp_at_sat_hist_1d = Hist1D(num_xbins = 250,min_xval = 0.0,max_xval = 50.0, name = f'W CCMP at {compare_sat}',units='W (m/s)')
    w_sat_hist_1d = Hist1D(num_xbins = 250,min_xval = 0.0,max_xval = 50.0, name = f'W {compare_sat}',units='W (m/s)')
    w_ccmp_sat_hist_2d = Hist2D(num_xbins=250,min_xval=0.0,max_xval=50.0,xname   = f'W CCMP at {compare_sat}',yname = f'W {compare_sat}',xunits='W (m/s)',yunits='W (m/s)')

    if has_vectors:
        u_ccmp_minus_sat_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                                min_lat = -89.875,max_lat = 89.875,
                                                min_lon = 0.125,max_lon=360.0-0.125)

        v_ccmp_minus_sat_maps       =    MapStat(num_lats = 720,num_lons = 1440,
                                                min_lat = -89.875,max_lat = 89.875,
                                                min_lon = 0.125,max_lon=360.0-0.125)

        u_ccmp_minus_sat_binned_stats = BinnedStat(num_bins = 60,x_rng=[-30.0,30.0])
        v_ccmp_minus_sat_binned_stats = BinnedStat(num_bins = 60,x_rng=[-30.0,30.0])

        #set up 1d histograms
        u_ccmp_at_sat_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'U CCMP at {compare_sat}',units='U (m/s)')
        u_sat_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'U {compare_sat}',units='U (m/s)')
    
        v_ccmp_at_sat_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'V CCMP at {compare_sat}',units='V (m/s)')
        v_sat_hist_1d = Hist1D(num_xbins = 500,min_xval = -50.0,max_xval = 50.0, name = f'V {compare_sat}',units='V (m/s)')
    
        #set up 2d histograms
        u_ccmp_sat_hist_2d = Hist2D(num_xbins=500,min_xval=-50.0,max_xval=50.0,xname = f'U CCMP at {compare_sat}',yname = f'U {compare_sat}',xunits='U (m/s)',yunits='U (m/s)')
        v_ccmp_sat_hist_2d = Hist2D(num_xbins=500,min_xval=-50.0,max_xval=50.0,xname = f'V CCMP at {compare_sat}',yname = f'V {compare_sat}',xunits='V (m/s)',yunits='V (m/s)')
    
    num_days_in_month = calendar.monthrange(year,month)[1]
    for day in range(1,num_days_in_month+1):
    #for day in range(1,3):
        if exclude_ccmp_with_sat:
            print(f"{month:02}/{day:02}/{year:04}, {compare_sat}, No Sat")
        else:
            print(f"{month:02}/{day:02}/{year:04}, {compare_sat}, All CCMP")
        
        try:
            if compare_sat == 'ASCAT-B':
                sat = read_rss_satellite_daily_xr(year=year, month=month, day=day,satellite=str.lower(compare_sat),base_dir='//ops1p-to-be-ren/q/')
            else:
                sat = read_rss_satellite_daily_xr(year=year, month=month, day=day,satellite=str.lower(compare_sat))
        except:
            print(f"Data for {compare_sat} missing for {year:04}-{month:02}-{day:02}")
            continue

        if compare_sat in ['ASCAT-A','ASCAT-B']:
            #apply rain flag to ascat
            sat_w, sat_winddir, sat_ok = apply_scat_flags(sat,use_scatflag = True,use_radrain=False)
            #decode u and v from w and direction
            sat_u = sat_w*np.sin(np.deg2rad(sat_winddir))
            sat_v = sat_w*np.cos(np.deg2rad(sat_winddir))
        else:
            sat_w = sat['wspd_lf'].data

        # add a little noise to sat_w to prevent problems with binning due to discretized  ASCAT values
        # this version only used for the histograms.  Anything  that uses differences, or U or V dooes 
        # not show the binning problem

        sat_w_plus_noise = sat_w + np.random.normal(loc=0.0,scale=0.05,size=(2,720,1440))

        version_path = f"P:/CCMP/MEaSUREs/{version}/{version_sh}/native/"
        if version == 'nrt':
            version_path = f"P:/CCMP_RT/MEaSUREs/{version_sh}/v2.0/native/"
        try:
            ccmp = load_ccmp_daily_arrays(year=year, 
                            month=month,
                            day=day,
                            verbose = True,
                            ccmp_path = version_path,
                            calc_wind_speed = True)
        except FileNotFoundError:
            continue       
        
        for node in [0,1]:
            sat_time_map = sat['time'].data[node,:,:]

            ccmp_w_at_sat,ccmp_n_at_sat =  time_interpolate_synoptic_maps(ccmp['w10'], ccmp['time'], sat_time_map,nobs_array =ccmp['nobs'], max_time_diff = 1.0,closest=True)
            w_diff_map = ccmp_w_at_sat - sat_w[node,:,:]
            if has_vectors:
                ccmp_u_at_sat,ignored = time_interpolate_synoptic_maps(ccmp['u10'], ccmp['time'], sat_time_map,nobs_array =ccmp['nobs'], max_time_diff = 1.0,closest=True)
                ccmp_v_at_sat,ignored = time_interpolate_synoptic_maps(ccmp['v10'], ccmp['time'], sat_time_map,nobs_array =ccmp['nobs'], max_time_diff = 1.0,closest=True)

                u_diff_map = ccmp_u_at_sat - sat_u[node,:,:]
                v_diff_map = ccmp_v_at_sat - sat_v[node,:,:]

            #accumulate ccmp-sat maps           
            #w_ccmp_minus_sat_maps.add_map(w_diff_map)

            #accumulate ccmp-sat binned stats
            ok = np.all([np.isfinite(sat_w[node,:,:]),np.isfinite(ccmp_w_at_sat)],axis=(0))
            if exclude_ccmp_no_sat:
                ok = np.all([np.isfinite(sat_w[node,:,:]),np.isfinite(ccmp_w_at_sat),ccmp_n_at_sat > 0],axis=(0))
            if exclude_ccmp_with_sat:
                ok = np.all([np.isfinite(sat_w[node,:,:]),np.isfinite(ccmp_w_at_sat),ccmp_n_at_sat == 0],axis=(0))

            if spatial_subset == 'south40':
                ok[200:719,:] = 0

            data_ok_mask = np.ones_like(ccmp_w_at_sat,dtype=np.int32)
            data_ok_mask[ok] = 0

            map_mask = np.ones_like(ccmp_w_at_sat,dtype=np.float)*np.nan
            map_mask[ok] = 1.0

            #accumulate ccmp-sat maps           
            w_ccmp_minus_sat_maps.add_map(w_diff_map*map_mask)
            if has_vectors:
                u_ccmp_minus_sat_maps.add_map(u_diff_map*map_mask)
                v_ccmp_minus_sat_maps.add_map(v_diff_map*map_mask)

            w_ccmp_minus_sat_binned_stats.add_data(0.5*(sat_w[node,:,:] + ccmp_w_at_sat),ccmp_w_at_sat - sat_w[node,:,:],mask=data_ok_mask)
            if has_vectors:
                v_ccmp_minus_sat_binned_stats.add_data(0.5*(sat_v[node,:,:] + ccmp_v_at_sat),ccmp_v_at_sat - sat_v[node,:,:],mask=data_ok_mask)
                u_ccmp_minus_sat_binned_stats.add_data(0.5*(sat_u[node,:,:] + ccmp_u_at_sat),ccmp_u_at_sat - sat_u[node,:,:],mask=data_ok_mask)

            #accumulate ccmp-sat 1-d histograms
            w_ccmp_at_sat_hist_1d.add_data(ccmp_w_at_sat[ok])
            w_sat_hist_1d.add_data(sat_w_plus_noise[node,:,:][ok])

            if has_vectors:
                u_ccmp_at_sat_hist_1d.add_data(ccmp_u_at_sat[ok])
                u_sat_hist_1d.add_data(sat_u[node,:,:][ok])
                v_ccmp_at_sat_hist_1d.add_data(ccmp_v_at_sat[ok])
                v_sat_hist_1d.add_data(sat_v[node,:,:][ok])
            
            #accumulate 2d hists   
            w_ccmp_sat_hist_2d.add_data(ccmp_w_at_sat[ok],sat_w_plus_noise[node,:,:][ok])
            if has_vectors:         
                u_ccmp_sat_hist_2d.add_data(ccmp_u_at_sat[ok],sat_u[node,:,:][ok])
                v_ccmp_sat_hist_2d.add_data(ccmp_v_at_sat[ok],sat_v[node,:,:][ok])

    # calculate mean and stdddev maps for the month    
    w_ccmp_minus_sat_mean    = w_ccmp_minus_sat_maps.mean()
    w_ccmp_minus_sat_stddev  = w_ccmp_minus_sat_maps.stddev()
    
    if has_vectors:
        u_ccmp_minus_sat_mean    = u_ccmp_minus_sat_maps.mean()
        u_ccmp_minus_sat_stddev  = u_ccmp_minus_sat_maps.stddev()
        v_ccmp_minus_sat_mean    = v_ccmp_minus_sat_maps.mean()
        v_ccmp_minus_sat_stddev  = v_ccmp_minus_sat_maps.stddev()

    ccmp_str = f"CCMP {version} {version_sh}"
    date_str = f"{month:02}_{year:04}"
    if exclude_ccmp_no_sat:
        date_str = f"{month:02}_{year:04}_sat_only"
    if exclude_ccmp_with_sat:
        date_str = f"{month:02}_{year:04}_no_sat"

    #date_str = date_str + '_radrain'

    if spatial_subset != 'global':
        date_str = f"{date_str}_{spatial_subset}"

    if plot_1d_hists:
        figw,ax_w = w_ccmp_at_sat_hist_1d.plot()
        ax_w = w_sat_hist_1d.plot(ax=ax_w)
        if save_1d_plots:
            plot_path_date = f"{plt_path}{year:04}/{month:02}/"
            os.makedirs(plot_path_date,exist_ok=True)
            png_file = f"{plot_path_date}ccmp_vs_{compare_sat}_w_1d_hists_{date_str}.png"
            figw.savefig(png_file)
            plt.close(figw)

        if has_vectors:
            figu,ax_u = u_ccmp_at_sat_hist_1d.plot()
            ax_u =      u_sat_hist_1d.plot(ax=ax_u)
            if save_1d_plots:
                
                png_file = f"{plot_path_date}ccmp_vs_{compare_sat}_u_1d_hists_{date_str}.png"
                figu.savefig(png_file)
                plt.close(figu)

            figv,ax_v = v_ccmp_at_sat_hist_1d.plot()
            ax_v = v_sat_hist_1d.plot(ax=ax_v)
            if save_1d_plots:
                png_file = f"{plot_path_date}ccmp_vs_{compare_sat}_v_1d_hists_{date_str}.png"
                figv.savefig(png_file)
                plt.close(figv)

    if save_1d_hists:
        if has_vectors:
            ds = xr.Dataset(data_vars = {f'n_u_{compare_sat}' : (('xcenters_vec'),u_sat_hist_1d.data['n']),
                                        f'n_u_ccmp_{compare_sat}':(('xcenters_vec'),u_ccmp_at_sat_hist_1d.data['n']),
                                        f'n_v_{compare_sat}' : (('xcenters_vec'),v_sat_hist_1d.data['n']),
                                        f'n_v_ccmp_{compare_sat}':(('xcenters_vec'),v_ccmp_at_sat_hist_1d.data['n']),
                                        f'n_w_{compare_sat}' : (('xcenters_spd'),w_sat_hist_1d.data['n']),
                                        f'n_w_ccmp_{compare_sat}':(('xcenters_spd'),w_ccmp_at_sat_hist_1d.data['n']),
                                        },
                            coords = {'xedges_vec'     : u_sat_hist_1d.data['xedges'].values,
                                    'xcenters_vec'   : u_sat_hist_1d.data['xcenters'].values,
                                    'xedges_spd'     : w_sat_hist_1d.data['xedges'].values,
                                    'xcenters_spd'   : w_sat_hist_1d.data['xcenters'].values})
        else:
            ds = xr.Dataset(data_vars = {f'n_w_{compare_sat}' : (('xcenters_spd'),w_sat_hist_1d.data['n']),
                                         f'n_w_ccmp{compare_sat}':(('xcenters_spd'),w_ccmp_at_sat_hist_1d.data['n']),
                                        },
                            coords = {'xedges_spd'     : w_sat_hist_1d.data['xedges'].values,
                                    'xcenters_spd'   : w_sat_hist_1d.data['xcenters'].values})

            
        nc_file = f"{nc_path}/ccmp_{compare_sat}_1d_hists_{date_str}.nc"
        ds.to_netcdf(nc_file)

    if plot_hist_2d:
        plot_path_date = f"{plt_path}{year:04}/{month:02}/"
        os.makedirs(plot_path_date,exist_ok=True)

        ccmp_sat_w_fig,ax = w_ccmp_sat_hist_2d.plot(title=f'CCMP - {compare_sat}, Wind Speed', ytitle='W CCMP (m/s)', xtitle=f'W {compare_sat} (m/s)', 
             aspect='equal', plot_diagonal=True)
        
        if has_vectors:
            ccmp_sat_v_fig,ax = v_ccmp_sat_hist_2d.plot(title=f'CCMP - {compare_sat}, Meridional Wind', ytitle='V CCMP (m/s)', xtitle=f'V {compare_sat} (m/s)', 
                aspect='equal', plot_diagonal=True)
            
            ccmp_sat_u_fig,ax = u_ccmp_sat_hist_2d.plot(title=f'CCMP - {compare_sat}, Zonal Wind', ytitle='U CCMP (m/s)', xtitle=f'U {compare_sat} (m/s)', 
                aspect='equal', plot_diagonal=True)
       

        if save_2d_plots:
            png_file = f"{plot_path_date}ccmp_{compare_sat}_w_hist_2d_{date_str}.png"
            ccmp_sat_w_fig.savefig(png_file)

            if has_vectors:
                png_file = f"{plot_path_date}ccmp_{compare_sat}_u_hist_2d_{date_str}.png"
                ccmp_sat_u_fig.savefig(png_file)

                png_file = f"{plot_path_date}ccmp_{compare_sat}_v_hist_2d_{date_str}.png"
                ccmp_sat_v_fig.savefig(png_file) 

            
        if has_vectors:
            plt.close(ccmp_sat_u_fig)
            plt.close(ccmp_sat_v_fig)
        plt.close(ccmp_sat_w_fig)

    if save_2d_hist:
        if has_vectors:
            ds = xr.Dataset(data_vars = {f'n_u_ccmp_{compare_sat}':(('xcenters_vec','ycenters_vec'),u_ccmp_sat_hist_2d.data['n']),
                                         f'n_v_ccmp_{compare_sat}':(('xcenters_vec','ycenters_vec'),v_ccmp_sat_hist_2d.data['n']),
                                         f'n_w_ccmp_{compare_sat}':(('xcenters_spd','ycenters_spd'),w_ccmp_sat_hist_2d.data['n']),
                                         },
                            coords = {'xedges_vec'     : u_ccmp_sat_hist_2d.data['xedges'].values,
                                    'xcenters_vec'   : u_ccmp_sat_hist_2d.data['xcenters'].values,
                                    'xedges_spd'     : w_ccmp_sat_hist_2d.data['xedges'].values,
                                    'xcenters_spd'   : w_ccmp_sat_hist_2d.data['xcenters'].values,
                                    'yedges_vec'     : u_ccmp_sat_hist_2d.data['yedges'].values,
                                    'ycenters_vec'   : u_ccmp_sat_hist_2d.data['ycenters'].values,
                                    'yedges_spd'     : w_ccmp_sat_hist_2d.data['yedges'].values,
                                    'ycenters_spd'   : w_ccmp_sat_hist_2d.data['ycenters'].values})
        else:
            ds = xr.Dataset(data_vars = {f'n_w_ccmp_{compare_sat}':(('xcenters_spd','ycenters_spd'),w_ccmp_sat_hist_2d.data['n']),
                                         },
                            coords = {'xedges_spd'     : w_ccmp_sat_hist_2d.data['xedges'].values,
                                    'xcenters_spd'   : w_ccmp_sat_hist_2d.data['xcenters'].values,
                                    'yedges_spd'     : w_ccmp_sat_hist_2d.data['yedges'].values,
                                    'ycenters_spd'   : w_ccmp_sat_hist_2d.data['ycenters'].values})

        nc_file = f"{nc_path}/ccmp_{compare_sat}_2d_hists_{date_str}.nc"
        ds.to_netcdf(nc_file)

    if save_binned_stats:
        nc_file = f"{nc_path}/w_ccmp_minus_{compare_sat}_binned_stats_{date_str}.nc"
        w_ccmp_minus_sat_binned_stats.as_DataArray().to_netcdf(nc_file)
        if has_vectors:
            nc_file = f"{nc_path}/u_ccmp_minus_{compare_sat}_binned_stats_{date_str}.nc"
            u_ccmp_minus_sat_binned_stats.as_DataArray().to_netcdf(nc_file)
            nc_file = f"{nc_path}/v_ccmp_minus_{compare_sat}_binned_stats_{date_str}.nc"
            v_ccmp_minus_sat_binned_stats.as_DataArray().to_netcdf(nc_file)
        

    if plot_binned_stats:
        plot_path_date = f"{plt_path}{year:04}/{month:02}/"
        os.makedirs(plot_path_date,exist_ok = True)

        fig9,ax = w_ccmp_minus_sat_binned_stats.plot(yrng=[-4.0,4.0], xrng=[0.0,40.0],xlab='Mean Wind', ylab=f'CCMP - {compare_sat}',title=f'{compare_sat} Wind Speed, {ccmp_str}, {date_str}')
        if save_binned_plots:
            png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_w_binned_diffs_{date_str}.png"
            fig9.savefig(png_file)
        plt.close(fig9)

        if has_vectors:
            fig5,ax = u_ccmp_minus_sat_binned_stats.plot(yrng=[-4.0,4.0], xrng=[-30.0,30.0],xlab='Mean Wind', ylab=f'CCMP - {compare_sat}',title=f'{compare_sat} U Wind, {ccmp_str}, {date_str}')
            if save_binned_plots:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_u_binned_diffs_{date_str}.png"
                fig5.savefig(png_file)
            plt.close(fig5)

            fig7,ax = v_ccmp_minus_sat_binned_stats.plot(yrng=[-4.0,4.0], xrng=[-30.0,30.0],xlab='Mean Wind', ylab=f'CCMP -{compare_sat}',title=f'{compare_sat} V Wind, {ccmp_str}, {date_str}')
            if save_binned_plots:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_v_binned_diffs_{date_str}.png"
                fig7.savefig(png_file)
            plt.close(fig7)



    if plot_maps:

        fig34,ax = global_map_w_zonal_mean(w_ccmp_minus_sat_mean, cmap='BrBG', vmin=-1.5, vmax=1.5, zmin=-1.0, zmax=1.0,plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, W Bias (m/s), {date_str}")
        if save_maps:
            png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_w_diff_map_{date_str}.png"
            fig34.savefig(png_file)
        plt.close(fig34)

        fig36,ax = global_map_w_zonal_mean(w_ccmp_minus_sat_stddev, vmin=0.0, vmax=3.0, zmin=0.0, zmax=2.5, plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, W StdDev (m/s), {date_str}")
        if save_maps:
            png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_w_sdev_map_{date_str}.png"
            fig36.savefig(png_file)
        plt.close(fig36)

        if has_vectors:
            fig14,ax = global_map_w_zonal_mean(u_ccmp_minus_sat_mean, cmap='BrBG', vmin=-1.5, vmax=1.5,  zmin=-1.0, zmax=1.0, plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, U Bias (m/s), {date_str}")
            if save_maps:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_u_diff_map_{date_str}.png"
                fig14.savefig(png_file)
            plt.close(fig14)

            fig16,ax = global_map_w_zonal_mean(u_ccmp_minus_sat_stddev, vmin=0.0, vmax=3.0, zmin=0.0, zmax=2.5, plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, U StdDev (m/s), {date_str}")
            if save_maps:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_u_sdev_map_{date_str}.png"
                fig16.savefig(png_file)
            plt.close(fig16)

            fig24,ax = global_map_w_zonal_mean(v_ccmp_minus_sat_mean, cmap='BrBG', vmin=-1.5, vmax=1.5,  zmin=-1.0, zmax=1.0, plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, V Bias (m/s), {date_str}")
            if save_maps:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_v_diff_map_{date_str}.png"
                fig24.savefig(png_file)
            plt.close(fig24)

            fig26,ax = global_map_w_zonal_mean(v_ccmp_minus_sat_stddev, vmin=0.0, vmax=3.0, zmin=0.0, zmax=2.5, plt_colorbar=True,title=f"{ccmp_str} - {compare_sat}, V StdDev (m/s), {date_str}")
            if save_maps:
                png_file = f"{plot_path_date}ccmp_minus_{compare_sat}_v_sdev_map_{date_str}.png"
                fig26.savefig(png_file)
            plt.close(fig26)




    #assemble and xarray dataset with all the needed maps

    if has_vectors:
        ds_ccmp_minus_sat = xr.Dataset(
            data_vars={'num'    : (('Latitude', 'Longitude'), u_ccmp_minus_sat_maps.num),
                        'u_tot'    : (('Latitude', 'Longitude'), u_ccmp_minus_sat_maps.tot),
                        'u_tot_sqr'    : (('Latitude', 'Longitude'),u_ccmp_minus_sat_maps.totsqr),
                        'v_tot'    : (('Latitude', 'Longitude'), v_ccmp_minus_sat_maps.tot),
                        'v_tot_sqr'    : (('Latitude', 'Longitude'),v_ccmp_minus_sat_maps.totsqr),                
                        'w_tot'    : (('Latitude', 'Longitude'), w_ccmp_minus_sat_maps.tot),
                        'w_tot_sqr'    : (('Latitude', 'Longitude'),w_ccmp_minus_sat_maps.totsqr),
                        },
            coords={'Latitude': u_ccmp_minus_sat_maps.lats,
                    'Longitude': u_ccmp_minus_sat_maps.lons},
            attrs={ 'satellite'      : compare_sat,
                    'CCMP_version'  : f"{version}_{version_sh}",
                    'date'          : f"{year:04}_{month:02}"}
        )
    else:
         ds_ccmp_minus_sat = xr.Dataset(
            data_vars={'num'    : (('Latitude', 'Longitude'), w_ccmp_minus_sat_maps.num),
                       'w_tot'    : (('Latitude', 'Longitude'), w_ccmp_minus_sat_maps.tot),
                       'w_tot_sqr'    : (('Latitude', 'Longitude'),w_ccmp_minus_sat_maps.totsqr),
                        },
            coords={'Latitude': w_ccmp_minus_sat_maps.lats,
                    'Longitude': w_ccmp_minus_sat_maps.lons},
            attrs={ 'satellite'      : compare_sat,
                    'CCMP_version'  : f"{version}_{version_sh}",
                    'date'          : f"{year:04}_{month:02}"}
        )

    
    nc_file = f"{nc_path}/CCMP_vs_{compare_sat}_{version}_{date_str}.nc"
    print(nc_file)
    ds_ccmp_minus_sat.to_netcdf(nc_file)

year_list = range(2013,2020)
#month_list = [1,2,4,5,7,8,10,11]
month_list = [3,6,9,12]


for year in year_list:
    for month in month_list: 
        compare_sat = 'ASCAT-B'   

        compare_ccmp_to_sat(  year = year,month=month,
                                version = 'era5_current_corrected_adj2',
                                version_sh = 'v2.0',
                                spatial_subset='global',
                                compare_sat = compare_sat,
                                exclude_ccmp_with_sat = True,
                                plt_root =f'B:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/',
                                nc_root = f'B:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/')
        compare_ccmp_to_sat(  year = year,month=month,
                                version = 'era5_current_corrected_adj2',
                                version_sh = 'v2.0',
                                spatial_subset='global',
                                compare_sat = compare_sat,
                                exclude_ccmp_with_sat = False,
                                plt_root =f'B:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/',
                                nc_root = f'B:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/')

        
