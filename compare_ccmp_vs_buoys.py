'''
Compares CCMP  to buoy observations.

Any CCMP version can be used by adjusting the nrt, ccmp_version and  ccmp_subversion parameters
Either the RSS or NAVY buoy dataset can be choosen by adjusting the compare_buoys = 'RSS'

Right now, not totally functional

'''


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import calendar
import copy

import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/python_packages_cam/rss_stats/')
sys.path.append('C:/job_CCMP/python/CCMP/')
sys.path.append('C:/python_packages_cam/wind_compare_tools/')
sys.path.append('c:/job_CCMP/python/buoys/')
sys.path.append('../')

from rss_maps import zonal_mean, global_mean
from plt_zonal_mean_array import plt_zonal_mean_array
from plotting import plot_scat
from ccmp import read_ccmp_native,load_ccmp_daily_arrays

from global_map import global_map
from binned_stats import BinnedStat
from hist_1d import Hist1D
from hist_2d import Hist2D
from map_stats import MapStat
from rss_maps import zonal_mean, global_mean

from compare_daily_map_to_buoys import colocate_daily_wind_maps_to_buoys
from Eq_Neutral_Wind import EN_wind_table

def convert_to_rss_grid(w):  
    import numpy as np  
    temp = np.zeros((721,1441))
    temp[:,0:1440] = w
    temp[:,1440] = w[:,0]
    return np.flipud(0.25*(temp[0:720,0:1440] + 
                           temp[1:721,0:1440] +
                           temp[0:720,1:1441] +
                           temp[1:721,1:1441]))


def compare_ccmp_to_buoy(   year = 2015,month=1,
                            ccmp_version = 'flk',
                            ccmp_subversion = 'V2.0',
                            spatial_subset='globe',
                            compare_buoys = 'RSS',
                            nrt = True,
                            plt_root = None,
                            nc_root = None,
                            save_maps = True, save_binned_plots = True, 
                            plot_maps = True, plot_binned_stats = True, 
                            save_binned_stats = True,
                            plot_1d_hists = True, save_1d_hists = True,save_1d_plots = True,
                            plot_hist_2d = True, save_2d_hist = True,save_2d_plots=True,
                            exclude_ccmp_with_sat = False, exclude_ccmp_no_sat = False):

    EN_wind = EN_wind_table()
    ccmp_wind_vars = ['w10','u10','v10']
    first_time_step = True
    num_days_in_month = calendar.monthrange(year,month)[1]
    #for day in range(1,num_days_in_month+1):
    for day in range(1,3):
                # read in CCMP with ERA5 background
        if nrt:
            version_path = f"P:/CCMP_RT/MEaSUREs/{ccmp_version}/{ccmp_subversion}/native/"
        else:
             version_path = f"P:/CCMP/MEaSUREs/{ccmp_version}/{ccmp_subversion}/native/"
            
        ccmp = load_ccmp_daily_arrays(year=year, 
                            month=month,
                            day=day,
                            verbose = True,
                            ccmp_path = version_path,
                            calc_wind_speed = True)

        for time_step in [0,1,2,3]:
            for wind_var in ccmp_wind_vars:
                date_str = f"{year:04}_{month:02}_{day:02}_{6*time_step:02}"
                print(date_str)
                
                d_cc = colocate_daily_wind_maps_to_buoys(wind_maps=ccmp,
                                                rain_map=None,
                                                current_map=None,
                                                buoy_set='RSS',
                                                max_time_diff=30,
                                                convert_buoy_to_neutral_wind=True,
                                                EN_wind=EN_wind,
                                                adjust_for_current=False,
                                                verbose=True)

            '''d_cc_navy = colocate_daily_wind_maps_to_buoys(wind_maps=ccmp_era5_cc,
                                                rain_map=None,
                                                current_map=None,
                                                buoy_set='NAVY',
                                                max_time_diff=30,
                                                convert_buoy_to_neutral_wind=True,
                                                EN_wind=EN_wind,
                                                adjust_for_current=False,
                                                verbose=True)'''

            '''d_era_cc = colocate_daily_wind_maps_to_buoys(wind_maps=era5_cc_maps,
                                                rain_map=None,
                                                current_map=None,
                                                buoy_set='RSS',
                                                max_time_diff=30,
                                                convert_buoy_to_neutral_wind=True,
                                                EN_wind=EN_wind,
                                                adjust_for_current=False,
                                                verbose=True)
            d_era = colocate_daily_wind_maps_to_buoys(wind_maps=era5_ns_maps,
                                                rain_map=None,
                                                current_map=None,
                                                buoy_set='RSS',
                                                max_time_diff=30,
                                                convert_buoy_to_neutral_wind=True,
                                                EN_wind=EN_wind,
                                                adjust_for_current=False,
                                                verbose=True)'''
            

compare_buoys = 'RSS'
ccmp_version = 'era5_current_corrected_adj2'
ccmp_subversion = 'v2.0'
compare_ccmp_to_buoy( year = 2016,month=1,
                            ccmp_version =ccmp_version,
                            ccmp_subversion = ccmp_subversion,
                            spatial_subset='globe',
                            compare_buoys = compare_buoys,
                            nrt = False,
                            plt_root =f'B:/job_CCMP/compare_ccmp_vs_{compare_buoys}_buoys/plots/',
                            nc_root = f'B:/job_CCMP/compare_ccmp_vs_{compare_buoys}_buoys/nc_files/',
                            save_maps = True, save_binned_plots = True, 
                            plot_maps = True, plot_binned_stats = True, 
                            save_binned_stats = True,
                            plot_1d_hists = True, save_1d_hists = True,save_1d_plots = True,
                            plot_hist_2d = True, save_2d_hist = True,save_2d_plots=True,
                            exclude_ccmp_with_sat = False, exclude_ccmp_no_sat = False)