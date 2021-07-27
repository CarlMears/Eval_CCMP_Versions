import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import calendar
import copy
import os

import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/python_packages_cam/rss_stats/')
sys.path.append('C:/python_packages_cam/')
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
sys.path.append('../')

from binned_stats import BinnedStat

from hist_2d import Hist2D
from hist_1d import Hist1D
from plot_2d_hist import plot_2d_hist,averages_from_histograms
import copy
import os


compare_sat = 'WindSAT'   
ccmp_version = 'era5_current_corrected_adj2'
version_sh = 'v2.0'

var_list = ['w']
var_name_list = ['Wind Speed']
for iw,wind_var in enumerate(var_list):
    wind_var_name = var_name_list[iw]

    start_year=2013
    end_year = 2019
    year_list = np.arange(start_year,end_year+1)
    verbose = True
    rad_rain = False
    for ccmp_subversion in ['v2.0']:
        for no_sat in [False,True]:
            first_month = True
            for  year in year_list:
                for month in np.arange(1,13,dtype=np.int):
                    #binned differences
                    nc_file = f"C:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/{ccmp_version}_{ccmp_subversion}/{wind_var}_ccmp_minus_{compare_sat}_binned_stats_{month:02}_{year:04}"
                    if no_sat:
                        nc_file = nc_file + "_no_sat"
                    if rad_rain:
                        nc_file = nc_file + "_radrain"
                    nc_file = nc_file + ".nc"
                    if verbose:
                        print(nc_file)
                    try:
                        da = xr.open_dataarray(nc_file)
                        w_ccmp_minus_ascat = BinnedStat.from_DataArray(da)
                    except:
                        w_ccmp_minus_ascat = BinnedStat.from_netcdf(nc_file)
                    print(f"{ccmp_subversion} {no_sat} {month:02d}/{year:04d}  {w_ccmp_minus_ascat.overall_num}")
                    if first_month:
                        w_ccmp_minus_ascat_tot = copy.copy(w_ccmp_minus_ascat)
                        
                    else:
                        w_ccmp_minus_ascat_tot.combine(w_ccmp_minus_ascat)

                    #2-D histograms
                    nc_file = f"C:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/{ccmp_version}_{ccmp_subversion}/CCMP_{compare_sat}_2d_hists_{month:02}_{year:04}"
                    if no_sat:
                        nc_file = nc_file + "_no_sat"
                    if rad_rain:
                        nc_file = nc_file + "_radrain"
                    nc_file = nc_file + ".nc"
                    if verbose:
                        print(nc_file)
                    
                    w_hist_2D = Hist2D.from_netcdf(nc_file = nc_file,var = wind_var,
                                                    compare_sat=compare_sat,
                                                    xname=f"{wind_var.upper()} {compare_sat}",xunits="m/s",
                                                    yname=f"{wind_var.upper()} CCMP",yunits="m/s")
                    if first_month:
                        w_hist_2D_tot = copy.copy(w_hist_2D)
                    else:
                        w_hist_2D_tot.combine(w_hist_2D)

                    #1_d histograms
                    nc_file = f"C:/job_CCMP/compare_ccmp_vs_{compare_sat}/nc_files/{ccmp_version}_{ccmp_subversion}/CCMP_{compare_sat}_1d_hists_{month:02}_{year:04}"
                    if no_sat:
                        nc_file = nc_file + "_no_sat"
                    if  rad_rain:
                        nc_file = nc_file + "_radrain"
                    nc_file = nc_file + ".nc"
                    if verbose:
                        print(nc_file)
                    w_hist_1D_sat  = Hist1D.from_netcdf(nc_file=nc_file, varname = f'n_{wind_var}_{compare_sat}', name = f"{wind_var.upper()} {compare_sat}")
                    w_hist_1D_ccmp = Hist1D.from_netcdf(nc_file=nc_file, varname = f'n_{wind_var}_ccmp{compare_sat}', name = f"{wind_var.upper()} CCMP at {compare_sat}")
                    
                    if first_month:
                        w_hist_1D_sat_tot = copy.copy(w_hist_1D_sat)
                        w_hist_1D_ccmp_tot = copy.copy(w_hist_1D_ccmp)
                    else:
                        w_hist_1D_sat_tot.combine(w_hist_1D_sat)
                        w_hist_1D_ccmp_tot.combine(w_hist_1D_ccmp)
                    first_month = False

            title_str = f"{wind_var.upper()}, CCMP {ccmp_version}_{ccmp_subversion} minus {compare_sat}, {start_year}-{end_year}"
            if no_sat:
                title_str = title_str + ', No Sat'
            if rad_rain:
                title_str = title_str + ', RadRain'
            xlab = f'Mean {wind_var_name} (ASCAT, CCMP)'

            fig,ax = w_ccmp_minus_ascat_tot.plot(yrng=[-5.0,5.0], xlab=xlab, ylab=f'{wind_var.upper()}, CCMP - ASCAT', title=title_str,fontsize=12)
            
            plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/binned_mean_summaries/'
            os.makedirs(plt_path,exist_ok = True)
            png_file = f"{plt_path}CCMP_{wind_var}_{ccmp_version}_{ccmp_subversion}_minus_{compare_sat}_{start_year}-{end_year}"
            if no_sat:
                png_file = png_file + '.no_sat'
            if rad_rain:
                png_file = png_file  + '.radrain'
            png_file = png_file + '.png'
            fig.savefig(png_file)
            title_str = f"CCMP {ccmp_version}_{ccmp_subversion} vs {compare_sat}, {start_year}-{end_year}"
            if no_sat:
                title_str = title_str + ', No Sat'
            
            fig,ax = w_hist_2D_tot.plot(title=title_str,xtitle=f'{wind_var.upper()} {compare_sat} (m/s)', ytitle=f'{wind_var.upper()} CCMP (m/s)',fontsize=12)
            
            plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/hist_2D_summaries/'
            os.makedirs(plt_path,exist_ok = True)
            png_file = f"{plt_path}CCMP_{wind_var}_{ccmp_version}_{ccmp_subversion}_minus_{compare_sat}_{start_year}-{end_year}"
            if no_sat:
                png_file = png_file + '.no_sat'
            if rad_rain:
                png_file = png_file  + '.radrain'
            png_file = png_file + '.png'
            fig.savefig(png_file)

            fig,ax = w_hist_1D_sat_tot.plot(title=title_str,fontsize=12)
            fig,ax = w_hist_1D_ccmp_tot.plot(fig=fig,ax=ax,fontsize=12)

            figlog,axlog = w_hist_1D_sat_tot.plot(semilog=True,title=title_str,fontsize=12)
            figlog,axlog = w_hist_1D_ccmp_tot.plot(fig=figlog,ax=axlog,semilog=True,fontsize=12)

            plt_path  = f'C:/job_CCMP/compare_ccmp_vs_{compare_sat}/plots/{ccmp_version}_{ccmp_subversion}/hist_1D_summaries/'
            os.makedirs(plt_path,exist_ok = True)
            png_file = f"{plt_path}CCMP_{wind_var}_{ccmp_version}_{ccmp_subversion}_minus_{compare_sat}_{start_year}-{end_year}"
            if no_sat:
                png_file = png_file + '.no_sat'
            if rad_rain:
                png_file = png_file  + '.radrain'
            png_file = png_file + '.png'

            fig.savefig(png_file)
            png_file = f"{plt_path}CCMP_{wind_var}_{ccmp_version}_{ccmp_subversion}_minus_{compare_sat}_{start_year}-{end_year}"
            if no_sat:
                png_file = png_file + '.no_sat'
            if rad_rain:
                png_file = png_file  + '.radrain'
            png_file = png_file + '.log.png'
            figlog.savefig(png_file)

            stats = w_ccmp_minus_ascat_tot.calc_stats()
            v_str = f"{wind_var_name}-{ccmp_version}_{ccmp_subversion}-{compare_sat}"
            if no_sat:
                v_str = v_str + ' No Sat'
            if rad_rain:
                v_Str =v_str + ' Rad Rain'
            print(f"{v_str:30}: {stats['overall_num']:9.0f}, {stats['overall_mean']:6.3f}, {stats['overall_rms']:6.3f}, {stats['overall_stddev']:6.3f}")
plt.show()
print