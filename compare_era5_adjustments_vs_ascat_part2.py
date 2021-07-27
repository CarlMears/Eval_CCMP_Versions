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
sys.path.append('../')

from binned_stats import BinnedStat
from hist_2d import Hist2D

import copy
import os

plt_path  = 'C:/job_CCMP/compare_era5_vs_ASCAT/plots/binned_diffs/w/'

wind_vars = ['u','v','w']
wind_var_names = ['Zonal Wind','Mer. Wind','Wind Speed']

wind_vars = ['w','u']
wind_var_names = ['Wind Speed','Zonal Wind']

era5_versions       = ['10m_ns_oscar']
era5_versions_short = ['ERA5 NS Oscar']
for ivar, wind_var in enumerate(wind_vars):

    wind_var_name = wind_var_names[ivar]
    plt_path  = f'C:/job_CCMP/compare_era5_vs_ASCAT/plots/binned_diffs/{wind_var}/'
    try:
        os.makedirs(plt_path)
    except FileExistsError:
            # directory already exists
        pass

    for j,era5_version in enumerate(era5_versions):
        for adj in ['no_adj','adj']:
            era5_ver_sh = era5_versions_short[j] 
            for no_sat in [False]:
                for compare_sat in ['ascat_a']:
                    first_month = True
                    for  year in [2015,2016,2017,2018]:
                        for month in np.arange(1,13,dtype=np.int):
                            nc_file = f"C:/job_CCMP/compare_era5_vs_ASCAT/nc_files/{era5_version}/{wind_var}_era5_minus_{compare_sat}_binned_stats_{month:02}_{year:04}"
                            nc_file = nc_file + ".1hr"
                            if adj == 'adj':
                                nc_file = nc_file + ".adj"
                            nc_file = nc_file + ".nc"

                            w_era5_minus_ascat = BinnedStat.from_netcdf(nc_file=nc_file)
                            try:
                                assert(isinstance(w_era5_minus_ascat,BinnedStat))
                            except:
                                print()
                            if first_month:
                                w_era5_minus_ascat_tot = copy.copy(w_era5_minus_ascat)
                            else:
                                w_era5_minus_ascat_tot.combine(w_era5_minus_ascat)

                            nc_file = f"C:/job_CCMP/compare_era5_vs_ASCAT/nc_files/{era5_version}/{wind_var}_era5_{compare_sat}_hist_2d_{month:02}_{year:04}"
                            nc_file = nc_file + '.1hr'
                            if adj == 'adj':
                                nc_file = nc_file + '.adj'
                            nc_file = nc_file + '.nc'
                            print(nc_file)

                            w_era2_ascat_2dhist = Hist2D.from_netcdf(nc_file = nc_file,xname='',xunits='',yname='',yunits='')
                            if first_month:
                                w_era2_ascat_2dhist_tot = copy.copy(w_era2_ascat_2dhist)
                            else:
                                w_era2_ascat_2dhist_tot.combine(w_era2_ascat_2dhist)
                            first_month = False

                    
                    title_str = f"{wind_var_name}, era5 {era5_ver_sh} minus {compare_sat}, 2015-2018" 

                    if no_sat:
                        title_str = title_str + ', No Sat'

                    ytitle = f'ERA5 {wind_var_name}'
                    if adj == 'adj':
                        ytitle = f'ERA5 {wind_var_name}, Adjusted'
                    fig,ax = w_era2_ascat_2dhist_tot.plot(plot_diagonal=True,
                                                          title=title_str, 
                                                          xtitle=f'ASCAT {wind_var_name}', 
                                                          ytitle=f'ERA5 {wind_var_name}')
                    
                    title_str = f"{wind_var_name},{era5_ver_sh} minus {compare_sat}, 2015-2018"
                    if no_sat:
                        title_str = title_str + ', No Sat'
                    if wind_var == 'w':
                        xlab = f'Mean {wind_var_name} (ASCAT, era5)'
                    else:
                        xlab = f'ASCAT {wind_var_name}'
                    fig,ax = w_era5_minus_ascat_tot.plt(yrng=[-5.0,5.0],xlab=xlab,ylab='era5 - ASCAT',title=title_str)
                    print()
                    png_file = f"{plt_path}era5_{wind_var}_{era5_version}_minus_{compare_sat}_2015-2018"
                    if no_sat:
                        png_file = png_file + '.no_sat'
                    png_file = png_file + '.png'
                    fig.savefig(png_file)
                    stats = w_era5_minus_ascat_tot.calc_stats()
                    v_str = f"{wind_var_name}-{era5_ver_sh}-{compare_sat}"
                    if no_sat:
                        v_str = v_str + ' No Sat'
                    print(f"{v_str:20}: {stats['overall_num']:9.0f}, {stats['overall_mean']:6.3f}, {stats['overall_rms']:6.3f}, {stats['overall_stddev']:6.3f}")
plt.show()
print