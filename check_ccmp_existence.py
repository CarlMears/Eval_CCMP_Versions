import sys
sys.path.append('C:/job_CCMP/python/eval_CCMP_monthly_bias/')
sys.path.append('C:/job_CCMP/python/plotting_and_analysis/')
sys.path.append('C:/job_CCMP/python/binned_stats/')
sys.path.append('C:/job_CCMP/python/map_stats/')
sys.path.append('C:/job_CCMP/python/CCMP/')
sys.path.append('C:/job_CCMP/python/eval_era5/')
sys.path.append('C:/job_CCMP/python/era5/')
sys.path.append('C:/job_CCMP/python/buoys/')
sys.path.append('../')


from ccmp import read_ccmp_native
import calendar

version = 'era5_current_corrected_scl_01'
version_path = f"P:/CCMP/MEaSUREs/{version}/v2.0/native/"

for year in range(1994,2005):
    for month in range(1,13):
        num_days_in_month = calendar.monthrange(year,month)[1]
        for dom  in range(1,num_days_in_month+1):
            found = True
            for time_step in [0,1,2,3]:
                exists = read_ccmp_native(year=year,month=month,day=dom,time_step=time_step,path=version_path,check_exists=True)
                if not exists:
                    found = False
            if not found:
                print(f"Data MISSING for {month}/{dom}/{year}")
