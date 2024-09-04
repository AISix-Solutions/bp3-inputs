#GC 06/24
#Downloads downscaled CMIP data (precip, temp, wet days) from CanDCS project using PCIC servers. Aggregates data so that we get 1 file per year (with data for each day) by combinding results from
#11 of the PCIC models. These models are the PCIC 12 blend (a selection of the best performing models for Canada) minus UKESM1-0-LL as it uses a 360 day year. Some commenting in/out will be needed to
#switch between climate vars

import urllib.request
import xarray as xr
import numpy as np
import pandas as pd
import calendar
from time import sleep
import os
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

#PCIC 12 models and versions. Note that we drop UKESM1-0-LL as it uses a 360 day year. 
mod_ls=["BCC-CSM2-MR","NorESM2-LM","MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","EC-Earth3-Veg","CMCC-ESM2","INM-CM5-0","FGOALS-g3","TaiESM1","IPSL-CM6A-LR"]
vers_ls=["r1i1p1f1_gn","r1i1p1f1_gn","r1i1p1f2_gn","r1i1p1f1_gn","r1i1p1f1_gn","r1i1p1f1_gr","r1i1p1f1_gn","r1i1p1f1_gr1","r1i1p1f1_gn","r1i1p1f1_gn","r1i1p1f1_gr"]

for var in ["pr"]:#["pr","tasmax"]:
    
    for ssp in ["historical","ssp126","ssp245","ssp585"]:
        if ssp=="historical":
            ssp_url="ssp126"
            yr_ls=np.arange(2000,2015)
        else:
            ssp_url=ssp
            yr_ls=np.arange(2015,2101)
        for year in yr_ls:
            #Wet days are done differently from precip and temp as we get the total number of model-days with precip above 1 mm for that DOY. This will range from
            #0 (i.e., no models have precip > 1 mm on that DOY) to 11 (all models have precip > 1 mm). 

            #if os.path.isfile('D:\\CanDCSensemble\\'+ssp+'\\'+var+'\\'+str(year)+'_mean.nc'):  #For precip and temp
            if os.path.isfile('D:\\CanDCSensemble\\'+ssp+'\\wet_days\\'+str(year)+'_tot.nc'): #For wet days
                continue
            num_mod=0
            for mod, vers in zip(mod_ls, vers_ls):
                print(mod,vers,var,year,ssp)

                #for gregorian (i.e. leap days)
                if mod in ["MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","EC-Earth3-Veg","IPSL-CM6A-LR"]:
                    #Calc start and end dates for the given year
                    strt_date=str(year)+"-01-01"
                    delta=int((np.datetime64(strt_date)-np.datetime64('2000-01-01'))/np.timedelta64(1,'D'))
                    time_strt=18262+delta
            
                    end_date=str(year)+"-12-31"
                    delta=int((np.datetime64(end_date)-np.datetime64('2000-01-01'))/np.timedelta64(1,'D'))
                    time_end=18262+delta

                    #Define download URL
                    download_url="https://services.pacificclimate.org/data/downscaled_cmip6_multi/"+var+"_day_MBCn+PCIC-Blend_"+mod+"_historical+"+ssp_url+"_"+vers+"_19500101-21001231.nc.nc?"+var+"["+str(time_strt)+":"+str(time_end)+"][0:510][0:1068]&"
                    
                    success=0
                    attempt=0
                    while attempt<=5 and success==0:
                        try:
                            urllib.request.urlretrieve(download_url,"C:\\Users\\GiovanniCorti\\Downloads\\"+mod+"_temp.nc")
                            success=1
                        except Exception as e:
                            print(e)
                            print("---------------------------------")
                            attempt+=1
                            print("Download failed. Trying again. Attempt "+str(attempt)+" of 5")
                            sleep(30)
                    
                    if success==0:
                        break
                    
                    #Wrangle time to drop leap day
                    nc_data=xr.open_dataset("C:\\Users\\GiovanniCorti\\Downloads\\"+mod+"_temp.nc")
                    nc_data=nc_data.where(nc_data<1,1).where(nc_data>1,0)
                    if calendar.isleap(year):
                        nc_data=nc_data.drop([np.datetime64(str(year)+'-02-29T12:00:00.000000000')],dim="time")
                    #print(nc_data["time"])
                    #nc_data["time"]=pd.to_datetime(nc_data["time"])

                    

                #for 365 day calender (i.e. no leap days)
                else: 
                    #Calc start and end dates for the given year
                    time_strt=365*(year-2000)+18250
                    time_end=time_strt+364
                    
                    #Define download URL
                    download_url="https://services.pacificclimate.org/data/downscaled_cmip6_multi/"+var+"_day_MBCn+PCIC-Blend_"+mod+"_historical+"+ssp_url+"_"+vers+"_19500101-21001231.nc.nc?"+var+"["+str(time_strt)+":"+str(time_end)+"][0:510][0:1068]&"
                    success=0
                    attempt=0
                    while attempt<=5 and success==0:
                        try:
                            urllib.request.urlretrieve(download_url,"C:\\Users\\GiovanniCorti\\Downloads\\"+mod+"_temp.nc")
                            success=1
                        except Exception as e:
                            print(e)
                            print("---------------------------------")
                            attempt+=1
                            print("Download failed. Trying again. Attempt "+str(attempt)+" of 5")
                            sleep(30)
                    
                    if success==0:
                        break
                    
                    #Wrangle time
                    nc_data=xr.open_dataset("C:\\Users\\GiovanniCorti\\Downloads\\"+mod+"_temp.nc",decode_times=False)
                    nc_data=nc_data.where(nc_data<1,1).where(nc_data>1,0) #For wet days only
                    
                    
                    time_ls=pd.date_range(str(year)+'-01-01',str(year)+'-12-31',freq='d')
                    if calendar.isleap(year):
                        time_ls=time_ls[time_ls!=np.datetime64(str(year)+'-02-29')]
                    nc_data["time"]=time_ls+np.timedelta64(12, 'h')
                    #nc_data["time"]=np.datetime64(str(year)+'-01-01T00:00:00.000000')+nc_data["time"].values.astype("timedelta64[D]")+np.timedelta64(12, 'h')
                
                if num_mod==0:
                    dsum=nc_data
                    #dsum2=nc_data**2  #For precip and temp, used to calc stds
                else:
                    dsum+=nc_data
                    #dsum2+=nc_data.astype('float32')**2  #For precip and temp, used to calc stds
                nc_data.close()
                num_mod+=1
            #avg=dsum/num_mod #For precip and temp, used to calc stds
            #std=np.sqrt((dsum2/num_mod)-(avg**2)) #For precip and temp, used to calc stds

            #os.makedirs('D:\\CanDCSensemble\\'+ssp+'\\'+var, exist_ok=True)
            os.makedirs('D:\\CanDCSensemble\\'+ssp+'\\wet_days', exist_ok=True)
            if num_mod==11:
                dsum.to_netcdf('D:\\CanDCSensemble\\'+ssp+'\\wet_days\\'+str(year)+'_tot.nc',engine='scipy') #For wet days
                #std.to_netcdf('D:\\CanDCSensemble\\'+ssp+'\\'+var+'\\'+str(year)+'_std.nc') #For precip and temp
                #avg.to_netcdf('D:\\CanDCSensemble\\'+ssp+'\\'+var+'\\'+str(year)+'_mean.nc') #For precip and temp
            else:
                print("Not enough models downloaded...")
        

    
    