#This file contains code for tuning SEDs for each ecozone in order to get
#fire sizes that match real-world results. 
#GC -- 05/24

import pysyncrosim as ps
import pandas as pd
import os
import rasterio as rio
import time
import rioxarray as rxr
import numpy as np
import subprocess
import multiprocessing as mp
from scipy.optimize import fsolve
import math
import geopandas as gpd


def run_BP3(NTS_code, it_list):
    """
    Wrapper function to setup and run NTS sheets. Rescales burn probability
    
    Parameters
    ----------
    NTS_code: string 
        NTS code for the desired sheet
    
    it_list: list of integers
        List of iteration numbers
    
    Returns
    -------
    """
    
    #Setup scenario
    ign_rescaling=_setup_scen(NTS_code,"C:\\Users\\GiovanniCorti\\Documents\\NTS_template.ssim")
    
    #Check if zero igns
    if ign_rescaling==None:
        return
    
    for it_num in it_list:
        #Run BP3 for all iteration values
        bp_dst_fp, it_num=_init_run(NTS_code,it_num,"C:\\Users\\GiovanniCorti\\Documents\\NTS_template.ssim")
        if ign_rescaling<1: #Rescale igns if needed. 
            raster=rio.open(bp_dst_fp)
            rd=raster.read(1)*ign_rescaling
            raster.close()
            with rio.open(bp_dst_fp, "w", **raster.meta) as dst:
                dst.write(rd,1)
                
def _setup_scen(NTS_code,ssim_fp,sess_fp='C:/Program Files/SyncroSim'):
    """
    Pulls in the NTS sheet specific parameters for a base Syncrosim scenario 
    (.ssim file) from the the the network drive (\\athena02) and saves the .ssim file. 

    Parameters
    ----------
    NTS_code : string
        Code/name for the NTS sheet that is being modified.
    ssim_fp : string
        File path for the .ssim file.
    sess_fp : string, optional
        Filepath for the Syncrosim install. The default is 
        'C:/Program Files/SyncroSim'.

    Returns
    -------
    ign_rescaling : float or None  
        rescaling factor needed for average ign numbers 
        less than 1. Returns none if the sheet has 0 igns. 

    """
    
    sess=ps.Session(sess_fp, silent=False)
    #Read in a library, project and Scenario that has already been created.
    #This scenario will then be modified in the loop to change the number of iterations.
    lib=ps.library(name = ssim_fp, session=sess, package='burnP3Plus', addons='burnP3PlusCell2Fire', use_conda=True)
    proj=lib.projects(pid=1)
    scen=proj.scenarios(sid=1)
    
    #Set landscape rasters
    ls_ds=scen.datasheets(name='burnP3Plus_LandscapeRasters')
    Fuel_fp='Y:/client-data/demo_projects/climate85/Working_data/NARR_weather_csvs/NTS_SNRC_'+NTS_code+'/Fuel_NTS_SNRC_'+NTS_code+'.tif'
    DEM_fp='Y:/client-data/demo_projects/climate85/Working_data/NARR_weather_csvs/NTS_SNRC_'+NTS_code+'/DEM_NTS_SNRC_'+NTS_code+'.tif'
    ls_ds.loc[0,['ElevationGridFileName','FuelGridFileName']]=[DEM_fp,Fuel_fp]
    try: #May fail if no fuel map exists (i.e., sheet is all water/Arctic)
        scen.save_datasheet(name='burnP3Plus_LandscapeRasters',data=ls_ds)
    except:
        print("No fuel map found. Moving to next NTS sheet") 
        return None
    #Set ignition number/distribution
    ign_dist=pd.read_csv("Y:/client-data/demo_projects/climate85/Working_data/ign_dist/ign_dist_"+NTS_code+".csv")
    
    #Needed to check if average ign number is less than 1 as BP3+ will not accept zero as an ign value. 
    ign_num=(ign_dist['ign_per_it']*(ign_dist['pct']/100)).sum()
   
    if ign_num>=1: #Normal scenario
        #Create ign distribution
        ign_dist['Name']='Igns'
        ign_dist=ign_dist.rename(columns={'ign_per_it': 'Value', 'pct': 'RelativeFrequency'})
        ign_rescaling=1
    elif ign_num==0: #No igns, exits function
        print('No Ignitions for this sheet.')
        return None
    else: #Sheets with ignition number between 0 and 1. Here we just run w/ 1 ign and then rescale BPs. 
        #A distribution (w/ 100% chance of 1) is used to keep formatting consistent
        val_list=[1]
        rf_list=[100]
        dict = {'Value': val_list, 'RelativeFrequency': rf_list} 
        ign_dist=pd.DataFrame(dict)
        ign_dist['Name']='Igns'
        ign_rescaling=ign_num

    
    #Check for probablistic ign grids
    file_pre = ["H_Spring_", "H_Summer_", "L_Spring_","L_Summer_"]
    ign_grid_fps=[]
    for f in file_pre:
        ign_grid_fps.append('Y:/client-data/demo_projects/climate85/Working_data/NARR_weather_csvs/NTS_SNRC_'
        +NTS_code+'/'+f+'NTS_SNRC_'+NTS_code+'.tif')
        
    ign_grid_ds=scen.datasheets(name="burnP3Plus_ProbabilisticIgnitionLocation")
    if all(list(map(os.path.isfile,ign_grid_fps))):
        #Set ign grids
        ign_grid_ds.loc[0,['Season','Cause','IgnitionGridFileName']]=[1,'Human',ign_grid_fps[0]]
        ign_grid_ds.loc[1,['Season','Cause','IgnitionGridFileName']]=[2,'Human',ign_grid_fps[1]]
        ign_grid_ds.loc[2,['Season','Cause','IgnitionGridFileName']]=[1,'Lightning',ign_grid_fps[2]]
        ign_grid_ds.loc[3,['Season','Cause','IgnitionGridFileName']]=[2,'Lightning',ign_grid_fps[3]]
        scen.save_datasheet(name="burnP3Plus_ProbabilisticIgnitionLocation",data=ign_grid_ds)
    else:
        #If no ign grids, ensure that the ign grid values are empty
        scen.save_datasheet(name="burnP3Plus_ProbabilisticIgnitionLocation",data=ign_grid_ds[0:0])
    
    #Set weather stream
    FWI_ds=scen.datasheets(name="burnP3Plus_WeatherStream")
    FWI=pd.read_csv('Y:/client-data/demo_projects/climate85/Working_data/NARR_weather_csvs/NTS_SNRC_'+NTS_code+'/FWI_NTS_SNRC_'+NTS_code+'.csv')
    cnames=['Season', 'Temperature', 'RelativeHumidity', 'WindSpeed', 'WindDirection', 
            'Precipitation', 'FineFuelMoistureCode', 'DuffMoistureCode', 'DroughtCode', 
            'InitialSpreadIndex', 'BuildupIndex', 'FireWeatherIndex']
    FWI.columns=cnames
    scen.save_datasheet(name='burnP3Plus_WeatherStream',data=FWI)

    
    #Get SED from csv set in SyncroSim
    SED=pd.read_csv("Y:/client-data/demo_projects/climate85/Working_data/sed_dist/sed_dist_"+NTS_code+".csv")
    SED['Name']='SED'
    SED=SED.rename(columns={'sp_ev_days': 'Value', 'pct': 'RelativeFrequency'})
    SED=pd.concat([SED, ign_dist], ignore_index=True)
    scen.save_datasheet(name="burnP3Plus_DistributionValue",data=SED)
    print('Setup complete for NTS sheet '+NTS_code)

    return ign_rescaling, scen

def _init_run(NTS_code,it_num,ssim_fp,sess_fp='C:/Program Files/SyncroSim'):
    '''
    Runs BP3+, with given number(s) of iterations, for an NTS sheet and saves 
    outputs to the C:\\BP3IO folder.
    
    Parameters
    ----------
    NTS_code : string
        Code/name for the NTS sheet being run.
    it_num : Integer
        Number of iteration to run.
    ssim_fp : string
        File path for the .ssim file.
    sess_fp : string, optional
        Filepath for the Syncrosim install. The default is 
        'C:/Program Files/SyncroSim'.

    Returns
    -------
    bp_dst_fp: String
        File path for the burn probability file
    
    it_num: Integer
        Number of iterations run

    '''
    #Connect to syncrosim
    sess=ps.Session(sess_fp, silent=False)
    #Read in a library, project and Scenario that has already been created.
    #This scenario will then be modified in the loop to change the number of 
    #iterations. Scenario here corresponds to Ft. Mac NTS sheet
    lib=ps.library(name = ssim_fp, session=sess, package='burnP3Plus', addons='burnP3PlusCell2Fire', use_conda=True)
    proj=lib.projects(pid=1)
    scen=proj.scenarios(sid=1)
    
    #Running scenario for it_num iterations
    df_time=pd.DataFrame(columns=['Iteration Number','Time (min)'])
    os.makedirs('C:\\Users\\GiovanniCorti\\Documents\\Stats', exist_ok=True)
    #os.makedirs('BP3IO/'+NTS_code+'/BurnMap', exist_ok=True)
    #os.makedirs('BP3IO/'+NTS_code+'/Stats', exist_ok=True)
    
    #Change iteration number
    RC=scen.datasheets(name="burnP3Plus_RunControl")
    RC.at[0,'MaximumIteration']=it_num
    scen.save_datasheet(name="burnP3Plus_RunControl",data=RC)
                 
    #Run simulation
    print('Running NTS Sheet '+NTS_code+ ' for '+ str(it_num)+ ' iterations.')
    st_time=time.time()

    out=scen.run(jobs=mp.cpu_count()-1)
        
    #Copy burn probability maps        
    bp_fp="C:\\NTS_template.ssim.temp\\summary\\"+str(out.datasheets(name="burnP3Plus_OutputBurnProbability")["FileName"][0])
    bp_dst_fp='C:\\Users\\GiovanniCorti\\Documents\\BP_maps\\BP_it'+str(it_num)+'.tif'
    os.system('copy /y '+bp_fp+' '+bp_dst_fp)
            
    #Output burn count
    #bc_fp="C:\\NTS_template.ssim.temp\\summary\\"+str(out.datasheets(name="burnP3Plus_OutputBurnCount")["FileName"][0])
    #bc_dst_fp=r"C:\\BP3IO\\"+NTS_code+'\\BurnMap\\BM_it'+str(it_num)+'.tif'
    #os.system('copy /y '+bc_fp+' '+bc_dst_fp)
            
    #Burn stats spreadsheet
    stat_tab=out.datasheets(name="burnP3Plus_OutputFireStatistic")
    stat_tab.to_csv("C:\\Users\\GiovanniCorti\\Documents\\Stats\\FireStats_"+NTS_code+"_it"+str(it_num)+'.csv')
        

    end_time=time.time()
    tt=(end_time-st_time)/60
    df_time.loc[len(df_time.index)] = [it_num,tt]
    print('Done running '+str(it_num)+' iterations.')
        
    #df_time.to_csv('C:/BP3IO/'+NTS_code+'/run_time.csv')
    bp_dst_fp=None
    return bp_dst_fp, it_num

#Get probs for zero-truncated poisson dist dists
def ztp(L,x):
    return np.float64(((L**x)*np.e**(-1*L))/(math.factorial(x)*(1-np.e**(-1*L))))

#Calc parameters to define ztp based on desired average. Will be fed into 
#fsolve, a numerical solver
def ztp_mean(L,mu):
    return (L*np.e**(L))/(np.e**(L)-1)-mu    

#Create ztp dist
def ztp_dist(name,mu):
    assert name in ['SED','Igns'], "Name must be either SED or Igns"
    rf_ls=[]
    val_ls=[]
    #Calc lambda param
    L=fsolve(ztp_mean,0.01,args=mu)
    x=1
    
    #Calc prob until prob is less that 1 percent
    prob=ztp(L,x)
    if mu==1:
        rf_ls.append(1)
        val_ls.append(1)
    else:
        rf_ls.append(prob)
        val_ls.append(x)
        while prob>.01 or x<mu:
            x=x+1
            prob=ztp(L,x)
            rf_ls.append(prob)
            val_ls.append(x)
    ztp_df=pd.DataFrame({'Value':val_ls,'RelativeFrequency':rf_ls})
    ztp_df['Name']=name
    return ztp_df

def run_test_nts(NTS_ls, SED_mu,ez_fs_mu):
    """
    Function to setup and run NTS sheets, adjusting the number of SEDs
    until the average. Here we use the ign values and probabilistic ign grids
    stored on the Y: drive
    
    Parameters
    ----------
    NTS_ls: list of string 
        NTS codes for the desired sheet
    
    SED_mu: float
        Average SED number to try
        
    ez_fs_mu: float
        Target average fire size
            
    Returns
    ---------
    
    SED_mu: float
        Final SED number attempted. This will yield a result close to the
        desired fire size
    
    """
    fs_delta=1
    while fs_delta>.05:
        print("Trying with SED avg of "+str(SED_mu))
        area_sum=0
        tot_fires=0
        for NTS_code in NTS_ls:
            #Setup NTS_sheet
            ign_rescaling,scen=_setup_scen(NTS_code,"C:\\Users\\GiovanniCorti\\Documents\\NTS_template.ssim")
            
            #Create ztp dist based of SED_mu
            SED_df=ztp_dist('SED',SED_mu)
            ign_df=pd.read_csv("Y:/client-data/demo_projects/climate85/Working_data/ign_dist/ign_dist_"+NTS_code+".csv")
            ign_df['Name']='Igns'
            ign_df=ign_df.rename(columns={'ign_per_it': 'Value', 'pct': 'RelativeFrequency'})
            
            dist_df=pd.concat([SED_df, ign_df], ignore_index=True)
            scen.save_datasheet(name="burnP3Plus_DistributionValue",data=dist_df)
            
            #Run NTS sheet for 500 its
            bp_dst_fp, it_num=_init_run(NTS_code,500,"C:\\Users\\GiovanniCorti\\Documents\\NTS_template.ssim")
            
            #Read stats csv and calc params for average fire size
            stats_df=pd.read_csv("C:\\Users\\GiovanniCorti\\Documents\\Stats\\FireStats_"+NTS_code+"_it"+str(it_num)+'.csv')
            area_sum=area_sum+stats_df["Area"].sum()
            tot_fires=tot_fires+len(stats_df["Area"])
            
        #Calc avg fire size and delta from desired target size
        mod_fs_mu=area_sum/tot_fires
        fs_delta=np.abs((mod_fs_mu-ez_fs_mu)/ez_fs_mu)
        print("Ecozone Avg fire size",ez_fs_mu)
        print("Model Avg fire size",mod_fs_mu)
        
        #Stop if SED falls below 1
        if SED_mu==1 and ez_fs_mu<mod_fs_mu:
            print('Hit min SED val. Stopping')
            exit()
            
        #Adjust SED_mu up or down depending on how far off we are from the
        #desired ecozone fire size. Step size varies based on how far off
        #we are from the target
        if ez_fs_mu>mod_fs_mu:
            SED_mu=SED_mu+((SED_mu*fs_delta)/2)
        if ez_fs_mu<mod_fs_mu:
            SED_mu=SED_mu-((SED_mu*fs_delta)/2)
            if SED_mu<1:
                SED_mu=1
    return SED_mu
        
def run_EZs():
    #Wrapper function to run run_test_nts for all ecozones
    
    #Test sheets for each ecozone. Here I attempt span the geographic extent of
    #the ecozone with a few sheets.  
    ts_dict={4.0:["084M","095P","106N"], 5.0:["086B","065D","023N"],
             6.1:["012L","032H","042C"], 6.2:["052H","053L","074K"],9.0:["094A","073M","063B"],
             11.0:["105P","116J"],12.0:["094L","115I"],13.0:["092K","103I"],
             14.0:["082E","083D","093N"],15.0:["054F","042O"]}
    #Target average fire size for each ecozone. 6.1 is BSE and 6.2 is BSW
    fs_dict={4.0: 2824.47, 5.0: 3740.222, 6.2: 3700.776, 6.1: 2783.094,
             9.0: 2121.775, 11.0: 2958.133, 
             12.0: 4296.312, 13.0: 465.003, 14.0: 1679.746, 15.0: 1254.873}
    
    tuned_sed_df=pd.DataFrame(columns=["Ecozone Code", "SED value"])
    #Run tunning for each ecozone
    for ez_code in fs_dict:
        print("Running Ecozone", ez_code)
        SED_mu=run_test_nts(ts_dict[ez_code],2,fs_dict[ez_code])
        tdf=pd.DataFrame({"Ecozone Code":[ez_code], "SED value":[SED_mu]})
        
        #Saved tuned SEDs to a csv
        tuned_sed_df=pd.concat([tuned_sed_df,tdf])
        tuned_sed_df.to_csv("C:\\Users\\GiovanniCorti\\Documents\\tuned_sed.csv")
        
    
    
if __name__ == "__main__":
    

    #Connect to Y: drive for NTS sheet parameter files
    subprocess.call(r"net use y: \\192.168.99.12\shared /user:gcorti 850Whastings")
    run_EZs()    
