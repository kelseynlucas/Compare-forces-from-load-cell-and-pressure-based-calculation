# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 17:13:39 2016

@author: Kelsey


These functions will compare pressure code and Flapper data, using the Flapper 
data as predicted values and the pressure code as measured.  First, the two will
be compared with a cross-correlation, which will show that the two time series 
have the same shape and no (limited?) phase shifts.  Also find the normalized 
correlation coefficent for standard comparison.  Next, root mean square error
(RMSE) will be used to compare the amplitudes of the signals and determine 
approximately how far off the pressure code is from the Flapper values.


"""

#Import useful packages
import pandas as pd
import numpy as np
import sklearn
import matplotlib.pyplot as plt


def rmse(predicted, measured):
    """
    calculates the RMSE for 2 datasets
    """
    resids = []    
    for i in range(0,len(predicted)):
        resids.append((predicted[i]-measured[i])**2)
    
    return (sum(resids)/float(len(resids)))**0.5







#load useful file operations
from os import listdir
from os.path import isfile, join
 
 
def get_paired_files(test_type):
    """
    Reads a test type - and finds the pairs of flapper and pcode datasets for
        all the parent vids in that test type.  Returns a dictionary of the filepaths 
        for these pairs
     
    Input
     
    -test type - name of the test type
         
    """
     
    #If have a subfolder (not a file), get name and add to list of subfolders
    parentVids = [f for f in listdir('G:/Nonunif_performance_PIV/Pressure code analysis/' + test_type) if not isfile(join('G:/Nonunif_performance_PIV/Pressure code analysis/' + test_type,f))]
 
    #Initialize storage for full stack folders
    pairs = {}
 
    #for each parent vid, get a list of the PIV data folders
    for vid in parentVids:
        #set filepath to the video's folder        
        
        if 'tail' in vid:
            pcode_path = 'G:/Nonunif_performance_PIV/Pressure code analysis/' + test_type + '/' + vid + '/results_3cyc/queen2_results_fulltimetrace_fixed_axisfix.xlsx'
        else:
            pcode_path = 'G:/Nonunif_performance_PIV/Pressure code analysis/' + test_type + '/' + vid + '/results_3cyc/queen2_results_fulltimetrace_axisfix.xlsx'
        
        flapper_path = 'G:/Nonunif_performance_PIV/Data_fromFlapper_analyzed/final_analysis_7Hzfiltercutoff/' + test_type + '/' + vid + '_downsampled_resfix.xlsx'
         
        #add the path to the dictionary
        pairs[vid] = [flapper_path,pcode_path]
 
    #return the dictionary of PIV folders
    return pairs




def compare(Flapper_path,Pcode_path):
    """
    Reads in corresponding Flapper and Pcode data, does a cross-correlation to
    find the best phase lag, reports the normalized correlation coefficient for 
    phase lag = 0, and finds the RMSE%.
    """
    
    

    #Read in pressure code data
    Pcode = pd.ExcelFile(Pcode_path)
    Pcode = Pcode.parse('Sheet1', index_col=None)
    
    #read in processed, downsampled Flapper data
    Flapper = pd.ExcelFile(Flapper_path)
    Flapper = Flapper.parse('Sheet1', index_col=None)
    
    #Below not needed after axis fix    
    '''
    #correct for pressure sensor (+)Fy definition change
    if 'heave' in Pcode_path:
        Flapper['Fy'] = -Flapper['Fy']
    if 'tail' in Pcode_path:
        Pcode['PTz'] = -Pcode['PTz']
    '''
    
    #Begin normalizing the time traces by subtracting the mean
    Flapper['nFx'] = Flapper['Fx'] - np.mean(Flapper['Fx'])
    Flapper['nFy'] = Flapper['Fy'] - np.mean(Flapper['Fy'])
    Flapper['nTz'] = Flapper['Tz'] - np.mean(Flapper['Tz'])
    
    Pcode['nPFx'] = Pcode['PFx'] - np.mean(Pcode['PFx'])
    Pcode['nPFy'] = Pcode['PFy'] - np.mean(Pcode['PFy'])
    Pcode['nPTz'] = Pcode['PTz'] - np.mean(Pcode['PTz'])
    
    
    
    #Generate normalization coefficient for correlation factor "std(A)*std(B)*(n-1)
    factX = np.sqrt(np.sum((Flapper['Fx']-np.mean(Flapper['Fx']))**2.)/(len(Pcode['PFx'])-1.))*np.sqrt(np.sum((Pcode['PFx']-np.mean(Pcode['PFx']))**2.)/(len(Pcode['PFx'])-1.))*(len(Flapper['Fx'])-1.)
    factY = np.sqrt(np.sum((Flapper['Fy']-np.mean(Flapper['Fy']))**2.)/(len(Pcode['PFy'])-1.))*np.sqrt(np.sum((Pcode['PFy']-np.mean(Pcode['PFy']))**2.)/(len(Pcode['PFy'])-1.))*(len(Flapper['Fy'])-1.)
    factZ = np.sqrt(np.sum((Flapper['Tz']-np.mean(Flapper['Tz']))**2.)/(len(Pcode['PTz'])-1.))*np.sqrt(np.sum((Pcode['PTz']-np.mean(Pcode['PTz']))**2.)/(len(Pcode['PTz'])-1.))*(len(Flapper['Tz'])-1.)
    
    #Create correlation factor sequence for all phase lags
    #Finish normalization by dividing by normalization coefficient
    corr_seq_Fx=np.correlate(Flapper['nFx'], Pcode['nPFx'], 'full')/factX
    corr_seq_Fy=np.correlate(Flapper['nFy'], Pcode['nPFy'], 'full')/factY
    corr_seq_Tz=np.correlate(Flapper['nTz'], Pcode['nPTz'], 'full')/factZ
    
    #Create a vector with the lags (shift number of spaces)
    n=(len(corr_seq_Fx)-1)/2
    lags = np.arange(-n,n+1)
    
    #Find the list index of the largest correlation coefficient in each case
    lag_Fx=np.argmax(corr_seq_Fx)
    lag_Fy=np.argmax(corr_seq_Fy)
    lag_Tz=np.argmax(corr_seq_Tz)
    
    #Find the largest correlation coefficient
    coef_Fx=corr_seq_Fx[lag_Fx]
    coef_Fy=corr_seq_Fy[lag_Fy]
    coef_Tz=corr_seq_Tz[lag_Tz]
    
    """
    #print the max coefficients and the associated lag
    print coef_Fx, lags[lag_Fx]
    print coef_Fy, lags[lag_Fy]
    print coef_Tz, lags[lag_Tz]
    """
    if 'tail' in Pcode_path:
        corr_path = Pcode_path.split('queen2')[0] + 'xcorr_seq_fixed_axisfix.xlsx'
    else:
        corr_path = Pcode_path.split('queen2')[0] + 'xcorr_seq_resfix.xlsx'
    corr_df = pd.DataFrame(lags,columns=['lags'])
    corr_df['Fx'] = corr_seq_Fx    
    corr_df['Fy'] = corr_seq_Fy
    corr_df['Tz'] = corr_seq_Tz
    corr_df.to_excel(corr_path)
    
    #Get the normalized correlation coeff.
    ncoef_Fx = np.corrcoef(Flapper['Fx'],Pcode['PFx'])[0,1]
    ncoef_Fy = np.corrcoef(Flapper['Fy'],Pcode['PFy'])[0,1]
    ncoef_Tz = np.corrcoef(Flapper['Tz'],Pcode['PTz'])[0,1]
    
    #Calculate the standard error of the correlation coeff at 0 phase shift
    ncoef_stderr_Fx = np.sqrt((1.-ncoef_Fx**2)/(len(Flapper['Fx'])-2.))
    ncoef_stderr_Fy = np.sqrt((1.-ncoef_Fy**2)/(len(Flapper['Fx'])-2.))
    ncoef_stderr_Tz = np.sqrt((1.-ncoef_Tz**2)/(len(Flapper['Fx'])-2.))
    
    
    
    #Calculate 95% confidence limits for corrcoef (based on Zar)
    
    #Step 1: Transform r to z.
    z_Fx = 0.5*np.log((1.+ncoef_Fx)/(1.-ncoef_Fx))
    z_Fy = 0.5*np.log((1.+ncoef_Fy)/(1.-ncoef_Fy))
    z_Tz = 0.5*np.log((1.+ncoef_Tz)/(1.-ncoef_Tz))
    
    #Step 2: Calculate sigma_z = sqrt(1/(n-3))
    sigma = np.sqrt(1./(len(Flapper['Fx'])-3.))
    
    #Step 3: Calculate lower and upper limits in terms of z.  Note 1.96 = Student's t value for alpha=0.05, 2-tailed, infinite samples
    low_Fx_z = z_Fx - 1.96*sigma
    low_Fy_z = z_Fy - 1.96*sigma
    low_Tz_z = z_Tz - 1.96*sigma
    high_Fx_z = z_Fx + 1.96*sigma
    high_Fy_z = z_Fy + 1.96*sigma
    high_Tz_z = z_Tz + 1.96*sigma
    
    #Step 4: Convert these back from z to r.
    low_Fx = (np.exp(2*low_Fx_z)-1.)/(np.exp(2*low_Fx_z)+1.)
    low_Fy = (np.exp(2*low_Fy_z)-1.)/(np.exp(2*low_Fy_z)+1.)
    low_Tz = (np.exp(2*low_Tz_z)-1.)/(np.exp(2*low_Tz_z)+1.)
    high_Fx = (np.exp(2*high_Fx_z)-1.)/(np.exp(2*high_Fx_z)+1.)
    high_Fy = (np.exp(2*high_Fy_z)-1.)/(np.exp(2*high_Fy_z)+1.)
    high_Tz = (np.exp(2*high_Tz_z)-1.)/(np.exp(2*high_Tz_z)+1.)
    
    
    
    #calculate RMSE
    RMSE_Fx= rmse(Flapper['Fx'],Pcode['PFx'])
    RMSE_Fy= rmse(Flapper['Fy'],Pcode['PFy'])
    RMSE_Tz= rmse(Flapper['Tz'],Pcode['PTz'])
    
    #Calculate RMSE percentage
    RMSE_Fx_per = RMSE_Fx/(max(Pcode['PFx'])-min(Pcode['PFx']))*100
    RMSE_Fy_per = RMSE_Fy/(max(Pcode['PFy'])-min(Pcode['PFy']))*100
    RMSE_Tz_per = RMSE_Tz/(max(Pcode['PTz'])-min(Pcode['PTz']))*100
        
    #print RMSE
    #print RMSE_per
    
    return lags[lag_Fx], lags[lag_Fy], lags[lag_Tz], coef_Fx, coef_Fy, coef_Tz, ncoef_Fx, ncoef_Fy, ncoef_Tz, ncoef_stderr_Fx, ncoef_stderr_Fy, ncoef_stderr_Tz, low_Fx, low_Fy, low_Tz, high_Fx, high_Fy, high_Tz, RMSE_Fx,RMSE_Fy,RMSE_Tz, RMSE_Fx_per,RMSE_Fy_per,RMSE_Tz_per





test_types = ['3D_2cm','3Dmid','3Dedge','dynamic','static','tail_gap','tail_mid']

vids = []
lags_Fx = []
lags_Fy = []
lags_Tz = []
lag_Fx_dT = []
lag_Fy_dT = []
lag_Tz_dT = []
lag_Fx_per = []
lag_Fy_per = []
lag_Tz_per = []
co_Fx = []
co_Fy = []
co_Tz = []
nc_Fx = []
nc_Fy = []
nc_Tz = []
se_Fx = []
se_Fy = []
se_Tz = []
lo_Fx = []
lo_Fy = []
lo_Tz = []
hi_Fx = []
hi_Fy = []
hi_Tz = []
R_Fx = []
R_Fy = []
R_Tz = []
Rp_Fx = []
Rp_Fy = []
Rp_Tz = []

for t in test_types:
    file_pairs = get_paired_files(t)
    
    for key in file_pairs.keys():
        lag_Fx, lag_Fy, lag_Tz, coef_Fx, coef_Fy, coef_Tz, ncoef_Fx, ncoef_Fy, ncoef_Tz, stderr_Fx, stderr_Fy, stderr_Tz, low_Fx, low_Fy, low_Tz, high_Fx, high_Fy, high_Tz, RMSE_Fx,RMSE_Fy,RMSE_Tz, RMSE_Fx_per,RMSE_Fy_per,RMSE_Tz_per = compare(file_pairs[key][0],file_pairs[key][1])
        
        vids.append(key)
        lags_Fx.append(lag_Fx)
        lags_Fy.append(lag_Fy)
        lags_Tz.append(lag_Tz)
        lag_Fx_dT.append(lag_Fx*0.01)
        lag_Fy_dT.append(lag_Fy*0.01)
        lag_Tz_dT.append(lag_Tz*0.01)
        co_Fx.append(coef_Fx)
        co_Fy.append(coef_Fy)
        co_Tz.append(coef_Tz)
        nc_Fx.append(ncoef_Fx)
        nc_Fy.append(ncoef_Fy)
        nc_Tz.append(ncoef_Tz)
        se_Fx.append(stderr_Fx)
        se_Fy.append(stderr_Fy)
        se_Tz.append(stderr_Tz)
        lo_Fx.append(low_Fx)
        lo_Fy.append(low_Fy)
        lo_Tz.append(low_Tz)
        hi_Fx.append(high_Fx)
        hi_Fy.append(high_Fy)
        hi_Tz.append(high_Tz)
        R_Fx.append(RMSE_Fx)
        R_Fy.append(RMSE_Fy)
        R_Tz.append(RMSE_Tz)
        Rp_Fx.append(RMSE_Fx_per)
        Rp_Fy.append(RMSE_Fy_per)
        Rp_Tz.append(RMSE_Tz_per)

        if t != 'static':
            print key
            freq = float(key.split('_')[3][0:3])
            print freq
            p = 1./freq
            
            lag_Fx_per.append(lag_Fx/p)
            lag_Fy_per.append(lag_Fy/p)
            lag_Tz_per.append(lag_Tz/p)
        
        else:
            lag_Fx_per.append('N/A')
            lag_Fy_per.append('N/A')
            lag_Tz_per.append('N/A')
            
            

result = pd.DataFrame(vids,columns=['parentVid'])
result['lags_Fx']=lags_Fx
result['lags_Fy']=lags_Fy
result['lags_Tz']=lags_Tz
result['lag_Fx_dT']=lag_Fx_dT
result['lag_Fy_dT']=lag_Fy_dT
result['lag_Tz_dT']=lag_Tz_dT
result['lag_Fx_per']=lag_Fx_per
result['lag_Fy_per']=lag_Fy_per
result['lag_Tz_per']=lag_Tz_per
result['co_max_Fx']=co_Fx
result['co_max_Fy']=co_Fy
result['co_max_Tz']=co_Tz
result['nc_0shift_Fx']=nc_Fx
result['nc_0shift_Fy']=nc_Fy
result['nc_0shift_Tz']=nc_Tz
result['nc_se_Fx']=se_Fx
result['nc_se_Fy']=se_Fy
result['nc_se_Tz']=se_Tz
result['low_lim_nc_Fx']=lo_Fx
result['low_lim_nc_Fy']=lo_Fy
result['low_lim_nc_Tz']=lo_Tz
result['high_lim_nc_Fx']=hi_Fx
result['high_lim_nc_Fy']=hi_Fy
result['high_lim_nc_Tz']=hi_Tz
result['RMSE_Fx']=R_Fx
result['RMSE_Fy']=R_Fy
result['RMSE_Tz']=R_Tz
result['RMSE_Fx_per']=Rp_Fx
result['RMSE_Fy_per']=Rp_Fy
result['RMSE_Tz_per']=Rp_Tz

result.to_excel('G:/Nonunif_performance_PIV/Pressure code analysis/comparison_Flapper_Pcode_v3_fixed.xlsx')

