# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 13:32:48 2015

@author: Kelsey
"""

def downsample_Flapper(FlapperData,PcodeData):
    """
    Reads flapper data and corresponding pressure code data and downsamples the
    flapper data to match the number of points in the pressure code data.
    
    Input:
    -FlapperData - file name for the flapper data file (post processing, rod
        resolving forces, rod subtraction, and filtering)
    -PcodeData - file containing the full time sequence for the queen2 results
    
    """
    #make sure that the dataframe package is loaded
    import pandas as pd
    
    #read in the flapper data
    Fdata = pd.read_excel(FlapperData, sheetname = 'Sheet1', index_col=0)
    
    #read in the pressure code data
    Pdata = pd.read_excel(PcodeData, sheetname = 'Sheet1', index_col=0)

    #initialize an index
    i=0
    
    #initialize storage for the downsampled Fx, Fy, and Tz sequences
    dsFx=[]
    dsFy=[]
    dsTz=[]
    
    #for each row in the pressure data table (ensures that the downsampled Flapper
    #sequence has the same number of samples as the pressure code data)
    for row in range(0,len(Pdata['PFx'])):
       
       #add the next data point, indicated by the index i
       dsFx.append(Fdata['Fx_noRod'][i])
       dsFy.append(Fdata['Fy_noRod'][i])
       dsTz.append(Fdata['Tz_noRod'][i])
       
       #Flapper data was sampled at 1000Hz, but pressure code at 100Hz.  Only 
       #need every 10th point in the Flapper data sequence.
       i+=10
       
    #make a data table.  The first column added is downsampled-Fx.
    df=pd.DataFrame(dsFx, columns=['Fx'])
    
    #add downsampled-Fy and downsampled-Tz to the table
    df['Fy']=dsFy
    df['Tz']=dsTz

    #make a savepath for the downsampled data.  Give it the same name as full
    #sequence, but with '_downsampled' at the end.
    savepath = FlapperData[0:-5] + '_downsampled_resfix.xlsx'

    #save out the downsampled data table.
    df.to_excel(savepath)
    


    

def get_Flapper_files(directory, testtypes):
    """
    Find all the flapper files in the directory for each test type in testtypes.
    Assumes Flapper_analysis_wrapper.py was ran to generate files.
    """
    
    #load file operations
    from os import listdir
    from os.path import isfile, join
    
    #initialize a list where the flapper data file names will be stored
    finalfilelist = []
    
    #for each test type,
    for test in testtypes:
        
        #create the full path for the test's data folder
        path = directory + '/' + test + '/'
        
        #get the full path for each file (but not subfolder) in the test's data folder
        files = [path + f for f in listdir(path) if isfile(join(path,f))]
        
        #for these paths,
        for f in files:
            
            #if they don't have '_3reps' in the name (added to the name during phase-averaging
            #and net value calculations to indicate that this file contains those
            #results),
            if not '_3reps' in f:
                
                if not '_downsample' in f:
                    
                    if 'resfix' in f:
                
                        #add the file to the final list
                        finalfilelist.append(f)
                

    #return all the files
    return finalfilelist
    
    
    
    

def get_queen2_files(directory, testtypes):
    """
    Find all the queen2 files in the directory for each test type in testtypes.
    Assumes queen2_fulltimetrace.py was ran to generate files.
    """
    
    #load useful file operations
    from os import listdir
    from os.path import isfile, join
    
    #initialize a list for the queen2_fulltimetrace file names & paths
    finalfilelist = []
    
    #for each test,
    for test in testtypes:
        #make the full path for that test's data folder
        path = directory + '/' + test + '/'
        
        #make the paths for the folder where the queen2 data are stored, for
        #each video analyzed
        parentvids = [path + f + '/results_3cyc/' for f in listdir(path) if not isfile(join(path,f))]
        
        #for each parent video's queen2 data folder,
        for f in parentvids:
            #make the full name and path for the queen2 results
            fullname = f + 'queen2_results_fulltimetrace_axisfix.xlsx'
            #add this to the list
            finalfilelist.append(fullname)
    
    #return the list
    return finalfilelist 
    
#ID all the test times
tests = ['dynamic','3Dmid','3D_2cm','3Dedge','static']

#get the flapper data files
FlapperFiles = get_Flapper_files('G:/Nonunif_performance_PIV/Data_fromFlapper_analyzed/final_analysis_7Hzfiltercutoff', tests)
#get the queen2 data files
Queen2Files = get_queen2_files('G:/Nonunif_performance_PIV/Pressure code analysis',tests)


#for each flapper data file,
for f in FlapperFiles:
    #find the name of the parent video for that file
    vid = f.split('/')[-1][0:-12]
    
    #look through all the queen2 data files
    for q in Queen2Files:
        
        #find the one for that parent video
        if vid in q:
            
            #downsample that flapper data file according to the number of 
            #samples in the queen2 data file
            downsample_Flapper(f,q)
            
            #once the queen2 file for the parent video is found, no need to keep
            #looking.  Stop and move on to the next flapper file.
            break
    
    
    
    
    
    
    
    
    
    
    
    