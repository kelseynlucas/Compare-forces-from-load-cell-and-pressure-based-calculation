# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 12:31:21 2017

@author: Kelsey N Lucas

Perform an FFT analysis on queen2 forces and torques calculated using 0.002s, 0.004s, and 0.01s time steps.


"""

import pandas as pd
import numpy as np

#Normalization factor for force
rhoSU = 1000*(0.08*0.18)*(0.3**2)

#set file name and paths for reading in data
#dT = 10ms - data were assembled into 1 file previously
dT10ms_file = 'F:/Nonunif_performance_PIV/Pressure code analysis/dynamic/3_3_0angle_2.0Hz_1.5cm_0.3ms/results_3cyc/queen2_results_fulltimetrace_axisfix.xlsx'
#Read in data
dT10ms = pd.read_excel(dT10ms_file, sheetname = 'Sheet1')


#dT = 4 or 2 ms - need to assemble data to one dataframe
folder = 'F:/Nonunif_performance_PIV/Pressure code analysis/smoothing_param_test/dynamic_0angle_2.0Hz'

#dT=4ms
PFx_file = folder + '/queen2_param0.05_rodsim_dT04ms_xforce.xlsx'
PFy_file = folder + '/queen2_param0.05_rodsim_dT04ms_yforce.xlsx'
PTz_file = folder + '/queen2_param0.05_rodsim_dT04ms_torque.xlsx'
#Read in Fx data
PFx = pd.read_excel(PFx_file, sheetname = 'Sheet1', header = None)
#Data is natively in one row.  Need it in one column, so transpose.
PFx = PFx.T
#Label the column something informative so it can be retrieved later
PFx.columns = ['PFx']
#Scale the pressure code data from N/m to m by multiplying by foil depth.
#Sign is assigned to match the definition of +/- given by the Flapper.
PFx['PFx'] = PFx['PFx']*-0.08
#repeat above for Fy
PFy = pd.read_excel(PFy_file, sheetname = 'Sheet1', header = None)
PFy = PFy.T
PFy.columns = ['PFy']
PFy['PFy'] = PFy['PFy']*0.08
#repeat above for Tz
PTz = pd.read_excel(PTz_file, sheetname = 'Sheet1', header = None)
PTz = PTz.T
PTz.columns = ['PTz']
PTz['PTz'] = PTz['PTz']*-80.
#Create a time sequence
time = range(0, len(PFx['PFx']))
#scale the time to ms
time = [t/100. for t in time]
#Store time sequence in final dataframe
dT04ms = pd.DataFrame(time, columns=['time'])
#Add Fx, Fy, and Tz to the final dataframe
dT04ms['PFx'] = PFx['PFx']
dT04ms['PFy'] = PFy['PFy']
dT04ms['PTz'] = PTz['PTz']


#dT=2ms
PFx_file = folder + '/queen2_param0.05_rodsim_dT02ms_xforce.xlsx'
PFy_file = folder + '/queen2_param0.05_rodsim_dT02ms_yforce.xlsx'
PTz_file = folder + '/queen2_param0.05_rodsim_dT02ms_torque.xlsx'
#Read in Fx data
PFx = pd.read_excel(PFx_file, sheetname = 'Sheet1', header = None)
#Data is natively in one row.  Need it in one column, so transpose.
PFx = PFx.T
#Label the column something informative so it can be retrieved later
PFx.columns = ['PFx']
#Scale the pressure code data from N/m to m by multiplying by foil depth.
#Sign is assigned to match the definition of +/- given by the Flapper.
PFx['PFx'] = PFx['PFx']*-0.08
#repeat above for Fy
PFy = pd.read_excel(PFy_file, sheetname = 'Sheet1', header = None)
PFy = PFy.T
PFy.columns = ['PFy']
PFy['PFy'] = PFy['PFy']*0.08
#repeat above for Tz
PTz = pd.read_excel(PTz_file, sheetname = 'Sheet1', header = None)
PTz = PTz.T
PTz.columns = ['PTz']
PTz['PTz'] = PTz['PTz']*-80.
#Create a time sequence
time = range(0, len(PFx['PFx']))
#scale the time to ms
time = [t/100. for t in time]
#Store time sequence in final dataframe
dT02ms = pd.DataFrame(time, columns=['time'])
#Add Fx, Fy, and Tz to the final dataframe
dT02ms['PFx'] = PFx['PFx']
dT02ms['PFy'] = PFy['PFy']
dT02ms['PTz'] = PTz['PTz']


#Normalize everything
dT10ms['PFx'] = dT10ms['PFx']/rhoSU
dT04ms['PFx'] = dT04ms['PFx']/rhoSU
dT02ms['PFx'] = dT02ms['PFx']/rhoSU

dT10ms['PFy'] = dT10ms['PFy']/rhoSU
dT04ms['PFy'] = dT04ms['PFy']/rhoSU
dT02ms['PFy'] = dT02ms['PFy']/rhoSU

dT10ms['PTz'] = dT10ms['PTz']/(rhoSU*0.18*1000)
dT04ms['PTz'] = dT04ms['PTz']/(rhoSU*0.18*1000)
dT02ms['PTz'] = dT02ms['PTz']/(rhoSU*0.18*1000)


#assign sample length
n10ms = len(dT10ms['PFy'])
n04ms = len(dT04ms['PFy'])
n02ms = len(dT02ms['PFy'])

#assign sampling frequency
fs10ms = 100.
fs04ms = 250.
fs02ms = 500.

#do fft
yFx = np.fft.fft(dT10ms['PFx'])
yFy = np.fft.fft(dT10ms['PFy'])
yTz = np.fft.fft(dT10ms['PTz'])
y10ms = pd.DataFrame(yFx.real, columns=['yFx'])
y10ms['yFy'] = yFy.real
y10ms['yTz'] = yTz.real
y10ms['freq'] = np.linspace(0,len(yFx)-1,len(yFx))*(fs10ms/n10ms)
y10ms['powerFx']=(abs(y10ms['yFx'])**2.)/n10ms
y10ms['powerFy']=(abs(y10ms['yFy'])**2.)/n10ms
y10ms['powerTz']=(abs(y10ms['yTz'])**2.)/n10ms

yFx = np.fft.fft(dT04ms['PFx'])
yFy = np.fft.fft(dT04ms['PFy'])
yTz = np.fft.fft(dT04ms['PTz'])
y04ms = pd.DataFrame(yFx.real, columns=['yFx'])
y04ms['yFy'] = yFy.real
y04ms['yTz'] = yTz.real
y04ms['freq'] = np.linspace(0,len(yFx)-1,len(yFx))*(fs04ms/n04ms)
y04ms['powerFx']=(abs(y04ms['yFx'])**2.)/n04ms
y04ms['powerFy']=(abs(y04ms['yFy'])**2.)/n04ms
y04ms['powerTz']=(abs(y04ms['yTz'])**2.)/n04ms


yFx = np.fft.fft(dT02ms['PFx'])
yFy = np.fft.fft(dT02ms['PFy'])
yTz = np.fft.fft(dT02ms['PTz'])
y02ms = pd.DataFrame(yFx.real, columns=['yFx'])
y02ms['yFy'] = yFy.real
y02ms['yTz'] = yTz.real
y02ms['freq'] = np.linspace(0,len(yFx)-1,len(yFx))*(fs02ms/n02ms)
y02ms['powerFx']=(abs(y02ms['yFx'])**2.)/n02ms
y02ms['powerFy']=(abs(y02ms['yFy'])**2.)/n02ms
y02ms['powerTz']=(abs(y02ms['yTz'])**2.)/n02ms


y10msSave = y10ms[0:int(n10ms/2.)]
y04msSave = y04ms[0:int(n04ms/2.)]
y02msSave = y02ms[0:int(n02ms/2.)]

y10msSave.to_excel('F:/Nonunif_performance_PIV/Pressure code analysis/smoothing_param_test/dynamic_0angle_2.0Hz/fft_dT10ms.xlsx')
y04msSave.to_excel('F:/Nonunif_performance_PIV/Pressure code analysis/smoothing_param_test/dynamic_0angle_2.0Hz/fft_dT04ms.xlsx')
y02msSave.to_excel('F:/Nonunif_performance_PIV/Pressure code analysis/smoothing_param_test/dynamic_0angle_2.0Hz/fft_dT02ms.xlsx')

