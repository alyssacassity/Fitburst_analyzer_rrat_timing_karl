# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:34:19 2024

@author: ktsan
"""
import glob
import os
import Fitburst_singlecut_mod as fbsc
import argparse
# Load files for analysis (Read in names of candidates)
#files = glob.glob(r'\\wsl.localhost\Ubuntu\home\ktsang45\*.fil')

#for file in files:


"""Decide on where to cut the blocks here, then pass the cutting parameters
as arguments into Fitburst_singlecut_mod.py to cut each file and generate
an .npz file for each. """
    
files = os.listdir(r'\\wsl.localhost\Ubuntu\home\ktsang45\positive_bursts_1')
filtime = []
fildm = []

for file in files:
    filparts = file.split('_')
    filtime.append(filparts[4])
    fildm.append(filparts[6])
    
filmjd = str(int(float(filparts[2])))

filfiles = glob.glob(r'\\wsl.localhost\Ubuntu\home\ktsang45\*.fil')
fils_to_run = []
for file in filfiles:
    if filmjd in file:
        fils_to_run.append(file)
        
for file in fils_to_run:
    for i in range(len(files)):
        fbsc.singlecut(file, float(filtime[i])-0.5, float(fildm[i]))
                
    
"""Some code here for calling  fitburst_pipeline.py on the .npz files and 
    iterate over them."""

npz_files = os.listdir(r'C:\Users\ktsan\Desktop\Research\NPZ_files')
for file in npz_files:
    os.system('python fitburst_pipeline.py ' + r'C:\Users\ktsan\Desktop\Research\NPZ_files' + r'\\'+ file)