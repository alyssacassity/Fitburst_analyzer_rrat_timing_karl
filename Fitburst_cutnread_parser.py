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
    
"""Use Argparse to enable command line inputs"""
parser = argparse.ArgumentParser()
parser.add_argument(
    "pulse_folder", 
    action="store", 
    type=str,
    help="Folder containing all the candidate pulses to be analyzed")
parser.add_argument(
    'fils_path', 
    action='store', 
    type=str,
    help='Path to all .fil files to be analyzed')
parser.add_argument(
    'npz_path', 
    action='store', 
    type=str,
    help='Path to all .npz files to be analyzed')

args = parser.parse_args()
pulse_folder = args.pulse_folder
fils_path = args.fils_path
npz_path = args.npz_path

"""Decide on where to cut the blocks here, then pass the cutting parameters
as arguments into Fitburst_singlecut_mod.py to cut each file and generate
an .npz file for each. """
    
#files = os.listdir(r'\\wsl.localhost\Ubuntu\home\ktsang45\positive_bursts_1')
files = os.listdir(pulse_folder)
filtime = []
fildm = []

for file in files:
    filparts = file.split('_')
    filtime.append(filparts[4])
    fildm.append(filparts[6])
    
filmjd = str(int(float(filparts[2])))

#filfiles = glob.glob(r'\\wsl.localhost\Ubuntu\home\ktsang45\*.fil')
filfiles = glob.glob(fils_path + r'/*.fil')
fils_to_run = []
for file in filfiles:
    if filmjd in file:
        fils_to_run.append(file)
print(fils_to_run)
    
"""Some code here for calling  fitburst_pipeline.py on the .npz files and 
    iterate over them."""

#npz_files = os.listdir(r'C:\Users\ktsan\Desktop\Research\NPZ_files')
npz_files = os.listdir(npz_path)
for file_run in fils_to_run:
    for i in range(len(files)):
        fbsc.singlecut(file_run, float(filtime[i])-0.5, float(fildm[i]), filtime[i])
                
for file in npz_files:
    print(file)
    if '.npz' in file:
    #os.system('python fitburst_pipeline.py ' + r'\\wsl.localhost\Ubuntu\home\ktsang45\NPZ_files' + r'\\'+ file + ' --outfile')
        os.system('python fitburst_pipeline.py '  +' --outfile '+ file )