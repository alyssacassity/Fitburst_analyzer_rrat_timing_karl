# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:34:19 2024

@author: ktsan
"""
import glob
import os
import Fitburst_singlecut_mod as fbsc
import argparse
import numpy as np
import json
# Load files for analysis (Read in names of candidates)
#files = glob.glob(r'\\wsl.localhost\Ubuntu\home\ktsang45\*.fil')

#for file in files:
    
"""Use Argparse to enable command line inputs"""
parser = argparse.ArgumentParser()
parser.add_argument(
    "pulse_path", 
    action="store", 
    type=str,
    help="Path to .csv file containing all the candidate pulses to be analyzed")
parser.add_argument(
    'fils_path', 
    action='store', 
    type=str,
    help='Path to all .fil files to be analyzed')

args = parser.parse_args()
pulse_path = args.pulse_path
fils_path = args.fils_path

"""Decide on where to cut the blocks here, then pass the cutting parameters
as arguments into Fitburst_singlecut_mod.py to cut each file and generate
an .npz file for each. """
    
#files = os.listdir(r'\\wsl.localhost\Ubuntu\home\ktsang45\positive_bursts_1')
datafile = pulse_path
data = np.genfromtxt(datafile, delimiter=",", dtype=str)
files = data[0]
tstart_list = []
filtime = []
fildm = []

for file in files:
    file = file[file.find('cand'):]
#filfiles = glob.glob(r'\\wsl.localhost\Ubuntu\home\ktsang45\*.fil')
filfiles = glob.glob(fils_path + r'/*.fil')
filparts = file.split('_')
tstart_list.append(filparts[2])
filtime.append(filparts[4])
fildm.append(filparts[6])
    
filmjd = str(int(float(filparts[2])))

fils_to_run = []
for file in filfiles:
    if filmjd in file:
        fils_to_run.append(file)
print(fils_to_run)

toa_list = []
print('test')
#for file_run in fils_to_run:
for i in range(len(files)):
    tstart_list.append(fbsc.singlecut(fils_to_run[0], float(filtime[i])-0.5, float(fildm[i]), filtime[i], float(tstart_list[i])))