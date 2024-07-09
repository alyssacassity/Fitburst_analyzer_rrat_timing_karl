# -*- coding: utf-8 -*-
"""
Created on Wed May 22 12:37:22 2024

@author: ktsan
"""
import os
import argparse
import numpy as np
import json
import sys # Needed for bash script interfacing

"""Use Argparse to enable command line inputs"""
parser = argparse.ArgumentParser()
parser.add_argument(
    'index', 
    action='store', 
    type=int,
    help='Index')
parser.add_argument(
    'json_path',
    action='store', 
    type=str,
    help='Path to all .npz files to be analyzed')
print('test')
args = parser.parse_args()
ind = args.index
npz_path = args.npz_path
print(npz_path)
filtime = []
tstart_list = []
#npz_files = os.listdir(r'C:\Users\ktsan\Desktop\Research\NPZ_files')
npz_files = [i for i in os.listdir(npz_path) if '.npz' in i]
#for i in range(len(npz_files)):
filparts = npz_files[ind].split('_')
filtime.append(filparts[3])
tstart_list.append(filparts[4])
print(npz_files[ind])

    #os.system('python fitburst_pipeline.py ' + r'\\wsl.localhost\Ubuntu\home\ktsang45\NPZ_files' + r'\\'+ file + ' --outfile')
os.system('python '+ npz_path + 'fitburst_pipeline.py '  +' --outfile '+ npz_path + npz_files[ind] + ' ' + npz_path )
