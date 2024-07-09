# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:07:48 2024

@author: ktsan
"""
import os
import argparse
import numpy as np
import json
from multiprocessing import Process

"""Use Argparse to enable command line inputs"""
parser = argparse.ArgumentParser()
parser.add_argument(
    'json_path',
    action='store', 
    type=str,
    help='Path to all .npz files to be analyzed')
args = parser.parse_args()
npz_path = args.npz_path
filtime = []
procs = []
npz_files = [i for i in os.listdir(npz_path) if '.npz' in i]

def fitpipe(i):
    os.system('python fitburst_pipeline.py '  +' --outfile '+ npz_files[i] )
"""Multiprocessing code"""
if __name__ == '__main__':
    for i in range(len(npz_files)):
        filparts = npz_files[i].split('_')
        filtime.append(filparts[3])
        proc = Process(target=fitpipe, args=(i))
        procs.append(proc)
        proc.start()
        
    for proc in procs:
        proc.join()
toa_list = []
tstart_list = []
results_files = [i for i in os.listdir(npz_path) if '.json' in i]
""" Some code for reading the TOA from the results json file"""
results_toa = []
ref_freqs = []
mjd_errors = []
proc = []
def make_tim(i):
    with open(results_files[i], 'r') as f:
        data = json.load(f)
        results_toa.append((data['model_parameters']['arrival_time'][0]-0.5)/86400)
        ref_freqs.append(800)
        if (isinstance(data['fit_statistics']['bestfit_uncertainties']['arrival_time'][0], float) and 
        (not np.isnan(data['fit_statistics']['bestfit_uncertainties']['arrival_time'][0]))) :
            mjd_errors.append(data['fit_statistics']['bestfit_uncertainties']['arrival_time'][0])
        else:
            mjd_errors.append(1e-6)
if __name__ == '__main__':
    for i in range(len(results_files)):
        proc = Process(target=make_tim, args=(i))
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join()
        
print("Results_TOA", results_toa)
filtime = [float(i) for i in filtime]
filtime = np.array(filtime)/86400
print(len(npz_files))
    
"""Some code here for calling  fitburst_pipeline.py on the .npz files and 
    iterate over them."""

toa_list.append(tstart_list+filtime+results_toa)
print("TOA_list", toa_list)
print(len(toa_list))


""" Save data to .tim file"""
res_file = open('pulsar_timing_results.tim', 'w')
txt_list = []
for i in range(len(npz_files)):
    txt_line = (npz_files[i].removesuffix('_'+ filtime[i]+'.npz')
                + ' ' + str(ref_freqs[i]) + ' ' + str(toa_list[0][i]) + ' ' 
                + str(mjd_errors[i]) + ' y  \n')
    txt_list.append(txt_line)
res_file.writelines(txt_list)
res_file.close()
    