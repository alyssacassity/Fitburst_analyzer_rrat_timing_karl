# -*- coding: utf-8 -*-
"""
Created on Tue May  7 12:45:27 2024

@author: Karl
"""

import numpy as np
from rich.pretty import Pretty
from sigpyproc.readers import FilReader
import matplotlib.pyplot as plt
from iqrm  import iqrm_mask
from scipy.stats import median_abs_deviation as mad

def singlecut(fil_name, t_start, disp_measure, fil_time, isddp=True):
    fbfile = FilReader(fil_name)
    sub = '.fil'
    # Remove suffix from filename for naming .npz file later
    fil_short_name = [i for i in fil_name.split("/") if sub in i][0].removesuffix(sub)
    fbh = fbfile.header

    # Define variables
    t_block = 1
    nsamps = int(t_block/fbh.tsamp)
    downsamp = 12
    dsampfreq = 8
    nsamps = nsamps-nsamps%downsamp
    #t_start = 396.3
    #print(nsamps)
    start_samp = round(t_start/fbh.tsamp)
    #print(start_samp)
    fblock = fbfile.read_block(start_samp, nsamps)
    #fbd = fbfile.read_dedisp_block(round(200/fbh.tsamp), 59000, 112.5)
    #print(fblock)
    

    
    # Identify bad channels for masking and dedisperse data if needed
    
    if isddp:
        fbt = fblock.dedisperse(dm=disp_measure)
    else: 
        fbt = fblock
    mn = np.mean(fbt, axis = 1)
    mask = (mn == 0)
    fbt[mask,:] = np.median(fbt[~mask,:])
    mad_mask = mad(fbt, axis=1) > 1
    mask[mad_mask] = True
    fbt[mask,:] = np.median(fbt[~mask,:])
    std_on = np.std(fbt, axis=1)
    
    iqrmmask,votes = iqrm_mask(std_on, threshold=1)
    mask[iqrmmask] = True
    
    # Remove bad channels using mask
    for i in range(len(mask)):
        if mask[i]:
            fbt[i,:] = np.median(fbt[i,:])
    fbt = fbt.downsample(downsamp, dsampfreq)
    fbt = fbt.normalise()

    
    # Create metadata and burst parameters dictionaries
    metadata = dict(bad_chans = 0, freqs_bin0 = fbh.fch1, is_dedispersed = isddp,
                    num_time = nsamps/downsamp, num_freq = fbh.nchans/dsampfreq,
                    times_bin0 = fbh.tstart+t_start/86400,  res_time = fbh.tsamp*downsamp, res_freq = fbh.foff*dsampfreq)
    burst_parameters = dict(ref_freq = [600.2], amplitude = [np.log10(np.max(fbt))], arrival_time = [0.5],
                            burst_width = [0.02], dm = [disp_measure], dm_index = [-2], 
                            scattering_index = [-4], scattering_timescale = [0.01], spectral_index = [0], 
                            spectral_running = [0])
                            
    
    # Plot data and save
    plt.imshow(fbt, aspect='auto')
    #plt.imshow(fblock, aspect='auto')
    data_full = fbt
    #print(fil_short_name)
    #print(r'~/NPZ_files/'+ fil_short_name + '_' + fil_time + '.npz')
    np.savez(
        (fil_short_name + '_' + fil_time + ".npz"), 
        data_full=data_full, 
        metadata=metadata, 
        burst_parameters=burst_parameters,
        )
        
    return fbh.tstart