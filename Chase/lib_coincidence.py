"""lib_coincidence.py -- 
    utility functions for calculating coincidence/acoincidence measurements"""

__author__ = "Chase Goddard"
__email__ = "cwg45@cornell.edu"

import numpy as np
import sys

sys.path.append('./..')
import xMAP

def process_fluor_first(f_num):
    """Process one netcdf file: Determine number of coincident 
    and accidental ("out of coincidence") pairs. 
    Coincident pairs assumed to begin with an event on the fluorescence detector.
    f_num: the number of the netcdf file to be processed"""
    
    data = xMAP.read_netcdf(directory, file_prefix + str(f_num) + file_ext)
    size = data['E'].size
    
    #set up result table
    res = np.array(list(zip(data['E'], np.zeros(size), np.zeros(size))),
        dtype = {'names':['E', 'in_coincidence', 'out_coincidence'],
                 'formats':['int32', bool, bool]})
    
    #determine coincidence
    for i in range(0,size):
        #check for fluorescence
        if data['channel'][i] == fluor_channel \
               and ROI_begin <= data['E'][i] <= ROI_end: 
            #time of fluorescence hit
            t_0 = data['time'][i]
            k = i
            #search for scattered photon within time window
            while k < size and data['time'][k] <= t_0 + time_window:
                #scattered photon must be on scattered photon detector
                if data['channel'][k] == scattered_channel:
                    res['in_coincidence'][k] = True #record coincidence
                    break
                k = k + 1
                
    #determine accidentals           
    for i in range(0,size):
        if data['channel'][i] == fluor_channel \
               and ROI_begin <= data['E'][i] <= ROI_end:
            t_0 = data['time'][i]
            t_1 = t_0 + time_offset #skip ahead time_offset seconds
            k = i
            while k < size and data['time'][k] <= t_1 + time_window:
                if data['time'][k] >= t_1 and \
                        data['channel'][k] == scattered_channel:
                    res['out_coincidence'][k] = True
                    break
                k = k + 1  
    #make histograms
    in_counts = np.histogram(
            res['E'][res['in_coincidence'] == True],num_bins)[0]
    out_counts = np.histogram(
            res['E'][res['out_coincidence'] == True],num_bins)[0]
        
    return (in_counts, out_counts)
 

def process_scattered_first(f_num):
    """Process one netcdf file: Determine number of coincident 
    and accidental ("out of coincidence") pairs. 
    Coincident pairs assumed to begin with an event on the scattered detector.
    f_num: the number of the netcdf file to be processed"""
    
    data = xMAP.read_netcdf(directory, file_prefix + str(f_num) + file_ext)
    size = data['E'].size
    
    #set up result table
    res = np.array(list(zip(data['E'], np.zeros(size), np.zeros(size))),
        dtype = {'names':['E', 'in_coincidence', 'out_coincidence'],
                 'formats':['int32', bool, bool]})
    
    #determine coincidence
    for i in range(0,size):
        #check for scattered photon
        if data['channel'][i] == scattered_channel:
            #time of scattered hit
            t_0 = data['time'][i]
            k = i
            #search for fluorescence photon within time window
            while k < size and data['time'][k] <= t_0 + time_window:
                #fluorescence photon must be on fluorescence photon detector
                if data['channel'][k] == fluor_channel \
                        and ROI_begin <= data['E'][k] <= ROI_end: 
                    res['in_coincidence'][i] = True #record coincidence
                    break
                k = k + 1
                
    #determine accidentals           
    for i in range(0,size):
        if data['channel'][i] == scattered_channel:
            t_0 = data['time'][i]
            t_1 = t_0 + time_offset #skip ahead time_offset seconds
            k = i
            while k < size and data['time'][k] <= t_1 + time_window:
                if data['time'][k] >= t_1 and data['channel'][k] == fluor_channel \
                        and ROI_begin <= data['E'][k] <= ROI_end:
                    res['out_coincidence'][i] = True
                    break
                k = k + 1  

    #make histograms
    in_counts = np.histogram(
            res['E'][res['in_coincidence'] == True],num_bins)[0]
    out_counts = np.histogram(
            res['E'][res['out_coincidence'] == True],num_bins)[0]
        
    return (in_counts, out_counts)
