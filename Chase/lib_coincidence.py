"""lib_coincidence.py -- 
    utility functions for calculating coincidence/acoincidence measurements"""

__author__ = "Chase Goddard"
__email__ = "cwg45@cornell.edu"

import numpy as np
import sys

sys.path.append('./..')
import xMAP

class Coincidence:
    """functions for coincidence/acoincidence calculation"""
    def __init__(self, directory, file_prefix, file_ext, time_window, time_offset, fluor_channel,
            scattered_channel, ROI_begin, ROI_end, num_bins):
        self.directory = directory
        self.file_prefix = file_prefix
        self.file_ext = file_ext
        self.time_window = time_window
        self.time_offset = time_offset
        self.fluor_channel = fluor_channel
        self.scattered_channel = scattered_channel
        self.ROI_begin = ROI_begin
        self.ROI_end = ROI_end
        self.num_bins = num_bins

    def process_fluor_first(self, f_num):
        """Process one netcdf file: Determine number of coincident 
        and accidental ("out of coincidence") pairs. 
        Coincident pairs assumed to begin with an event on the fluorescence detector.
        f_num: the number of the netcdf file to be processed"""
        
        data = xMAP.read_netcdf(self.directory, self.file_prefix + str(f_num) + self.file_ext)
        size = data['E'].size
        
        #set up result table
        res = np.array(list(zip(data['E'], np.zeros(size), np.zeros(size))),
            dtype = {'names':['E', 'in_coincidence', 'out_coincidence'],
                     'formats':['int32', bool, bool]})
        
        #determine coincidence
        for i in range(0,size):
            #check for fluorescence
            if data['channel'][i] == self.fluor_channel \
                   and self.ROI_begin <= data['E'][i] <= self.ROI_end: 
                #time of fluorescence hit
                t_0 = data['time'][i]
                k = i
                #search for scattered photon within time window
                while k < size and data['time'][k] <= t_0 + self.time_window:
                    #scattered photon must be on scattered photon detector
                    if data['channel'][k] == self.scattered_channel:
                        res['in_coincidence'][k] = True #record coincidence
                        break
                    k = k + 1
                    
        #determine accidentals           
        for i in range(0,size):
            if data['channel'][i] == self.fluor_channel \
                   and self.ROI_begin <= data['E'][i] <= self.ROI_end:
                t_0 = data['time'][i]
                t_1 = t_0 + self.time_offset #skip ahead time_offset seconds
                k = i
                while k < size and data['time'][k] <= t_1 + self.time_window:
                    if data['time'][k] >= t_1 and \
                            data['channel'][k] == self.scattered_channel:
                        res['out_coincidence'][k] = True
                        break
                    k = k + 1  
        #make histograms
        in_counts = np.histogram(
                res['E'][res['in_coincidence'] == True],self.num_bins)[0]
        out_counts = np.histogram(
                res['E'][res['out_coincidence'] == True],self.num_bins)[0]
            
        return (in_counts, out_counts)
     

    def process_scattered_first(self, f_num):
        """Process one netcdf file: Determine number of coincident 
        and accidental ("out of coincidence") pairs. 
        Coincident pairs assumed to begin with an event on the scattered detector.
        f_num: the number of the netcdf file to be processed"""
        
        data = xMAP.read_netcdf(self.directory, self.file_prefix + str(f_num) + self.file_ext)
        size = data['E'].size
        
        #set up result table
        res = np.array(list(zip(data['E'], np.zeros(size), np.zeros(size))),
            dtype = {'names':['E', 'in_coincidence', 'out_coincidence'],
                     'formats':['int32', bool, bool]})
        
        #determine coincidence
        for i in range(0,size):
            #check for scattered photon
            if data['channel'][i] == self.scattered_channel:
                #time of scattered hit
                t_0 = data['time'][i]
                k = i
                #search for fluorescence photon within time window
                while k < size and data['time'][k] <= t_0 + self.time_window:
                    #fluorescence photon must be on fluorescence photon detector
                    if data['channel'][k] == self.fluor_channel \
                            and self.ROI_begin <= data['E'][k] <= self.ROI_end: 
                        res['in_coincidence'][i] = True #record coincidence
                        break
                    k = k + 1
                    
        #determine accidentals           
        for i in range(0,size):
            if data['channel'][i] == self.scattered_channel:
                t_0 = data['time'][i]
                t_1 = t_0 + self.time_offset #skip ahead time_offset seconds
                k = i
                while k < size and data['time'][k] <= t_1 + self.time_window:
                    if data['time'][k] >= t_1 and data['channel'][k] == self.fluor_channel \
                            and self.ROI_begin <= data['E'][k] <= self.ROI_end:
                        res['out_coincidence'][i] = True
                        break
                    k = k + 1  

        #make histograms
        in_counts = np.histogram(
                res['E'][res['in_coincidence'] == True],self.num_bins)[0]
        out_counts = np.histogram(
                res['E'][res['out_coincidence'] == True],self.num_bins)[0]
            
        return (in_counts, out_counts)
