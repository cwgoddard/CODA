"""lib_coincidence.py -- 
    utility functions for calculating coincidence/acoincidence measurements"""

__author__ = "Chase Goddard"
__email__ = "cwg45@cornell.edu"

import numpy as np
import sys
import random
import collections

sys.path.append('./..')
import read_data

class Coincidence:
    """functions for coincidence/acoincidence calculation"""
    def __init__(self, directory, file_prefix, file_ext, time_window, time_offset, fluor_channel,
            scattered_channel, ROI_begin, ROI_end, num_bins, evt_offset):
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
        self.evt_offset = evt_offset

    def process_fluor_first(self, f_num):
        """Process one netcdf file: Determine number of coincident 
        and accidental ("out of coincidence") pairs. 
        Coincident pairs assumed to begin with an event on the fluorescence detector.
        f_num: the number of the netcdf file to be processed"""
        
        data = read_data.read_netcdf(self.directory, self.file_prefix + str(f_num) + self.file_ext)
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
        
        #data = read_data.read_netcdf(self.directory, self.file_prefix + str(f_num) + self.file_ext)
        data = read_data.read_csv('/home/cwg45/simulate/dump/short.csv')
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
                    if data['time'][k] > t_1 and data['channel'][k] == self.fluor_channel \
                            and self.ROI_begin <= data['E'][k] <= self.ROI_end:
                        res['out_coincidence'][i] = True
                        break
                    k = k + 1  

        #make histograms
        in_counts = np.histogram(
                res['E'][res['in_coincidence'] == True],self.num_bins)[0]
        out_counts = np.histogram(
                res['E'][res['out_coincidence'] == True],self.num_bins)[0]
        print(np.sum(in_counts) , np.sum(out_counts))

        return (in_counts, out_counts)

    def evt_scattered_first(self, f_num, offset_channel):
        """Process one netcdf file: Determine number of coincident 
        and accidental ("out of coincidence") pairs. 
        Coincident pairs assumed to begin with an event on the scattered detector.
        Offset by number of events instead of time
        f_num: the number of the netcdf file to be processed
        offset_channel: channel to use for evt offsets"""
        
        #data = read_data.read_netcdf(self.directory, self.file_prefix + str(f_num) + self.file_ext)
        data = read_data.read_csv('/home/cwg45/simulate/dump/short.csv')

        size = data['E'].size
        
        addl = np.array(list(zip(data['time'] + data['time'][-1], data['E'], data['channel'])),
            dtype = data.dtype)       
        data = np.append(data , addl)

        
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
                while data['time'][k] <= t_0 + self.time_window:
                    #fluorescence photon must be on fluorescence photon detector
                    if data['channel'][k] == self.fluor_channel \
                            and self.ROI_begin <= data['E'][k] <= self.ROI_end: 
                        res['in_coincidence'][i] = True #record coincidence
                        break
                    k = k + 1
                    
        #determine accidentals           
        dups = np.zeros(size)
        for i in range(0,size):
            if data['channel'][i] == self.scattered_channel:
                t_0 = data['time'][i]
                k = i

                ct = 0
                while ct < self.evt_offset: 
                    if data['channel'][k] == offset_channel:
                        ct += 1
                        if ct == self.evt_offset:
                            break
                    k += 1

                assert data['channel'][k] == offset_channel
                
                dither = random.uniform(0,self.time_window/2)
                t_1 = data['time'][k] + dither
                k += 1

                o = k

                while data['time'][k] <= t_1 + self.time_window:
                     #fluorescence photon must be on fluorescence photon detector
                     if data['channel'][k] == self.fluor_channel \
                             and self.ROI_begin <= data['E'][k] <= self.ROI_end: 
                         res['out_coincidence'][i] = True #record coincidence
                         break
                     k = k + 1 

        #make histograms
        in_counts = np.histogram(
                res['E'][res['in_coincidence'] == True],self.num_bins)[0]
        out_counts = np.histogram(
                res['E'][res['out_coincidence'] == True],self.num_bins)[0]

        print(np.sum(in_counts) , np.sum(out_counts))
            
        return (in_counts, out_counts)

    def bin_then_and(self, f_num, bin_width):

        #binning parameter
        div = int(bin_width/20)

        #read in data
        data = read_data.read_netcdf_raw(self.directory,
                self.file_prefix + str(f_num) + self.file_ext)

        #get fluorescence channel events
        fluor = np.extract(data['channel'] == self.fluor_channel, data)

        #get scattered channel events
        scattered = np.extract(data['channel'] == self.scattered_channel, data)

        #get events in fluorescence region
        cond = np.logical_and(fluor['E'] > self.ROI_begin, fluor['E'] < self.ROI_end)
        fluor = np.extract(cond, fluor)

        #store fluorescence events in hash tables
        #really used as hash sets: O(1) add, lookup & remove, associated data is ignored
        coin_ht = {}
        acoin_ht = {}
        for f in fluor:
            coin_ht[int(f[0]/div)] = 1 #1 is just a placeholder value
            acoin_ht[int(f[0]/div)] = 1

        offset = div*2560 #synchrotron period = 2560ns

        #store coincidence/acoincidence energies
        coin = collections.deque() 
        acoin = collections.deque() 

        #look for coincident/acoincident pairs
        for s in scattered:
            if int(s[0]/div) in coin_ht: 
                coin.append(s[1]) #record coincidence
                coin_ht.pop(int(s[0]/div)) #remove fluor event -- no duplicates
            if int((s[0] + offset)/div) in acoin_ht:
                acoin.append(s[1]) #record acoincidence 
                acoin_ht.pop(int((s[0] + offset)/div))

        #make histograms
        in_counts = np.histogram(np.array(list(coin)),self.num_bins, (0,2048))[0]
        out_counts = np.histogram(np.array(list(acoin)),self.num_bins, (0,2048))[0]

        return (in_counts, out_counts)


