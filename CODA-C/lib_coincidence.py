"""lib_coincidence.py -- 
    utility functions for calculating coincidence/acoincidence measurements"""

__author__ = "Chase Goddard"
__email__ = "cwg45@cornell.edu"

import numpy as np
import sys
import random
import collections

sys.path.append('/home/cwg45/CODA/')
import read_data

class Coincidence:
    """functions for coincidence/acoincidence calculation"""
    def __init__(self, directory, file_prefix, file_ext, time_window, 
            time_offset, fluor_channel, scattered_channel, ROI_begin, ROI_end,
            num_bins, evt_offset):
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
        Coincident pairs assumed to begin with an event on fluor detector.
        f_num: the number of the netcdf file to be processed"""
        
        data = read_data.read_netcdf(self.directory,
                self.file_prefix + str(f_num) + self.file_ext)
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
         
        in_counts = np.histogram(
                res['E'][res['in_coincidence'] == True],self.num_bins)[0]
           
        out_c = np.zeros(self.num_bins)
        #determine accidentals           
        for j in range(1,7): 
            res = np.array(list(zip(data['E'], np.zeros(size), np.zeros(size))),
            dtype = {'names':['E', 'in_coincidence', 'out_coincidence'],
                     'formats':['int32', bool, bool]})

            for i in range(0,size):
                if data['channel'][i] == self.fluor_channel \
                       and self.ROI_begin <= data['E'][i] <= self.ROI_end:
                    t_0 = data['time'][i]
                    t_1 = t_0 + self.time_offset*j 
                    k = i
                    while k < size and data['time'][k] <= t_1 + self.time_window:
                        if data['time'][k] >= t_1 and \
                                data['channel'][k] == self.scattered_channel:
                            res['out_coincidence'][k] = True
                            break
                        k = k + 1
            out_c +=  np.histogram(
                res['E'][res['out_coincidence'] == True],self.num_bins)[0]
 
        #make histograms
        out_counts = np.zeros(self.num_bins)
           
        return (in_counts, out_counts)
     

    def process_scattered_first(self, f_num):
        """Process one netcdf file: Determine number of coincident 
        and accidental ("out of coincidence") pairs. 
        Coincident pairs assumed to begin with an event on scattered detector.
        f_num: the number of the netcdf file to be processed"""
        
        #data = read_data.read_netcdf(self.directory, 
            #self.file_prefix + str(f_num) + self.file_ext)
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
                    if data['time'][k] > t_1 \
                            and data['channel'][k] == self.fluor_channel \
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
        Coincident pairs assumed to begin with an event on scattered detector.
        Offset by number of events instead of time
        f_num: the number of the netcdf file to be processed
        offset_channel: channel to use for evt offsets"""
        
        #data = read_data.read_netcdf(self.directory,
            #self.file_prefix + str(f_num) + self.file_ext)
        data = read_data.read_csv('/home/cwg45/simulate/dump/short.csv')

        size = data['E'].size
        
        addl = np.array(list(zip(data['time'] + data['time'][-1],
            data['E'], data['channel'])), dtype = data.dtype)       
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
                     #fluorescence photon must be on fluorescence detector
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

        print(f_num)

        #binning parameter
        div = int(bin_width/20)
	
        #read in data
        try:
            data = read_data.read_netcdf_raw(self.directory,
                    self.file_prefix + str(f_num) + self.file_ext)
        except:
            return ((np.zeros(self.num_bins, dtype=np.int32), np.zeros(self.num_bins, 
                dtype = np.int32)), 0)

        #get fluorescence channel events
        fluor = np.extract(data['channel'] == self.fluor_channel, data)

        #get scattered channel events
        scattered = np.extract(data['channel'] == self.scattered_channel, data)

        #get events in fluorescence region
        cond = np.logical_and(fluor['E'] > self.ROI_begin,
                fluor['E'] < self.ROI_end)
        fluor = np.extract(cond, fluor)

        #store fluorescence events in hash tables
        #really used as hash sets: O(1) add, lookup & remove,
        #associated data is ignored
        coin_ht = {}
        acoin_ht = {}
        for f in fluor:
            coin_ht[int(f[0]/div)] = 1 #1 is just a placeholder value
            acoin_ht[int(f[0]/div)] = 1 #all that matters is the key

        offset = div*95. #synchrotron period = 2560ns

        #store coincidence/acoincidence energies
        coin = collections.deque() 
        acoin = collections.deque() 
        listofaccidentals = []	

        #look for coincident/acoincident pairs
        for i in range(1,11):
            #acoin_ht = {}
            #for f in fluor:
            #    acoin_ht[int(f[0]/div)] = 1 #all that matters is the key
            for s in scattered:
                if i == 1:
                    if int(s[0]/div) in coin_ht: 
                        coin.append(s[1]) #record coincidence
                        coin_ht.pop(int(s[0]/div)) #remove fluor event -- no duplicates
                if int((s[0] + offset*i)/div) in acoin_ht:
                    acoin.append(s[1]) #record acoincidence 
                    acoin_ht.pop(int((s[0] + offset*i) / div))
            temp = np.histogram(np.array(list(acoin)),self.num_bins,(0,2048))[0]
            listofaccidentals.append(sum(temp)-sum(listofaccidentals))

        #make histograms
        in_counts = np.histogram(np.array(list(coin)),self.num_bins, 
                (0,2048))[0]
        #in_counts = np.empty(self.num_bins)
        #in_counts.fill(np.std(listofaccidentals))
        #listofaccidentals = np.array(listofaccidentals)
        #in_counts = np.append([np.std(listofaccidentals)],np.zeros(99))
        out_counts = np.histogram(np.array(list(acoin)),self.num_bins, 
                (0,2048))[0]

        return ((in_counts, out_counts), fluor.size)


    def inner_bootstrap(self, data, bin_width):

        #binning parameter
        div = int(bin_width/20)

        #get fluorescence channel events
        fluor = np.extract(data['channel'] == self.fluor_channel, data)

        #get scattered channel events
        scattered = np.extract(data['channel'] == self.scattered_channel, data)

        #get events in fluorescence region
        cond = np.logical_and(fluor['E'] > self.ROI_begin,
                fluor['E'] < self.ROI_end)
        fluor = np.extract(cond, fluor)
        
        #get a random sample 
        fluor = np.random.choice(fluor, fluor.size, replace=True)


        #store fluorescence events in hash tables
        #increment value by 1 for each event (bootstrapping allows duplicates)
        coin_ht = {}
        acoin_ht = {}
        for f in fluor:
            #if the key is already associated with a value, just increment
            if int(f[0]/div) in coin_ht:
                coin_ht[int(f[0]/div)] += 1
                acoin_ht[int(f[0]/div)] += 1
            else: #otherwise, initialize the key to 1
                coin_ht[int(f[0]/div)] = 1 
                acoin_ht[int(f[0]/div)] = 1


        offset = div*128 #synchrotron period = 2560ns/20ns*bin^-1 = 128 bins

        #store coincidence/acoincidence energies
        coin = collections.deque() 
        acoin = collections.deque() 

        #look for coincident/acoincident pairs
        for s in scattered:
            if int(s[0]/div) in coin_ht: 
                coin.append(s[1]) #record coincidence
                coin_ht[int(s[0]/div)] -= 1 #subtract one fluor event
                if coin_ht[int(s[0]/div)] <= 0:
                    coin_ht.pop(int(s[0]/div)) #remove fluor event if there are none left
            if int((s[0] + offset)/div) in acoin_ht:
                acoin.append(s[1]) #record acoincidence 
                acoin_ht[int((s[0] + offset)/div)] -= 1
                if acoin_ht[int((s[0]+offset)/div)] <= 0:
                    acoin_ht.pop(int((s[0] + offset) / div))

        #make histograms
        in_counts = np.histogram(np.array(list(coin)),self.num_bins, 
                (0,2048))[0]
        out_counts = np.histogram(np.array(list(acoin)),self.num_bins, 
                (0,2048))[0]

        return ((in_counts, out_counts), fluor.size)


    def bootstrap(self, f_num, bin_width, num_samples):
        #read in data
        print(f_num)

        #read in data, or return 0's if the file is broken
        try:
            data = read_data.read_netcdf_raw(self.directory,
                    self.file_prefix + str(f_num) + self.file_ext)
        except:
            return ((np.zeros(self.num_bins, dtype=np.int32), np.zeros(self.num_bins, 
                dtype = np.int32)), 0)

        #store bootstrap results
        tbl = collections.deque()
        for b_id in range(0,num_samples):
            #perform one measurement
            res = self.inner_bootstrap(data, bin_width)
            #store results
            tbl.append(res)

        return tbl
