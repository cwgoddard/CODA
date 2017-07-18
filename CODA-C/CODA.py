"""OSPS.py: Processing script for IAB data. Reads in netcdf files from xMAP, 
performs coincidence calculation, and sums the result."""

__author__ = "Chase Goddard"
__email__ = "cwg45@cornell.edu"

#imports
import numpy as np
from functools import *
import lib_coincidence

#import joblib -- for parallelization
#joblib isn't on CLASSE servers by default, so some trickery is required
import sys
#modify path to include joblib package
sys.path.append('/home/cwg45/local/lib/python3.3/site-packages') 

from joblib import Parallel, delayed

#to parse command-line arguments:
import argparse
import ast

parser = argparse.ArgumentParser()
parser.add_argument("data_directory", help="xMAP data directory")
parser.add_argument("file_prefix", help="netCDF file prefix, e.g list_")
parser.add_argument("file_start", help="number of first file in directory", type = int)
parser.add_argument("file_end", help="number of last file in directory", type = int)
parser.add_argument("file_ext", help="netCDF file extension, e.g. .nc") 

parser.add_argument("write_directory", help="directory to write output .csv to") 
parser.add_argument("write_name", help="name of output .csv") 

parser.add_argument("time_bin", help="length of time bin (ns)", 
        type = int)
parser.add_argument("fluor_begin", help="beginning of fluorescence region in 10s of eV", 
        type = int)
parser.add_argument("fluor_end", help="end of fluorescence region in 10s of eV", 
        type = int)
#parser.add_argument("time_offset", help="offset to to determine out of coincidence hits (s)", 
#        type = float)
parser.add_argument("num_bins", help="number of bins in final histogram", type = int) 
parser.add_argument("fluor_channel", help="channel of fluorescence detector", type = int)
parser.add_argument("scattered_channel", help="channel of scattered detector", type = int)
parser.add_argument("num_cores", help="number of cores to parallelize over", type = int)
#parser.add_argument("scattered_first", help="if True, coincident pairs start with scattered photons"
#        , type = ast.literal_eval)
#parser.add_argument("evt_offset", help="number of events to offset by", type = int)
#parser.add_argument("use_evt_offset", help="if true use evt_offset else use time_offset"
#        , type = ast.literal_eval)
args = parser.parse_args();

#where to get data from
directory = args.data_directory 
file_prefix = args.file_prefix
file_start = args.file_start
file_end = args.file_end
file_ext = args.file_ext

#where to write output to
name = args.write_name
writedir = args.write_directory

#processing parameters
time_window = args.time_bin
ROI_begin = args.fluor_begin
ROI_end = args.fluor_end
#time_offset = args.time_offset
num_bins = args.num_bins
fluor_channel = args.fluor_channel
scattered_channel = args.scattered_channel
num_cores = args.num_cores
#scattered_first = args.scattered_first
#evt_offset = args.evt_offset
#use_evt_offset = args.use_evt_offset


#########################################
#       End of User Input Region        #
#########################################

file_end+=1 #python loops are over the range [start, end)

def process_file(f_num, evt_offset):
    """Process one file according to the parameters scattered_first and use_evt_offset"""
    lc = lib_coincidence.Coincidence(directory, file_prefix, file_ext, time_window, 0,
        fluor_channel, scattered_channel, ROI_begin, ROI_end, num_bins, 0)
    #if use_evt_offset:
    #    return lc.evt_scattered_first(f_num, fluor_channel) if scattered_first \
    #            else lc.evt_fluor_first(f_num, fluor_channel)
    #else:
    #    return lc.process_scattered_first(f_num) if scattered_first \
    #            else lc.process_fluor_first(f_num)
    return lc.bin_then_and(f_num, time_window)

#map over all files
r = range(file_start, file_end)

#reduce all data
#result is stored in numpy array of tuples (in, out) where in and out are numpy arrays of 
#size num_bins. Each tuple contains data from one file.
to_reduce = Parallel(n_jobs = num_cores)(delayed(process_file)(i, i) for i in r)

#sum up in coincidence and out of coincidence counts
#returns tuple of np.array: (in, out)
def sum_data(acc_in, acc_out, in_cts, out_cts):
    return (acc_in + in_cts, acc_out + out_cts)

#def sum_data(acc_in, acc_out, in_cts, out_cts):
#    return (np.append(acc_in, np.sum(in_cts)) , np.append(acc_out, np.sum(out_cts)))

#sum up results
totals = reduce(lambda x, y: sum_data(*x+y[0]), to_reduce, (0,0)) 
fluor_cts = reduce(lambda x, y: x + y[1], to_reduce, 0)

def write_totals(totals):
    """Write out a .csv file of the summed distribution.
    Totals is a tuple of numpy arrays, where the first element is in coincidence
    and the second is out of coincidence"""

    write_file = open(writedir + name + '.csv', 'w')
    write_file.write('in_coincidence, out_coincidence\n') #column headers
    for i in range(0,totals[0].size):
        write_file.write(str(totals[0][i]) + ',' + str(totals[1][i]) + '\n')
    write_file.close()

    fluor = open(writedir + name + '-fluor.csv', 'w')
    fluor.write(str(fluor_cts))
    fluor.close()

def write_log():
    """write a log file for the run"""

    f = open(writedir + name + '.log', 'w')
    f.write('OSPS Log File \n')
    f.write('Log file for ' + writedir + name + '.csv \n')
    f.write('\n')
    f.write('Running on folder ' + directory + ', files ' + 
            str(file_start) + ' to ' + str(file_end - 1) + '\n')
    f.write('\n')
    f.write('Runtime parameters: \n')
    f.write('time bin: ' + '{:.3e}'.format(time_window) + '\n')
    #f.write('time offset: ' + '{:.3e}'.format(time_offset) + '\n')
    f.write('fluorescence region: [' + str(ROI_begin) + ',' + str(ROI_end) +
            ']\n')
    f.write('scattered channel: ' + str(scattered_channel) + '\n')
    f.write('fluorescence channel: ' + str(fluor_channel) + '\n')
    f.write('number of histogram bins: ' + str(num_bins) + '\n')
    #f.write('scattered_first: ' + str(scattered_first) + '\n')

    f.close()

write_totals(totals)
write_log()
