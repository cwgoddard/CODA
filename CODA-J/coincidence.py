import xMAP
import numpy as np
import csv

def get_data(file):
# Read nc file and get data

    data = xMAP.read_netcdf(file)
    return data

def get_sim_data(file):
# Read csv file containing simulation and get data

    time = []
    energy = []
    channel = []

    with open(file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            time.append(float(row[0]))
            energy.append(int(float(row[1])))
            channel.append(int(row[2]))

    return np.array(zip(time, energy, channel), 
                    dtype = {'names':['time', 'E', 'channel'],
                    'formats':['float64', 'int32', 'int32']})


def get_channel_data(channelID, data):
# Get energies and times associated with a given channel

    energy = []
    time = []
    for i in range(data.shape[0]):
        if (channelID == data['channel'][i]):
            energy.append(data['E'][i])
            time.append(data['time'][i])
    return energy, time          

def set_time_window(event_time, time_window, time_window_type):
# Sets time window around a given event of a certain type (right-sided, left-sided, or centered)

    lower_time = 0 # lower limit of window
    upper_time = 0 # upper limit of window

    if (time_window_type == 0):
    # Set up right-sided time window
         lower_time = event_time
         upper_time = lower_time + time_window

    elif (time_window_type == 1):
    # Set up left-sided time window
         upper_time = event_time
         lower_time = upper_time - time_window

    elif (time_window_type == 2):
    # Set up centered time window
         offset = time_window/2.0
         lower_time = event_time - offset
         upper_time = event_time + offset

    return lower_time, upper_time

def measure_in_coincidence(data, scattered_channel, fluor_channel, time_window_info, ROI, analysis_type):
# Measure the number of in coincidence counts

    # Get energies and times in scattered and fluorescence channels
    scattered_energy, scattered_time = get_channel_data(scattered_channel, data)
    fluor_energy_raw, fluor_time_raw = get_channel_data(fluor_channel, data)
    
    # Get time window information
    time_window = time_window_info[0]
    time_window_type = time_window_info[1]

    # Get fluorescence region
    ROI_start = ROI[0]
    ROI_end = ROI[1]

    # Filter out events in fluorescence channel that are not in the fluorescence
    # region of inerest

    fluor_time = []

    for i in range(len(fluor_energy_raw)):
         if (fluor_energy_raw[i] >= ROI_start and fluor_energy_raw[i] <= ROI_end):
             fluor_time.append(fluor_time_raw[i])

    # Set up measurement variables
    coincidences = 0 # total number coincidences
    energy_coincidences = [] # energy of the scattered photon in the coincidences
    first_index = 0 # index of first event that was checked for coincidence
    last_index = -1 # index of last event that was checked for coincidence

    first = True

    # Get start and end times of exposure
    start = data['time'][0]
    end = data['time'][-1]

    last_position = 0

    if (analysis_type == 0):
        # Look for scattered photon
        for i in range(len(scattered_time)):
            found = False
            lower_time, upper_time = set_time_window(scattered_time[i], time_window, time_window_type)
            if (lower_time >= start and upper_time <= end):

                if (first):
                    first_index = i
                    first = False

                last_index = i
                for j in range(last_position, len(fluor_time)):
                    if (fluor_time[j] > lower_time and fluor_time[j] < upper_time):
                        found = True
                        last_position = j-1
                        break
            if (found):
                 coincidences += 1
                 energy_coincidences.append(scattered_energy[i])

    if (analysis_type == 1):
        # Look for fluorescence photon
        for i in range(len(fluor_time)):
            found =  False
            lower_time, upper_time = set_time_window(fluor_time[i], time_window, time_window_type)
            if (lower_time >= start and upper_time <= end):
               
                if (first):
                    first_index = i
                    first = False

                last_index = i
                for j in range(last_position, len(scattered_time)):
                    if (scattered_time[j] > lower_time and scattered_time[j] < upper_time):
                        found = True
                        last_position = j-1
                        break
            if (found):
                coincidences += 1
                energy_coincidences.append(scattered_energy[i])
    return coincidences, energy_coincidences, first_index, last_index


def get_map(channel_1_time, channel_2_time):
# Generates a list of indicies that correspond to the next closest event in channel 2
# given the time of the ith event in channel 1. Note that this includes events in the two
# channels that occur at the exact same time

    map = [] # List of indicies
    loc = 0 # Next index to search from in channel 2

    for i in range(len(channel_1_time)):
        for j in range(loc, len(channel_2_time)):
            if (channel_2_time[j] >= channel_1_time[i]):
                loc = j
                break
        map.append(j)

    return map



def measure_accidentals_with_time_offset(data, scattered_channel, fluor_channel, time_window_info, ROI, analysis_type, offset, first_index, last_index):
# Estimate the number of accidental counts using either the time offset algorithm

    # Get energies and times in scattered and fluorescence channels
    scattered_energy, scattered_time = get_channel_data(scattered_channel, data)
    fluor_energy_raw, fluor_time_raw = get_channel_data(fluor_channel, data)
    
    # Get time window information
    time_window = time_window_info[0]
    time_window_type = time_window_info[1]

    # Get fluorescence region
    ROI_start = ROI[0]
    ROI_end = ROI[1]

    # Filter out events in fluorescence channel that are not in the fluorescence
    # region of inerest

    fluor_time = []

    for i in range(len(fluor_energy_raw)):
         if (fluor_energy_raw[i] >= ROI_start and fluor_energy_raw[i] <= ROI_end):
             fluor_time.append(fluor_time_raw[i])

    # Set up measurement variables
    acoincidences = 0 # total number coincidences
    energy_acoincidences = [] # energy of the scattered photon in the coincidences

    # Get start and end times of exposure
    start = data['time'][0]
    end = data['time'][-1]

    
    # Estimate accidental counts using time offsets
    time_offset = offset # time offset in seconds
    last_position = 0

    if (analysis_type == 0):
        # Look for scattered photon
        for i in range(first_index, last_index+1):
            found = False
            lower_time, upper_time = set_time_window(scattered_time[i] + time_offset, time_window, time_window_type)
            if (lower_time >= start and upper_time <= end):
                for j in range(last_position, len(fluor_time)):
                    if (fluor_time[j] > lower_time and fluor_time[j] < upper_time):
                        found = True
                        last_position = j-1
                        break
                
                # If time offset places time window outside exposure, offset backwards instead
            else:
                lower_time, upper_time = set_time_window(scattered_time[i] - time_offset, time_window, time_window_type)
                for j in range(len(fluor_time)):
                    if (fluor_time[j] > lower_time and fluor_time[j] < upper_time):
                        found = True
                        break
            if (found):
                acoincidences += 1
                energy_acoincidences.append(scattered_energy[i])

    if (analysis_type == 1):
        # Look for fluorescence photon
        for i in range(first_index, last_index+1):
            found =  False
            lower_time, upper_time = set_time_window(fluor_time[i] + time_offset, time_window, time_window_type)
            if (lower_time >= start and upper_time <= end):
                for j in range(last_position, len(scattered_time)):
                    if (scattered_time[j] > lower_time and scattered_time[j] < upper_time):
                        found = True
                        last_position = j-1
                        break
                
            # If time offset places time window outside exposure, offset backwards instead
            else:
                lower_time, upper_time = set_time_window(fluor_time[i] - time_offset, time_window, time_window_type)
                for j in range(len(scattered_time)):
                    if (scattered_time[j] > lower_time and scattered_time[j] < upper_time):
                        found = True
                        break
            if (found):
                acoincidences += 1
                energy_acoincidences.append(scattered_energy[i])
    return acoincidences, energy_acoincidences


def measure_accidentals_with_event_offset(data, scattered_channel, fluor_channel, offset_channel, time_window_info, ROI, analysis_type, offset, first_index, last_index):
# Estimate accidentals with event offset algorithm

    # Get energies and times in scattered and fluorescence channels
    scattered_energy, scattered_time = get_channel_data(scattered_channel, data)
    fluor_energy_raw, fluor_time_raw = get_channel_data(fluor_channel, data)
#    _, offset_time = get_channel_data(offset_channel, data)

    # Get time window information
    time_window = time_window_info[0]
    time_window_type = time_window_info[1]

    # Get fluorescence region
    ROI_start = ROI[0]
    ROI_end = ROI[1]

    # Filter out events in fluorescence channel that are not in the fluorescence
    # region of inerest

    fluor_time = []

    for i in range(len(fluor_energy_raw)):
         if (fluor_energy_raw[i] >= ROI_start and fluor_energy_raw[i] <= ROI_end):
             fluor_time.append(fluor_time_raw[i])

    # If offset channel is the fluorescence channel, filter out flourescence events
    offset_time = []
    if (offset_channel == fluor_channel):
         for i in range(len(fluor_energy_raw)):
             if (fluor_energy_raw[i] < ROI_start or fluor_energy_raw[i] > ROI_end):
                 offset_time.append(fluor_time_raw[i])
    else:
        _, offset_time = get_channel_data(offset_channel, data)

    # Set up measurement variables
    acoincidences = 0 # total number of acoincidences
    energy_acoincidences = [] # energy of the scattered photon in the acoincidences

    # Get start and end times of exposure
    start = data['time'][0]
    end = data['time'][-1]

    
    # Estimate accidental counts using event offsets
    event_offset = offset # event offset
    
    if (analysis_type == 0):
        # Look for scattered photon
        map = get_map(scattered_time, offset_time)
        last_position = 0

        for i in range(first_index, last_index+1):
            found = False
            loc = map[i]
            if (loc + event_offset < len(offset_time)):
                lower_time, upper_time = set_time_window(offset_time[loc + event_offset], time_window, time_window_type)
            else:
                last_position = 0
                lower_time, upper_time = set_time_window(offset_time[loc - event_offset], time_window, time_window_type)
                
            if (lower_time < start or upper_time > end):
                last_position = 0
                lower_time, upper_time = set_time_window(offset_time[loc - event_offset], time_window, time_window_type)
            
            for j in range(last_position, len(fluor_time)):
                if (fluor_time[j] > lower_time and fluor_time[j] < upper_time):
                    found = True
                    last_position = j-1
                    break
            if (found):
                acoincidences += 1
                energy_acoincidences.append(scattered_energy[i])

    if (analysis_type == 1):
        # Look for fluorescence photon
        map = get_map(fluor_time_time, offset_time)
        last_position = 0

        for i in range(first_index, last_index+1):
            found = False
            loc = map[i]
            if (loc + event_offset < len(offset_time)):
                lower_time, upper_time = set_time_window(offset_time[loc + event_offset], time_window, time_window_type)
            else:
                last_position = 0
                lower_time, upper_time = set_time_window(offset_time[loc - event_offset], time_window, time_window_type)
                
            if (lower_time < start or upper_time > end):
                last_position = 0
                lower_time, upper_time = set_time_window(offset_time[loc - event_offset], time_window, time_window_type)
            
            for j in range(last_position, len(scattered_time)):
                if (scattered_time[j] > lower_time and scattered_time[j] < upper_time):
                    found = True
                    last_position = j-1
                    break
            if (found):
                acoincidences += 1
                energy_acoincidences.append(scattered_energy[i])
    return acoincidences, energy_acoincidences

def calculate_raw_spectra(energies, numbins, lowest_bin, largest_bin):
    # Calculates raw spectra and returns the bin centers and counts
    
    hist_energy = hist(energies, numbins, [lowest_bin, largest_bin])

    raw_spectra = hist_energy[0]
    bin_edges = hist_energy[1]
    bin_centers = bin_edges[-1] + 0.5*(bin_edges[1:] - bin_edges[:-1])

    return bin_centers, raw_spectra

def calculate_spectra_after_subtraction(energies, energies_background, numbins, lowest_bin, largest_bin):
    # Calculates background subtracted spectra
    
    bin_centers, raw_spectra = calculate_raw_spectra(energies, numbins, lowest_bin, largest_bin)
    _, background = calculate_raw_spectra(energies, numbins, lowest_bin, largest_bin)

    spectra = raw_spectra - spectra

    return bin_centers, spectra

def calculate_error(spectra):
    # Calculate bin error

    error = np.zeros(len(spectra))
    for i in range(len(spectra)):
        error[i] = np.sqrt(spectra[i])

    return error

