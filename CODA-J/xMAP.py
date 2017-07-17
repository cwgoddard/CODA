def read_netcdf(directory, f_name):
    """Read in a netcdf (.nc) file from disk and return the data contained
    in the form of a numpy structured array. 
    The file to be read is (directory + f_name) 
    The result is in the format:
        array([(time0, E0, channel0), (time1, E1, channel1), ...])
    The result is a structured array, so (e.g.) Energy data can be accessed by
    result['E']. The same holds for other columns of the data"""
    
    f = netcdf.netcdf_file(directory+f_name, 'r', mmap=False) #Load in netcdf file
    tmp_data = f.variables['array_data'] #get relevant data
    data = tmp_data[:].copy().astype('uint16') #copy data into memory from disk
    f.close() #close netcdf file
    
    data = data.ravel() #flatten data to 1 dimension 
    
    #Get number of events. 2^16 scales word 67 according to xMAP file format
    #(i.e. word 67 contains more significant digits of the num_events field)
    #word 66 and 67 are the header words that contain the number of events
    num_events =  data[66].astype('int32') + \
                      (2 ** 16) *(data[67].astype('int32'))
    
    time_constant = 20e-9 #conversion factor from xMAP time units to seconds
    
    #size of header, in words
    offset = np.array([256]).astype('int32') 
    
    #set up vectors to store data
    E = np.zeros(num_events)
    channel = np.zeros(num_events)
    time = np.zeros(num_events, dtype='float64')
    
    #keep track of how much data we have processed
    #start by skipping the header
    dynamic_offset = offset
    
    for i in range(0,num_events):
        dynamic_offset = offset + 3*i #we process 3 words each iteration
        
        #xMAP stores data in specific bits of these words
        word1 = data[dynamic_offset]
        word2 = data[dynamic_offset + 1]
        word3 = data[dynamic_offset + 2]
        
        #extract channel bits (13-15 of word 1)
        channel[i] = np.bitwise_and(np.right_shift(word1, 13), 3)
        #extract energy bits (0-12 of word 1)
        E[i] = np.bitwise_and(word1, 0x1fff)
        #extract time bits (word3 contains MSB)
        time[i] = (word3 * (2**16) + word2) * time_constant
                 
    #package data into table format (as numpy structured array)
    return np.array(list(zip(time, E, channel)),
            dtype = {'names':['time', 'E', 'channel'],
            'formats':['float64', 'int32', 'int32']})
