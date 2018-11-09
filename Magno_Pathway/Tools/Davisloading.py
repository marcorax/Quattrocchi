# -*- coding: utf-8 -*-
# @Author: marcorax93

''' 
This file contains functions used for Davis file loading 
'''
import numpy as np
import struct
import itertools


#These fanctions are supposed to deal with AEDAT 3.1 files.
# To know more and to understand what is going on in the lines underneat 
# please refer to page https://inivation.com/support/software/fileformat/#aedat-31

def skip_header(file_read):
    ''' skip header '''
    line = file_read.readline()
    while line.startswith(b'#'):
        if ( line == b'#!END-HEADER\r\n'):
            break
        else:
            line = file_read.readline()
            
            

def read_events(file_read, xdim, ydim):
    """ A simple function that read events from cAER tcp"""
    
    data = file_read.read(28)
    
    #raise exception if reached the end of the file
    if(len(data) == 0 ):
        return [-1], [-1], [-1], [-1], [-1], [-1]

    # read header
    eventtype = struct.unpack('H', data[0:2])[0]
    eventsource = struct.unpack('H', data[2:4])[0]
    eventsize = struct.unpack('I', data[4:8])[0]
    eventoffset = struct.unpack('I', data[8:12])[0]
    eventtsoverflow = struct.unpack('I', data[12:16])[0]
    eventcapacity = struct.unpack('I', data[16:20])[0]
    eventnumber = struct.unpack('I', data[20:24])[0]
    eventvalid = struct.unpack('I', data[24:28])[0]
    next_read = eventcapacity * eventsize  # we now read the full packet
    data = file_read.read(next_read)    
    counter = 0  # eventnumber[0]
    #return arrays
    x_addr_tot = []
    y_addr_tot = []
    pol_tot = []
    ts_tot =[]
    spec_type_tot =[]
    spec_ts_tot = []

    if(eventtype == 1):  # something is wrong as we set in the cAER to send only polarity events
        while(data[counter:counter + eventsize]):  # loop over all event packets
            aer_data = struct.unpack('I', data[counter:counter + 4])[0]
            timestamp = struct.unpack('I', data[counter + 4:counter + 8])[0]
            x_addr = (aer_data >> 17) & 0x00007FFF
            y_addr = (aer_data >> 2) & 0x00007FFF
            x_addr_tot.append(x_addr)
            y_addr_tot.append(y_addr)
            pol = (aer_data >> 1) & 0x00000001
            pol_tot.append(pol)
            ts_tot.append(timestamp)
            # print (timestamp, x_addr, y_addr, pol)
            counter = counter + eventsize
    
    #the special events are ultimately discarded (here they are loaded for debugging pourpuse)
    elif(eventtype == 0):
        spec_type_tot =[]
        spec_ts_tot = []
        while(data[counter:counter + eventsize]):  # loop over all event packets
            special_data = struct.unpack('I', data[counter:counter + 4])[0]
            timestamp = struct.unpack('I', data[counter + 4:counter + 8])[0]
            spec_type = (special_data >> 1) & 0x0000007F
            spec_type_tot.append(spec_type)
            spec_ts_tot.append(timestamp)
            if(spec_type == 6 or spec_type == 7 or spec_type == 9 or spec_type == 10):
                print (timestamp, spec_type)
            counter = counter + eventsize        


    return (np.array(x_addr_tot), np.array(y_addr_tot), np.array(pol_tot), np.array(ts_tot), np.array(spec_type_tot), np.array(spec_ts_tot))


def read_frames(file_read, xdim, ydim):
    """ A simple function that read events from cAER tcp"""
    
    data = file_read.read(28)
    
    #raise exception if reached the end of the file
    if(len(data) == 0 ):
        return [-1], [-1], [-1], [-1], [-1]

    # read header
    eventtype = struct.unpack('H', data[0:2])[0]
    eventsource = struct.unpack('H', data[2:4])[0]
    eventsize = struct.unpack('I', data[4:8])[0]
    eventoffset = struct.unpack('I', data[8:12])[0]
    eventtsoverflow = struct.unpack('I', data[12:16])[0]
    eventcapacity = struct.unpack('I', data[16:20])[0]
    eventnumber = struct.unpack('I', data[20:24])[0]
    eventvalid = struct.unpack('I', data[24:28])[0]
    next_read = eventcapacity * eventsize  # we now read the full packet
    data = file_read.read(next_read)    
    counter = 0  # eventnumber[0]
    pixelcounter = 36 
    #return arrays
    start_ts_tot = []
    end_ts_tot = []
    frame = []
    frames_tot = []
    spec_type_tot =[]
    spec_ts_tot = []
    if(eventtype == 2):  # something is wrong as we set in the cAER to send only frames
        while(data[counter:counter + eventsize]):  # loop over all frames packets
            while(pixelcounter<eventsize):
                pixel = struct.unpack('H', data[counter+pixelcounter:counter+ pixelcounter + 2])[0]
                frame.append(pixel)
                pixelcounter += 2
            pixelcounter = 36     
            frames_tot.append(np.array(frame).reshape(ydim,xdim))
            frame.clear()    
            timestamp = struct.unpack('I', data[counter + 4:counter + 8])[0]
            start_ts_tot.append(timestamp)
            timestamp = struct.unpack('I', data[counter + 8:counter + 12])[0]
            end_ts_tot.append(timestamp)
            # print (timestamp, x_addr, y_addr, pol)
            counter = counter + eventsize
        
    #the special events are ultimately discarded (here they are loaded for debugging pourpuse)
    elif(eventtype == 0):
        spec_type_tot =[]
        spec_ts_tot = []
        while(data[counter:counter + eventsize]):  # loop over all event packets
            special_data = struct.unpack('I', data[counter:counter + 4])[0]
            timestamp = struct.unpack('I', data[counter + 4:counter + 8])[0]
            spec_type = (special_data >> 1) & 0x0000007F
            spec_type_tot.append(spec_type)
            spec_ts_tot.append(timestamp)
            if(spec_type == 6 or spec_type == 7 or spec_type == 9 or spec_type == 10):
                print (timestamp, spec_type)
            counter = counter + eventsize        


    return (np.array(frames_tot), np.array(start_ts_tot), np.array(end_ts_tot), np.array(spec_type_tot), np.array(spec_ts_tot))

#import events from aedat 3.1 recorded with caer to Py
#It's less accessorate than the version on tili, but it's way quicker
#All the magic appens in read_events
def eventsToPy(filename, xdim, ydim):
    file_read = open(filename, "rb")
    skip_header(file_read)

    ts_events_tmp = []
    x_events_tmp = []
    y_events_tmp = []
    p_events_tmp = []
    while(1):
        x, y, p, ts_tot, spec_type, spec_type_ts = read_events(file_read, xdim, ydim)
        if(len(ts_tot) > 0 and ts_tot[0] == -1):   
            break
        ts_events_tmp.append(ts_tot)
        x_events_tmp.append(x)
        y_events_tmp.append(y)
        p_events_tmp.append(p)
    Events = np.zeros([4, len(list(itertools.chain(*ts_events_tmp)))])
    Events[0, :] = list(itertools.chain(*x_events_tmp))
    Events[1, :] = list(itertools.chain(*y_events_tmp))
    Events[2, :] = list(itertools.chain(*ts_events_tmp))
    Events[3, :] = list(itertools.chain(*p_events_tmp))
    file_read.close()
    return  (Events)

#import frames from aedat 3.1 recorded with caer to Py
#It's less accessorate than the version on tili, but it's way quicker
#All the magic appens in read_frames
def framesToPy(filename, xdim, ydim):
    file_read = open(filename, "rb")
    skip_header(file_read)

    frames_tmp = []
    start_ts_tmp = []
    end_ts_tmp = []
    
    while(1):
        frames, start_ts, end_ts, spec_type, spec_type_ts = read_frames(file_read, xdim, ydim)
        if(len(end_ts) > 0 and end_ts[0] == -1):   
            break
        frames_tmp.append(frames)
        start_ts_tmp.append(start_ts)
        end_ts_tmp.append(end_ts)
    frames = list(itertools.chain(*frames_tmp))    
    start_ts = list(itertools.chain(*start_ts_tmp))
    end_ts = list(itertools.chain(*end_ts_tmp))
    file_read.close()
    return  (frames,start_ts,end_ts)



#Final conversion of the events before being used as imputs of Brian2/teili nets
def fromNPtoBrian(Events, xRes, uniqueflag, milliflag):
    linEvents = np.zeros((3,len(Events[2])))
    linEvents[0] = Events[0]+Events[1]*xRes
    linEvents[1:3] = Events[2:4]
    posEvents = linEvents[:,linEvents[2]>0]
    negEvents = linEvents[:,linEvents[2]<1]
    if milliflag==1:
        posEvents[1]=(posEvents[1]/1000).astype(int)
        negEvents[1]=(negEvents[1]/1000).astype(int)
    if uniqueflag==1:
        posEvents=np.unique(np.transpose(posEvents),axis=0)
        negEvents=np.unique(np.transpose(negEvents),axis=0)
    return (np.transpose(posEvents), np.transpose(negEvents))
