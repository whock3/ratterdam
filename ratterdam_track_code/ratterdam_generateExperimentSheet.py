# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:27:49 2019

@author: whockei1

Script to generate (most of) an experimental session record sheet

Configuration, log file, and events nev files are used.
Info on TTs, CSCs, refs, and events are printed in sections
Weird ts jump on first acquisition is accounted for and removed

Some blank lines printed at end for notes

To enter which refs you use, search this script for "EDIT ME" and edit the line
below with what refs you used. 

Run script entering the following in a python terminal:
    %run ratterdam_generateExperimentSheet.py [full path to your data directory]
Saves in the data directory you provide as a .txt file.
    
You must have python installed as well as numpy (datetime, os, and sys are builtin but also used)

The sheet gets pretty long so I recommend printing it on both sides of a sheet,
two pages per sheet, on your printer options. That way it's only like 3pages printed
"""

from sys import argv
import os, datetime, numpy as np

#script, data_directory = argv
data_directory  =  "E:\\Ratterdam\\R781\\R781BRD2"


def readNEV(datafile):
    "Reads a Neuralynx Events File (.nev) and returns the header and events separately"
    NEVfile = datafile + 'Events.nev'
    f = open(NEVfile,'rb')
    hdr = f.read(16*1024)
    dt =  np.dtype([('nstx',np.int16),('npkt_id',np.int16),('npkt_data_size',np.int16),('qwTimeStamp',np.uint64),('nevent_id',np.int16),('nttl',np.int16),('ncrc',np.int16),('ndummy1',np.int16),
               ('ndummy2',np.int16),('dnExtra',np.int32,8),('EventString',np.string_,128)])
    rec = np.fromfile(f,dtype=dt)
    f.close()
    return hdr,rec


def getFileName(directory, extension):
    """
    Given a directory, find all files with the given extension and return their name(s)
    """
    files = []
    for f in os.listdir(data_directory+"\\ConfigurationLog\\"):
        if f.endswith(extension):
            files.append(f)
            
    if len(files) == 1:
        return files[0]
    else:
        return files
    
def getTimeZero():
    found_time = False
    with open(data_directory+"\\CheetahLogFile.txt", "r") as logFile:
        while not found_time:
            line = next(logFile)
            if "Time Opened" in line:
                found_time = True
                ts_str = line.split(" ")[-1].rstrip() # removes trailing newline and ts is last entry in line
                timeZero = datetime.datetime.strptime(ts_str, "%H:%M:%S") # timeZero year and day will be Jan 1 1900 so obv dont use
    return timeZero
                
def getTsJumpValue():
    """
    Apparently the 1st time you acquire data in a cheetah session the NL timestamp
    jumps by minutes to days (!) in value. So find that jump and return it so
    it can be corrected later.
    
    Returns the jump in us. So subtract it off observed NL times to get 'real' one. Probably some slight offset
    """
    found_first_acquisition = False
    with open(data_directory+"\\CheetahLogFile.txt", "r") as logFile:
        lines = logFile.readlines()

    for i,line in enumerate(lines):
        if "Acquisition Started" in line and found_first_acquisition == False:
            found_first_acquisition = True # there can be, and usually are, multiple acquisitions but by looping and toggling you get 1st
            idx = i
    real_jump = (datetime.datetime.strptime(lines[idx].split(" ")[5], "%H:%M:%S") - datetime.datetime.strptime(lines[idx-1].split(" ")[5], "%H:%M:%S")).total_seconds()*1e6 # cheetah 6 only has second resolution anyway
    fictive_jump = (int(lines[idx].split(" ")[7]) - int(lines[idx-1].split(" ")[7])) - real_jump # so the big jump on 1st acquisition isn't real but there's usually some real gap btwn last log entry and 1st acq so take that into account
    return fictive_jump
    

def convertToClock(timeZero, ts_micros):
    """
    Convert a time in micros to h/m/s based on timeZero
    """
    ts_new = timeZero + datetime.timedelta(microseconds=np.float64(ts_micros))
    ts_new_str = f"{ts_new.hour}:{ts_new.minute}:{ts_new.second}"
    return ts_new_str

   
##############################
# Parameters To Be Extracted # 
##############################
        
tt_params = ["SetChannelNumber",
             "SetInputRange",
             "SetDspLowCutFrequency",
             "SetDspLowCutNumberTaps",
             "SetDspHighCutFrequency",
             "SetDspHighCutNumberTaps",
             "SetSpikeThreshold",
             "SetSpikeRetriggerTime",
             "SetSpikeAlignmentPoint",
             "SetSubChannelEnabled" #there are four copies of this parameter - one for each wire. B/c program scans lines it gets all but this isn't explicit
            ]

csc_params = ["SetChannelNumber",
              "SetInputRange",
              "SetDspLowCutFrequency",
              "SetDspLowCutNumberTaps",
              "SetDspHighCutFrequency",
              "SetDspHighCutNumberTaps"        
            ]

###########################################
# Functions to extract each piece of data #
###########################################

def extractTT(tt):
    data = []
    config_file.seek(0)
    found_block = False
    end_block = False
    while not end_block:
        while not found_block:
            line = next(config_file)
            if "Acquisition Entity" in line and f"TT{tt}" in line:
                found_block = True
        line = next(config_file)
        for param in tt_params:
            if param in line:
                data.append(line)
        if line == '\n':
            end_block = True # block of tt info ends with a newline
    return data
    

def extractCSC(csc):
    data = []
    config_file.seek(0)
    found_block = False
    end_block = False
    while not end_block:
        while not found_block:
            line = next(config_file)
            if "Acquisition Entity" in line and f"CSC{tt}" in line:
                found_block = True
        line = next(config_file)
        for param in csc_params:
            if param in line:
                data.append(line)
        if line == '\n':
            end_block = True # block of tt info ends with a newline
    return data

def extractRef(ref):
    data = []
    config_file.seek(0)
    found_block = False
    end_block = False
    while not end_block:
        while not found_block:
            line = next(config_file)
            if "Acquisition Entity" in line and ref in line:
                found_block = True
        line = next(config_file)
        for param in tt_params: #use same params at tt, the names are same they just often have 1 val bc ref wires are shorted together
            if param in line:
                data.append(line)
        if line == '\n':
            end_block = True # block of tt info ends with a newline
    return data

def extractCamera():
    data = []
    config_file.seek(0)
    found_block = False
    end_block = False
    while not end_block:
        while not found_block:
            line = next(config_file)
            if "Video Capture" in line:
                found_block = True
        line = next(config_file)
        if line == '\n':
            end_block = True # block of tt info ends with a newline
        else:
            data.append(line)
    return data


with open(data_directory+"\\{}_experimentSheet.txt".format(data_directory.split('\\')[-1]), 'w') as expSheet:
    
    
    ###############################
    ## Cheetah Configuration File #
    ###############################
    
    config_file_location = getFileName(data_directory, ".cfg")
    config_file = open(data_directory+f"\\ConfigurationLog\\{config_file_location}", 'r')
    
    nTT = 16 # number tetrodes
    nCSCs = 16 # number CSC channels. Should be same as nTT in all cases?
    
    #EDIT ME
    refs = ["HS1R1", "HS2R1"]
    
    expSheet.write("Tetrode, CSC, and Camera Configuration\n\n")
    
    for tt in range(1,nTT+1):
        
        tt_lines = extractTT(tt)
        
        expSheet.write(f"Tetrode {tt}\n")
        for line in tt_lines:
            expSheet.write(line)
            
        expSheet.write("----------------------------------\n")
        
    for csc in range(1,nCSCs+1):
        
        csc_lines = extractCSC(csc)
        
        expSheet.write(f"CSC {csc}\n")
        for line in csc_lines:
            expSheet.write(line)
        
        expSheet.write("----------------------------------\n")
        
    for ref in refs:
        
        ref_lines = extractRef(ref)
        
        expSheet.write(f"Ref {ref}\n")
        for line in ref_lines:
            expSheet.write(line)
        
        expSheet.write("----------------------------------\n")
        
    
    camera_lines = extractCamera()
    
    expSheet.write("Camera Data\n")
    for line in camera_lines:
        expSheet.write(line)
    
    expSheet.write("----------------------------------\n")
    
    config_file.close()
    
    ##########
    # Events #
    ##########
    fictive_jump = getTsJumpValue()
    timeZero = getTimeZero()
    hdr, rec = readNEV(data_directory+"\\")    
    
    expSheet.write("Events Data\n")
    for line in rec:
        if b"TTL" not in line[10]:
            expSheet.write(f"{convertToClock(timeZero, line[3]-fictive_jump)}: {line[10]}\n")
    expSheet.write("---------------------------------\n")
    
    expSheet.write("\n\n")
    expSheet.write("NOTES\n")
    for i in range(25):
        expSheet.write("_____________________________________\n")
    
                       
print("Experiment Sheet for {} Created".format(data_directory.split('\\')[-1]))
          
    
            
            
