# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:56:01 2019
@author: jansenc3

- This python script converts binary files from the mass spectrometer to .txt files
- Works for "MID - vs time" measurement data. (NO rga files (yet))
- Converts all files in a single folder (and keeps the original)

HOW TO USE:
    
1. Copy/paste path to folder with your files to variable 'folder'.
2. Change backslashes to forward slashes, and make sure to add one at the end.
    Example: folder = 'C:/Users/jansenc3/surfdrive/Python/Convert_bin_to_txt/'
3. Change file extension of the input file if needed.
    Example: file_extension = '.mdc'
4. Run script
"""

folder = 'C:/Users/jansenc3/surfdrive/DATA/2022/07 Jul/220728/KW/'
file_extension = '.mdc'

"""
--------------------------------------------------------
DON'T CHANGE ANYTHING BELOW, UNLESS SCRIPT DOES NOT WORK
--------------------------------------------------------
"""
import numpy as np
import struct
import glob

def main():  
    filelist = glob.glob(folder+'*'+file_extension) #Find all files with given extension in folder
    if len(filelist) == 0:
        print('No files found in folder!')
    for filename in filelist:
        convert_file(filename)        
#End main

def convert_file(filename):
    #Variables
    nCycles =       {'pos': 100, 'val':-1, 'len':4, 'typ':'i'} #number of measurement cycles
    nDatablocks =   {'pos': 112, 'val':-1, 'len':2, 'typ':'h'} #number of datablocks (store different types of data)
    nChannels1 =    {'pos': 203, 'val':0, 'len':2, 'typ':'h','unit':''} #in TPD file: stores external (temperature) channel. In KW file: stores mass channels
    nChannels2 =    {'pos': 353, 'val':0, 'len':2, 'typ':'h','unit':''} #in TPD file: stores mass channels
    StartTime =     {'pos': 194, 'val':-1, 'len':4, 'typ':'i'} #Start time of measurement, in seconds after 1-1-1970
#    DataFormat =    {'pos': 201, 'val':-1, 'len':1} #I don't know what this is
    CycleLength =   6 #apparently 6 bytes for standard columns present in file
    DataPosition =  200 #starting position of data without considering data blocks. will increase according to # of channels
    ChannelNames = ['Reltime'] #names of the data columns. other names will be appended later
    
    #load file into python
    file_object = open(filename,'rb')
    dataset = file_object.read()
    file_object.close()
    filename = filename.replace(file_extension,'')#remove the extension from the filename
    
    #read several values from file
    for var in [nCycles, nDatablocks, nChannels1, StartTime]:
        #var['val'] = int.from_bytes(dataset[var['pos']:var['pos']+var['len']],byteorder='little', signed=True) #different way of reading data, only in python 3
        var['val'] = struct.unpack('<'+var['typ'], 
                       dataset[var['pos']:var['pos']+var['len']])[0] #[0] because function returns a tuple
    
    print("nCycles=", nCycles['val'])
    if nDatablocks['val'] > 2: 
        print('Too many datablocks - not supported (yet)!')
        return
    elif nDatablocks['val'] == 2: 
        #nChannels2['val'] = int.from_bytes(dataset[nChannels2['pos']:nChannels2['pos']+nChannels2['len']],byteorder='little', signed=True)
        nChannels2['val'] = struct.unpack('<'+nChannels2['typ'], 
                       dataset[nChannels2['pos']:nChannels2['pos']+nChannels2['len']])[0]      

    #Read names of the channels (four characters). 
    for nChannels in [nChannels1, nChannels2]: 
        nChannels['unit'] = (chr(dataset[nChannels['pos']+15])
                            +chr(dataset[nChannels['pos']+16]).replace('\x00','')) #reads two bytes but replaces null bytes if there are any
        for k in range(nChannels['val']): 
            pos = nChannels['pos']+114+33*k
            ChannelName = ''
            for j in range(5):
                ChannelName += chr(dataset[pos+j])
            ChannelNames.append(ChannelName.replace('\x00','')+'['+nChannels['unit']+']')         
    
    #Calculate start position of actual data from previous values, and length of each line/cycle       
    DataPosition += 117*nDatablocks['val'] + 33*(nChannels1['val']+nChannels2['val'])  
    CycleLength += 4*(nChannels1['val']+nChannels2['val']) #4 bytes per datapoint
    
    #make data array of the right size 
    data = np.full([nCycles['val'], len(ChannelNames)],np.nan) 
    
    for k in range(nCycles['val']): #read every cycle
        pos = DataPosition + CycleLength * k #this is redundant, but clear. can be replaced with pos = DataPosition before this for loop
        #Read time data, is stored in two bytes
        data[k][0] = (struct.unpack('<i',dataset[pos:pos+4])[0] 
                    + struct.unpack('<h',dataset[pos+4:pos+6])[0]/10000)
        pos += 6 #to get new position
        #Read channel data
        for l in range(len(ChannelNames)-1): #-1 because reltime is read separately
            data[k][l+1] = struct.unpack('<f',dataset[pos:pos+4])[0]
            pos += 4
    
    #Save data to file
    np.savetxt(filename+'.txt',data, comments='',
               header='Start time (s) =\n'+str(StartTime['val'])+'\n'+str(ChannelNames)[1:-1])
    print ('Converted '+filename.replace(folder[:-1]+'\\','')) #show which files are converted        
#End function convert_file   
    
main()