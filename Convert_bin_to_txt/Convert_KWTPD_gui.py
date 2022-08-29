# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
IMPORTANT:
I. This program only works in Python 3.
II. If using spyder (and possibly other editors): the program may do nothing 
every second time it's started. To prevent this (in spyder), choose to "execute
in a dedicated console" in the run menu (ctrl+f6 or f6, or in the menu bar)

If something does not work, please contact c.jansen@lic.leidenuniv.nl.
-------------------------------------------------------------------------------

@author: jansenc3

Made with the help of the tutorial:
http://zetcode.com/gui/pyqt5/
"""

import sys
import numpy as np
import struct
import glob

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


class Window(QWidget):
    
    def __init__(self):

        QApplication.setQuitOnLastWindowClosed(True)
#        QApplication.quit() 
        super().__init__()
        
#        self.initUI()
        
        vbox = QVBoxLayout()
        
        #Explanation text
        vbox.addWidget(QLabel('This program converts quadstar time domain \n' 
                              'data from .mdc to .txt. It works for files in which \n'
                              'only QMS data are recorded, but also for TPD\n'
                              'measurements with an external data channel.'))

        #Open button
        self.open_button = QPushButton('Select folder')
        self.open_button.setToolTip('All files with the given extension in this '
                                    'folder will be converted')
        self.open_button.clicked.connect(self.get_filelist)
        vbox.addWidget(self.open_button)

        #File extension box
        vbox.addWidget(QLabel('File extension'))
        self.ext = '.mdc'
        self.ext_text = QLineEdit('.mdc')
        self.ext_text.setToolTip('File extension')
        self.ext_text.setFixedSize(250, 20)
        self.ext_text.setDisabled(True)
        vbox.addWidget(self.ext_text)
        
        #List of files that can be converted       
        self.open_text = QLabel()
        self.open_text.setText('No folder selected')
        vbox.addWidget(self.open_text)

        #Convert button        
        self.convert_button = QPushButton('Convert files')
        self.convert_button.setToolTip('Files will be converted to .txt')
        self.convert_button.clicked.connect(self.convert_files)
        self.convert_button.setDisabled(True)
        vbox.addWidget(self.convert_button)
        
        #Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximum(100)
#        self.progress_timer = QBasicTimer()
#        self.step = 0
        vbox.addWidget(self.progress_bar)

        #Set window layout
        mainbox = QHBoxLayout()
        mainbox.addLayout(vbox)        
        self.setLayout(mainbox)        
        self.setGeometry(300, 300, 300, 220)
        self.setWindowTitle('Convert .mdc files')
        self.setWindowIcon(QIcon('web.png'))        

        self.show()
        
    def get_filelist(self):
        
        self.fn = QFileDialog.getExistingDirectory() #get folder name
        
        self.filelist = glob.glob(self.fn+'/*'+self.ext) #get all files in folder with the given extension

        filestr = ''
        if len(self.filelist)==0:
            filestr = 'No files with the given extension in this folder!'
            self.convert_button.setDisabled(True) #disable convert  button if there are no convertible files
        else:
            self.convert_button.setDisabled(False) #enable convert button if there are convertible files
            for i in range(len(self.filelist)):
                filestr+= self.filelist[i] + '\n'
        self.open_text.setText(filestr)    #show the list of files in the window
        
        self.progress_bar.setValue(0) #reset progress bar
        
    def convert_files(self):
        progress = 0
        for filename in self.filelist:
            progress+=1
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
            filename = filename.replace(self.ext,'')#remove the extension from the filename
            
            #read several values from file
            for var in [nCycles, nDatablocks, nChannels1, StartTime]:
                #var['val'] = int.from_bytes(dataset[var['pos']:var['pos']+var['len']],byteorder='little', signed=True) #different way of reading data, only in python 3
                var['val'] = struct.unpack('<'+var['typ'], 
                               dataset[var['pos']:var['pos']+var['len']])[0] #[0] because function returns a tuple
            
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
            self.progress_bar.setValue(int(progress*100/len(self.filelist)))
        #End function convert_file   

   
        
        
if __name__ == '__main__':    
    app = QApplication(sys.argv)
    ShowWindow = Window()
    app.exec_()
