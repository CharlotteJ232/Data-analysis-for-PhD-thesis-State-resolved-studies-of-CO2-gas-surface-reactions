# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:18:43 2019

@author: Charlotte
"""
import struct
import datetime
# Position of number of cycles: 100. RGA file

ind = 201
leng = 4

typ = {1:'b', 2:'h', 4:'i'} #so struct.unpack uses right datatype for number of bytes

file_object = open('tpd06.mdc','rb')
all_text = file_object.read()
print all_text[ind:ind+leng]
print ord(all_text[ind])

inte = struct.unpack('<'+typ[leng], all_text[ind:ind+leng])
inte = inte[0]

if ind == 194:
    print datetime.datetime.fromtimestamp(inte)

print inte