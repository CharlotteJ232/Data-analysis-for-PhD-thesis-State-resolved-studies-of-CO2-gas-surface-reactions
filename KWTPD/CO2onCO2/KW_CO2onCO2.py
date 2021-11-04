import numpy as np
from matplotlib import pyplot as plt
import datetime as dt

################ General Parameters #################

folderstart = "C:/Users/jansenc3/surfdrive/"
# folderstart = "C:/Users/Werk/Surfdrive/"

#Measurement parameters
folder = '2021/10 Oct/211015/KW/'
file = 'kw01.txt'
onresonance = True
modulation_frequency = 1 #Hz


#Script
def main():
    
    remove_beginend()
    remove_background()
    remove_tilt()
    stack_data()

    
def remove_beginend():
    print('remove_beginend')

def remove_background():
    print('remove_background')

def remove_tilt():
    print('remove_tilt')

def stack_data():
    print('stack_data')

main()








