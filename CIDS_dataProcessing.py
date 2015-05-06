'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import CIDS_func

###############################################################


# For drag-and-drop capacity on windows systems
# Ask for file name otherwise
if (len(sys.argv) < 2):
    FILENAME = raw_input('Name of the file to parse ? ')
else:
    FILENAME = sys.argv[1]

FILE = open(FILENAME, 'r')
FILENAME = FILENAME.rstrip('.txt')


# TODO: get a filePath from user
filePath='/Users/Max/Box Sync/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014_justMax_precise.csv'
print 'The file we\'re processing is: \n' + filePath

samples=CIDS_func.CIDS_parser(filePath)

samples =CIDS_func.CIDS_cleaner(samples)
samples=CIDS_func.D47_calculations(samples)

pblChoice = raw_input('Would you like to do a pressure baseline correction? (y/n) ')

if pblChoice == 'y':
    print 'Pressure baseline correction selected'
    print ''
