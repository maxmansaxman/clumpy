'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import CIDS_func

###############################################################

print('This script turns raw isodat files into a FLATLIST ')

analyses=[]
imported=[]

while True:
    sampleName = raw_input('Name of new sample, or press RETURN to stop: ')
    if len(sampleName) == 0:
        break
    analyses.append(CIDS_func.CI())
    analyses[-1].name=sampleName
    while True:
        acqName = raw_input('Drag an acq file for this sample, or press RETURN to stop: ')
        acqName=acqName.strip()
        if len(acqName) == 0:
            break
        acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0] #finds the acquision number from the file name
        acqNum=int(acqNum)
        if acqNum in imported:
            print('You already imported this file')
        else:
            imported.append(acqNum)
            analyses[-1].acqs.append(CIDS_func.ACQUISITION(acqNum))
            voltRef,voltSam,d13C,d18O = CIDS_func.Isodat_File_Parser(acqName)
            analyses[-1].acqs[-1].voltRef = voltRef
            analyses[-1].acqs[-1].voltSam = voltSam
            analyses[-1].acqs[-1].d13C_sample = d13C
            analyses[-1].acqs[-1].d18O_sample = d18O
            print('Acquisition '+str(acqNum)+ ' successfully imported.')

print('Acquisition imports complete')
print(str(len(analyses)) + ' analyses were imported')
print('Cleaning up analyses...')
analyses=CIDS_func.CIDS_cleaner(analyses)
print('Performing raw data reductions...')
analyses=CIDS_func.D47_calculations(analyses)
print('Exporting analyses to a flatlist...')
exportName = 'pythonFlatListExport'
CIDS_func.FlatList_exporter(analyses,exportName)
print('Analyses successfully exported')









#
#
# # For drag-and-drop capacity on windows systems
# # Ask for file name otherwise
# if (len(sys.argv) < 2):
#     fileName = raw_input('Name of the file to parse ? ')
# else:
#     fileName = sys.argv[1]
#
# FILE = open(FILENAME, 'r')
# FILENAME = FILENAME.rstrip('.txt')
#
#
# samples=CIDS_func.CIDS_parser(filePath)
#
# samples =CIDS_func.CIDS_cleaner(samples)
# samples=CIDS_func.D47_calculations(samples)
#
# pblChoice = raw_input('Would you like to do a pressure baseline correction? (y/n) ')
#
# if pblChoice == 'y':
#     print 'Pressure baseline correction selected'
#     print ''
