'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import CIDS_func
import os

###############################################################

print('This script turns raw isodat files into a FLATLIST ')

# list containing instances of CI class, each with their own acq objects
analyses=[]
# list containing numbers of acquisitions already imported
imported=[]

modeChoice=raw_input('(a)utomatic mode or (m)anual mode? ').lower()


if modeChoice == 'm':
    print('Manual mode selected')
    while True:
        sampleName = raw_input('Name of new sample, or press RETURN to stop: ')
        if len(sampleName) == 0:
            break
        analyses.append(CIDS_func.CI())
        analyses[-1].name=sampleName
        analyses[-1].num=len(analyses)
        while True:
            acqName = raw_input('Drag an acq file for sample ' + analyses[-1].name +', or press RETURN to stop: ')
            acqName=acqName.strip()

            if len(acqName) == 0:
                break
            acqName = acqName.strip('"')
            acqName = os.path.abspath(acqName)
            # acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
            acqNum=re.findall('[0-9]{4}',os.path.basename(acqName))[0] #finds the acquision number from the file name, no matter how long the path nor what it contains
            acqNum=int(acqNum)
            if acqNum in imported:
                print('You already imported this file')
            else:
                voltRef,voltSam,d13C,d18O,d13C_ref,d18O_ref,rawSampleName,lastAcq = CIDS_func.Isodat_File_Parser(acqName)
                if sampleName != rawSampleName :
                    print('Sample name: ' + analyses[-1].name + ' does not match name in file: ' + rawSampleName + ' ')
                    nameErrorChoice = raw_input('Are you sure you want to include this acquisition (y/n)? ')
                    if nameErrorChoice.lower() == 'n':
                        continue
                imported.append(acqNum)
                analyses[-1].acqs.append(CIDS_func.ACQUISITION(acqNum))
                analyses[-1].acqs[-1].voltRef = voltRef
                analyses[-1].acqs[-1].voltSam = voltSam
                analyses[-1].acqs[-1].d13C_sample = d13C
                analyses[-1].acqs[-1].d18O_sample = d18O
                analyses[-1].acqs[-1].d13Cref = d13C_ref
                analyses[-1].acqs[-1].d18Oref = d18O_ref


                print('Acquisition '+str(acqNum)+ ' successfully imported.')
        # Catching situation where sample name used in error, so no acqs imported
        if len(analyses[-1].acqs) == 0:
            print('No acqs imported for sample ' + analyses[-1].name)
            print('Deleting sample ' + analyses[-1].name +'...')
            del analyses[-1]

elif modeChoice == 'a':
    print('Automatic mode selected')
    while True:
        acqFolder = raw_input('Drag a folder containing all acquisitions: ').strip()
        acqFolder = acqFolder.strip('"')
        acqFolder = os.path.abspath(acqFolder)
        if os.path.exists(acqFolder):
            break
        else:
            print('Invalid folder, please try again ')

    acqList = [i for i in os.listdir(acqFolder) if i.endswith('.did')]
    if len(acqList) == 0:
        print("No acquisiton files ('.did') found in folder ")
        quit()
    acqList=sorted(acqList)
    startNum = raw_input('Number of first acquisition to be processed: ')
    stopNum = raw_input('Number of last acquisition to be processed: ')
    #convert to int first to remove any leading zeros
    startNum = int(startNum)
    stopNum = int(stopNum)

    #now, adding the proper number of leading zeros in order to get an exact match
    startName = 'Acquisition-' + (4-len(str(startNum)))*'0' + str(startNum) + '.did'
    stopName = 'Acquisition-' + (4-len(str(stopNum)))*'0' + str(stopNum) + '.did'

    startNumIndex = acqList.index(startName)
    stopNumIndex = acqList.index(stopName)

    lastAcq = False
    nextAcqIsNewSample = True

    for i in range(startNumIndex,stopNumIndex+1):
        acqName = acqFolder +'/'+ acqList[i]
        # Catches files with a size less than 123 kb and skips them
        if os.path.getsize(acqName) < 123000:
            print('Skipping acq num ' + str(acqList[i]) + 'because file too small')
            continue
        # Finds the acquision number from the file name, no matter how long the path nor what it contains
        acqNum=re.findall('[0-9]{4}',acqName.split('/')[-1])[0]
        acqNum=int(acqNum)
        # Actually processing the file
        print('Importing acq num ' + str(acqNum) + ' ')
        voltRef,voltSam,d13C,d18O,d13C_ref,d18O_ref,rawSampleName,lastAcq = CIDS_func.Isodat_File_Parser(acqName)
        # Creates a new sample if acquisition is not a 'CO2_multiply method'

        if nextAcqIsNewSample:
            analyses.append(CIDS_func.CI())
            analyses[-1].name = rawSampleName
            analyses[-1].num = len(analyses)
            print('Found new sample, with name: ' + rawSampleName)
            nextAcqIsNewSample = False
        # If acq is a continue_start, prepare for next acq to be a new sample
        # This must come after previous condition to account for situation where we start on a continue_start
        if lastAcq:
            nextAcqIsNewSample = True
        # Catches the rare case where an acquisition block starts with a 'CO2_multiply'
        # And checks whether this was intentional
        # if len(analyses) == 0:
        #     print("Acquisition block starts with a 'CO2_multiply' method ")
        #     blockStartChoice = raw_input('Are you sure you meant to start here? (y/n) ')
        #     if blockStartChoice.lower() == 'n':
        #         break
        #     elif blockStartChoice.lower() == 'y':
        #         analyses.append(CIDS_func.CI())
        #         analyses[-1].name = rawSampleName
        # Catching case where name in file does not match current acq name
        if analyses[-1].name != rawSampleName :
            print('Sample name: ' + analyses[-1].name + ' does not match name in file: ' + rawSampleName + ' ')
            nameErrorChoice = raw_input('(s)kip acquisition, (i)nclude it, or make a (n)ew sample from it? ')
            if nameErrorChoice.lower() == 's':
                print('Skipping acquisition ')
                continue
            elif nameErrorChoice.lower() == 'n':
                print('Making a new sample with name: ' + rawSampleName)
                analyses.append(CIDS_fund.CI())
                analyses[-1].name = rawSampleName
                analyses[-1].num = len(analyes)
            else:
                print('Including acquisition ')

        # if no errors caught above, actually add acquisition to analyses
        imported.append(acqNum)
        analyses[-1].acqs.append(CIDS_func.ACQUISITION(acqNum))
        analyses[-1].acqs[-1].voltRef = voltRef
        analyses[-1].acqs[-1].voltSam = voltSam
        analyses[-1].acqs[-1].d13C_sample = d13C
        analyses[-1].acqs[-1].d18O_sample = d18O
        analyses[-1].acqs[-1].d13Cref = d13C_ref
        analyses[-1].acqs[-1].d18Oref = d18O_ref
        print('Acquisition '+str(acqNum)+ ' successfully imported.')

else :
    print('Not a valid mode choice')
    print('Goodbye')

if len(analyses) != 0:
    print('Acquisition imports complete')
    print(str(len(analyses)) + ' analyses were imported')

    incFirstAcq = raw_input('Do you want to ignore the first acq of every sample (y/n)? ')
    if incFirstAcq.lower() == 'y':
        for i in range(len(analyses)):
            analyses[i].skipFirstAcq = True
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
