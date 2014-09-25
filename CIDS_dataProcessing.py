'''Program to process raw isdat intensities, and turn them into meaningful isotope ratios'''

import csv
import re

class CI:
  "A class for all the attributes of a single clumped isotope measurement"

  def __init__(self, name):
        self.name=name



class ACQUISITION:
  "A class for all the attributes of a single clumped isotope acquision"


  def __init__(self,acqNum):
      self.acqNum=acqNum
      self.voltSam=[]
      self.voltRef=[]
      self.background=[]
      self.d13Csample = 0
      self.d18Osample = 0
      self.d13Cref = 0
      self.d18Oref = 0

  def date(self,date):
    '''sets the date of the acquisiton'''
    self.date=date

  def time(self,time):
    '''sets the time of the acquition'''
    self.time=time




acqs=[]
CIDS=[]
acqIndex=0





fileReader = csv.reader(open(
  '/Users/Max/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014_justMax.csv', 'rU'), dialect='excel')
for line in fileReader:
  #storing all lines in a list for reference. Won't need this eventually
    CIDS.append(line)
    #finds a new acqusition, assuming 'pre' is in the line
    if 'pre' in line:
      cycle=0
      acqIndex += 1
      acqNum=int(line[line.index('pre')-1])
      acqs.append(ACQUISITION(acqNum))
      preIndex=line.index('pre')

    if cycle < 8:
      acqs[-1].voltSam.append(line[preIndex+1:preIndex+7])
      acqs[-1].voltRef.append(line[preIndex+7:preIndex+13])
      cycle += 1

    else:
      for s in line:
        if "Background:" in s:
          backgrounds=re.findall('[0-9]{0,20}\.[0-9]{0,10}', s)
          acqs[-1].background = [float(t) for t in backgrounds]

        
