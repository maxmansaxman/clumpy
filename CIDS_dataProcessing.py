'''Program to process raw isdat intensities, and turn them into meaningful isotope ratios'''

import csv
import re
import numpy as np

class CI:
  "A class for all the attributes of a single clumped isotope measurement"

  def __init__(self):
        self.name=''
        self.acqs=[]




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

def d13Ccalculator (samGas,refGas,refd13C):
  '''calculates d13C in VPDB'''
  # ideally, we'll write something that mimics isodat to correct raw ratios
  # to vpdb scale. Isodat is configured to do this using a method called
  # CO2_SSH, though CO2_craig is also an option.
  #For now, let's just take isodat's calculations
  vpdb=0.011237


filePath='/Users/Max/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014_justMax.csv'
def CIDS_parser(filePath)
samples=[]
CIDS=[]
acqIndex=0
sampleIndex=0

 # Initializing so that the first acquisitions are placed in the first sample
 # For now, naming the string by its number
samples.append(CI())



fileReader = csv.reader(open(
  '/Users/Max/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014_justMax.csv', 'rU'), dialect='excel')
for line in fileReader:
  #storing all lines in a list for reference. Won't need this eventually
    CIDS.append(line)
    if 'Counter' and 'Comments' in line:
      sampleIndex += 1
      samples.append(CI())

    #finds a new acqusition, assuming 'pre' is in the line
    if 'pre' in line:
      cycle=0
      acqIndex += 1
      acqNum=int(line[line.index('pre')-1]) # Not the same as acqIndex because acq are often skipped
      samples[-1].acqs.append(ACQUISITION(acqNum))
      preIndex=line.index('pre') #need to save this so that following lines know where to start pulling voltages
      if not samples[-1].name:
        try:
          samples[-1].name = line[line.index('SAMPLE')+1]
        except ValueError:
          samples[-1].name = line[line.index('SAMPLE:')+1]

    if cycle < 8:
      try:
        voltSams=[float(v) for v in line[preIndex+1:preIndex+7]]
        samples[-1].acqs[-1].voltSam.append(voltSams)
      except ValueError:
        pass
      voltRefs = [float(v2) for v2 in line[preIndex+7:preIndex+13]]
      samples[-1].acqs[-1].voltRef.append(voltRefs)
      cycle += 1

    elif 'Sample d13C (VPDB)' in line:
      samples[-1].acqs[-1].d13Csample=float(line[line.index('Sample d13C (VPDB)')+1])

    elif 'Sample d18O (SMOW)' in line:
      samples[-1].acqs[-1].d18Osample=float(line[line.index('Sample d18O (SMOW)')+1])

    elif 'Ref Gas d13C (VPDB)' in line:
      samples[-1].acqs[-1].d13Cref=float(line[line.index('Ref Gas d13C (VPDB)')+1])

    elif 'Ref Gas d18O (VSMOW)' in line:
      samples[-1].acqs[-1].d18Oref=float(line[line.index('Ref Gas d18O (VSMOW)')+1])

    else:
      for s in line:
        if "Background:" in s:
          backgrounds=re.findall('[0-9]{0,20}\.[0-9]{0,10}', s)
          samples[-1].acqs[-1].background = [float(t) for t in backgrounds]




# because of the way this code is written, an extra emtpy sample is added to the end of the list
if not samples[-1].name:
  samples.pop()

# for i in range(len(acqs)):
    # acqs[i].voltSam=np.asarray(acqs[i].voltSam)
    # acqs[i].voltRef=np.asarray(acqs[i].voltRef)
    # acqs[i].background=np.asarray(acqs[i].background)
