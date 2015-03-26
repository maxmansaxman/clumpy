'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''

import csv
import re
import numpy as np
import matplotlib.pyplot as plt

class CI:
  "A class for all the attributes of a single clumped isotope measurement"

  def __init__(self):
        self.name=''
        self.acqs=[]
        self.date=''
        self.type=''
        self.num=0




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
      self.date=''
      self.time=''
      self.D47_raw=0
      self.d47=0
      self.D48=0
      self.d48=0
      self.D47_excel=0
      self.D48_excel=0
      self.d45_excel=0
      self.d46_excel=0
      self.d47_excel=0
      self.d48_excel=0
      self.d46=0
      self.d45=0

def d13Ccalculator (samGas,refGas,refd13C):
  '''calculates d13C in VPDB'''
  # ideally, we'll write something that mimics isodat to correct raw ratios
  # to vpdb scale. Isodat is configured to do this using a method called
  # CO2_SSH, though CO2_craig is also an option.
  #For now, let's just take isodat's calculations
  vpdb=0.011237

def CIDS_parser(filePath):
  samples=[]
  CIDS=[]
  acqIndex=0
  sampleIndex=0
  cycle=0

   # Initializing so that the first acquisitions are placed in the first sample
   # For now, naming the string by its number
  samples.append(CI())



  fileReader = csv.reader(open(filePath, 'rU'), dialect='excel')
  for line in fileReader:
    #storing all lines in a list for reference. Won't need this eventually
      CIDS.append(line)
      if 'Counter' and 'Comments' in line:
        sampleIndex += 1
        samples.append(CI())
        samples[-1].num=sampleIndex

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

      if 'Date:' in line:
        samples[-1].date= line[line.index('Date:')+1]
        samples[-1].acqs[-1].D47_excel=float(line[line.index('D47=')+1])

      elif 'Time:' in line:
        samples[-1].acqs[-1].time=line[line.index('Time:')+1]
        samples[-1].acqs[-1].D48_excel=float(line[line.index('D48=')+1])

      elif 'Analyst' in line:
        samples[-1].acqs[-1].d45_excel=float(line[line.index('d45')+1])

      elif 'Type' in line:
        samples[-1].type=line[line.index('Type')+1].lower()

      elif 'Sample d13C (VPDB)' in line[0:10]:
        samples[-1].acqs[-1].d13Csample=float(line[line.index('Sample d13C (VPDB)')+1])
        samples[-1].acqs[-1].d46_excel=float(line[line.index('d46')+1])

      elif 'Sample d18O (SMOW)' in line[0:10]:
        samples[-1].acqs[-1].d18Osample=float(line[line.index('Sample d18O (SMOW)')+1])
        samples[-1].acqs[-1].d47_excel=float(line[line.index('d47')+1])

      elif 'Ref Gas d13C (VPDB)' in line[0:10]:
        samples[-1].acqs[-1].d13Cref=float(line[line.index('Ref Gas d13C (VPDB)')+1])
        try:
          samples[-1].acqs[-1].d48_excel=float(line[line.index('d48')+1])
        except:
          print 'couldn\'t find d48 in acquistions %d' % acqIndex
          print line

      elif 'Ref Gas d18O (VSMOW)' in line[0:10]:
        samples[-1].acqs[-1].d18Oref=float(line[line.index('Ref Gas d18O (VSMOW)')+1])

      else:
        for s in line:
          if "Background:" in s:
            backgrounds=re.findall('[0-9]{0,20}\.[0-9]{0,10}', s)
            samples[-1].acqs[-1].background = [float(t) for t in backgrounds]

  return samples


def CIDS_cleaner(samples):
  '''function for cleaning up a parsed CIDS file and alerting to any corrupted samples'''

  # Cleaning up the parsed file
  # because of the way the parser is written, an extra empty sample is added to the end of the list
  if not samples[-1].acqs:
    samples.pop()

  # checking for samples with not enough acqs, alerting if there are some
  temp=8
  lowAcqs= [k for k in samples if len(k.acqs)<temp]
  if not lowAcqs:
    print 'All samples have at least %d acquistions' % temp

  else:
    print 'Samples with too few acqs are:'
    print '\n'.join(['Sample '+ str(k.num)+ ' has ' + str(len(k.acqs)) + ' acquisitions' for k in samples])

  # converting the voltages and backgrounds to arrays so that they can be easily used to do calculations
  for i in range(len(samples)):
    for j in range(len(samples[i].acqs)):
      samples[i].acqs[j].voltSam=np.asarray(samples[i].acqs[j].voltSam)
      samples[i].acqs[j].voltRef=np.asarray(samples[i].acqs[j].voltRef)
      samples[i].acqs[j].background=np.asarray(samples[i].acqs[j].background)

  print 'All samples are cleaned, and voltages converted to arrays'

  return samples

def CIDS_subtractBackground(samples):
  '''Function to subtract the backgrounds off of all voltages'''

  for i in range(len(samples)):
    for j in range(len(samples[i].acqs)):
      background=samples[i].acqs[j].background
      samples[i].acqs[j].voltSam = samples[i].acqs[j].voltSam-background
      samples[i].acqs[j].voltRef = samples[i].acqs[j].voltRef-background


def D47_calculations(samples):
  '''Performs all the clumped isotope calculations on the raw voltages and bulk d13C and d18O values'''

  vpdb_13C=0.0112372 # values copied from CIDS spreadsheet
  vsmow_18O=0.0020052
  vsmow_17O=0.0003799
  lambda_17=0.5164

  # TODO: remove this enrichment stuff
  enrichments=np.zeros((len(samples)*9,6))
  counter=0

  for i in range(len(samples)):
    for j in range(len(samples[i].acqs)):
      k=samples[i].acqs[j]
      R13_sa=(k.d13Csample/1000+1)*vpdb_13C
      R18_sa=(k.d18Osample/1000+1)*vsmow_18O
      R17_sa=np.power((R18_sa/vsmow_18O),lambda_17)*vsmow_17O

      R13_ref=(k.d13Cref/1000+1)*vpdb_13C
      R18_ref=(k.d18Oref/1000+1)*vsmow_18O
      R17_ref=np.power((R18_ref/vsmow_18O),lambda_17)*vsmow_17O

      # calculating stochastic ratios
      # to keep things organized and avoid repetetive code lines,
      # organizing calculations in arrays where item 1 is sample gas, item 2 is ref gas
      R13=np.array([R13_sa, R13_ref])
      R17=np.array([R17_sa, R17_ref])
      R18=np.array([R18_sa, R18_ref])


      R45_stoch=R13+2*R17
      R46_stoch=2*R13*R17+R17*R17+2*R18
      R47_stoch=2*R13*R18+2*R17*R18+R13*R17*R17
      R48_stoch=2*R17*R18*R13+R18*R18
      R49_stoch=R13*R18*R18

      R_stoch=np.array([R45_stoch, R46_stoch, R47_stoch, R48_stoch, R49_stoch]).T

      # calculating measured voltage ratios, with sample/ref bracketing
      R_measured_sample=(k.voltSam[:,1:6]/(np.tile(k.voltSam[:,0],(5,1)).T))
      R_measured_ref=(k.voltRef[:,1:6]/(np.tile(k.voltRef[:,0],(5,1)).T))

      delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

      for l in range(len(R_measured_sample)):
        delta_measured[l,:]=(R_measured_sample[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
        # couches ratios in sample/std bracketing, put in delta notation

      delta_measured_mean=np.zeros((3,delta_measured.shape[1]))

      delta_measured_mean[0,:]=np.mean(delta_measured, axis=0) # averaging for all cycles
      delta_measured_mean[1,:]=np.std(delta_measured, axis=0)  # standard deviation among all cycles
      delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

      R_calculated_sample = (delta_measured_mean[0,:]/1000 + 1)*R_stoch[1,:]
      D_raw=(R_calculated_sample/R_stoch[0,:] - 1)*1000

      k.D47=D_raw[2]-D_raw[0]-D_raw[1]
      k.D48=D_raw[3]-D_raw[0]-D_raw[1]
      k.d47=delta_measured_mean[0,2]
      k.d48=delta_measured_mean[0,3]

      k.d45=delta_measured_mean[0,0]
      k.d46=delta_measured_mean[0,1]

      excels=np.array([k.d45_excel, k.d46_excel, k.d47_excel, k.d48_excel, k.D47_excel, k.D48_excel])
      deltas=np.array([k.d45, k.d46, k.d47, k.d48, k.D47, k.D48])
      enrichments[counter,:]=(deltas-excels)

      counter += 1



  return samples, enrichments








###############################################################
# Here's where the actual I/O happens


# TODO: get a filePath from user
filePath='/Users/Max/Box Sync/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014_justMax.csv'
print 'The file we\'re processing is: \n' + filePath

samples=CIDS_parser(filePath)

samples =CIDS_cleaner(samples)
(samples,enrichments)=D47_calculations(samples)

bckChoice=raw_input('Would you like to correct the voltages for background? (y/n) ')
pblChoice = raw_input('Would you like to do a pressure baseline correction? (y/n) ')
