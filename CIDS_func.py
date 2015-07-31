''' Stores all the functions and classes used by various clumped isotope scripts'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import struct

class CI_VALUE(object):
    '''subclass defining how important isotopic ratios are calculated, and d18O_mineral'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if len(instance.voltRef)>6:
            if self.name in ['d45', 'd46', 'd47', 'd48', 'D47_raw', 'D48_raw']:
                return np.around(D47_calculation_valued(instance, self.name),3)
            elif self.name in ['d18O_min']:
                return carb_gas_oxygen_fractionation_acq(instance)

    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI calculation scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI value')

class CI_AVERAGE(object):
    '''subclass defining how isotopic ratios are averaged'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if len(instance.acqs)>=1:
            if self.name in ['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
            'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48','d48_stdev','D48_raw','D48_stdev']:
                return np.around(CI_averages_valued_individual(instance, self.name),3)

        else:
            raise ValueError('Sample has no acquisitions to average')


    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI averaging scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI average')


class CI(object):
    "A class for all the attributes of a single clumped isotope measurement"

    def __init__(self):
        self.name=''
        self.acqs=[]
        self.date=''
        self.type=''
        self.num=np.nan
        self.user=''
        self.skipFirstAcq = False
        self.TCO2 = np.nan
    #
    # def __getattr__(self,name):
    #     if not self.acqs:
    #         print 'No acquisitions for current sample'
    #
    #     elif name in ['d45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O',
    #     'd45_stdev','d46_stdev','d47_stdev','d48_stdev','D47_stdev','D47_sterr','D48_stdev','d13C_stdev','d18O_stdev','d18O_min']:
    #         CI_averages(self)
        #
        #
        # else:
        #   raise AttributeError, name
    d45 = CI_AVERAGE('d45')
    d46 = CI_AVERAGE('d46')
    d47 = CI_AVERAGE('d47')
    d48 = CI_AVERAGE('d48')
    D47_raw = CI_AVERAGE('D47_raw')
    D48_raw = CI_AVERAGE('D48_raw')
    d13C = CI_AVERAGE('d13C')
    d18O_gas = CI_AVERAGE('d18O_gas')
    d18O_min = CI_AVERAGE('d18O_min')

    d45_stdev = CI_AVERAGE('d45_stdev')
    d46_stdev = CI_AVERAGE('d46_stdev')
    d47_stdev = CI_AVERAGE('d47_stdev')
    d48_stdev = CI_AVERAGE('d48_stdev')
    D47_stdev = CI_AVERAGE('D47_stdev')
    D47_sterr = CI_AVERAGE('D47_sterr')
    D48_stdev = CI_AVERAGE('D48_stdev')
    d13C_stdev = CI_AVERAGE('d13C_stdev')
    d18O_stdev = CI_AVERAGE('d18O_stdev')



class ACQUISITION(object):
    "A class for all the attributes of a single clumped isotope acquision"


    def __init__(self,acqNum):
        self.acqNum=acqNum
        self.voltSam=[]
        self.voltRef=[]
        self.background=[]
        self.d13C = 0
        self.d18O_gas = 0
        self.d13Cref = 0
        self.d18Oref = 0
        self.date=''
        self.time=''

        self.D47_excel=0
        self.D48_excel=0
        self.d45_excel=0
        self.d46_excel=0
        self.d47_excel=0
        self.d48_excel=0


    d46=CI_VALUE('d46')
    d45=CI_VALUE('d45')
    D47_raw=CI_VALUE('D47_raw')
    d47=CI_VALUE('d47')
    D48_raw=CI_VALUE('D48_raw')
    d48=CI_VALUE('d48')
    d18O_min = CI_VALUE('d18O_min')


    # def __getattr__(self,name):    #the first time this is called, calculates all relevant values, and stores them. So, next time it's not called
    #     if name in ['d45','d46','d47','d48','D47_raw','D48_raw', 'd45_stdev',
    #             'd46_stdev','d47_stdev','d48_stdev','D47_stdev','D47_sterr','D48_stdev', 'd18O_min']:
    #         D47_calculation(self)
    #         carb_gas_oxygen_fractionation(self)
    #
    #     else:
    #         raise AttributeError, name
    #



def d13Ccalculator (samGas,refGas,refd13C):
    '''calculates d13C in VPDB'''
  # ideally, we'll write something that mimics isodat to correct raw ratios
  # to vpdb scale. Isodat is configured to do this using a method called
  # CO2_SSH, though CO2_craig is also an option.
  #For now, let's just take isodat's calculations
    vpdb=0.011237

def carb_gas_oxygen_fractionation(sample):
    '''calculates the d18O of a carbonate mineral from which the CO2 was digested'''
    # Done properly, this should be a function of the d18O of the gas, the reaction T,
    # and the cabonate phase of the sample, but we're sticking to 90C calcite for now'''

    vsmow_18O=0.0020052
    vpdb_18O = 0 #don't know this right now
    d18O_vpdb = (sample.d18O_gas-30.86)/1.03086
    sample.d18O_min = ((d18O_vpdb+1000)/1.00821)-1000
    return sample

def carb_gas_oxygen_fractionation_acq(acq):
    '''calculates the d18O of a carbonate mineral from which the CO2 was digested'''
    # Done properly, this should be a function of the d18O of the gas, the reaction T,
    # and the cabonate phase of the sample, but we're sticking to 90C calcite for now'''

    vsmow_18O=0.0020052
    vpdb_18O = 0 #don't know this right now
    d18O_vpdb = (acq.d18O_gas-30.86)/1.03086
    d18O_min = ((d18O_vpdb+1000)/1.00821)-1000
    return d18O_min


def Isodat_File_Parser(fileName):
    '''Reads in a .did file (Isodat acquisition file), returns the raw voltages for
    ref gas, sample gas, and the isodat-calculated d13C and d18O of the Acquisition'''

    f=open(fileName,'rb')
    try:
        buff = f.read()
    finally:
        f.close()

    #1. Getting raw voltage data
    #Searching for the start of the raw voltages
    start=buff.find('CIntensityData')
    keys=[]
    voltRef=[]
    voltSam=[]
    #Slightly elegant method: pulling out based on spacing after 'CIntensityData' ascii key
    # TODO: make this more flexible
    # TODO: Catch errors if wrong voltage sequence found
    startPreVolt=start+2304+168 #observed location of 'pre' cycle on ref gas side
    voltRef.append(struct.unpack('6d',buff[startPreVolt:(startPreVolt+6*8)]))
    for i in range(7):
        startRefVolt=start+52+i*164 #observed location of ref gas voltage cycles
        voltRef.append(struct.unpack('6d',buff[startRefVolt:(startRefVolt+6*8)]))

        startSamVolt=start+1344+i*160 #observed location of sample gas voltage cycles
        voltSam.append(struct.unpack('6d',buff[startSamVolt:(startSamVolt+6*8)]))

    #2. Getting d13C and d18O data for each cycle
    startEval=buff.find('CDualInletEvaluatedDataCollect') #rough guess of starting position
    d13C=[]
    d18O=[]
    # Exact position is not consistent, so running searching over a 200 byte range for the right start point
    # This hinges on recognition of pattern where an alternating sequence of 8-bit doubles of cycle number (starting with zero) and bulk isotopic composition for that number
    found13Cstart=False
    found18Ostart=False
    while i < 200 and not (found13Cstart and found18Ostart) :

        if not found13Cstart:
            start13C = startEval+720+i
            testList1=struct.unpack('ddd',buff[start13C:start13C+8] + buff[start13C+16:start13C+16+8] + buff[start13C+32:start13C+32+8])
            if testList1 == (1.0, 2.0, 3.0):
                found13Cstart=True

        if not found18Ostart:
            start18O = startEval+1000+i
            testList2=struct.unpack('ddd',buff[start18O:start18O+8] + buff[start18O+16:start18O+16+8] + buff[start18O+32:start18O+32+8])
            if testList2 == (1.0, 2.0, 3.0):
                found18Ostart=True
        i+=1

    # Alerting if one of the search sequences failed
    if not found13Cstart:
        print('Failed to find an appropriate byte sequence for d13C')

    if not found18Ostart:
        print('Failed to find the appropriate byte sequence for d18O')


    # Now, actually pulling bulk isotope data based on start13C and start18O variables

    for i in range(7):
        # start13C=newstartEval-1461+16*i #observed location of ref gas voltage cycles
        d13C.append(struct.unpack('d',buff[start13C-8+16*i:(start13C+16*i)])[0])

        # start18O=newstartEval-1188+16*i
        d18O.append(struct.unpack('d',buff[start18O-8+16*i:(start18O+16*i)])[0])

    # Averaging d13C and d18O for each cycle
    d13C_final = sum(d13C)/len(d13C)
    d18O_final = sum(d18O)/len(d18O)

    #3. Pulling out other auxiliary info
    # 3.1 Ref gas isotope composition
    startRefGas=buff.find('CEvalDataSecStdTransferPart')

    d13C_ref=struct.unpack('d',buff[startRefGas+203:startRefGas+203+8])[0]
    d18O_ref=struct.unpack('d',buff[startRefGas+423:startRefGas+423+8])[0]

    # 3.2 whether or not method is a CO2_multiply or a *_start
    firstAcq = False
    startMethod = buff.find('CDualInletBlockData')
    methodBlock = buff[startMethod-120:startMethod-20].decode('utf-16')
    # if 'CO2_multiply_16V' in methodBlock:
    #     firstAcq = False
    # if 'Contin_Start' in methodBlock:
    #     lastAcq = True
    #     #TODO: make this a more robust test
    if 'AL_Pump_Trans' in methodBlock:
        firstAcq = True

    # 3.3 sample name
    # Find start of block with sample name
    startName = buff.find('CSeqLineIndexData')
    # Rough guess of range where sample name is, accounting for a large variation in length
    nameBlock = buff[startName+200:startName+400].decode('utf-16')
    #Exact name based on locations of unicode strings directly before and after
    sampleName = nameBlock[(nameBlock.find('Background')+18):(nameBlock.find('Identifier 1')-2)]
    # Encode as ascii for consistency
    sampleName = sampleName.encode('ascii')

    # 3.4 background values
    #find start of block with background values
    # startBackground = buff.find('CISLScriptMessageData')
    # stopBackground = buff.find('CMeasurmentErrors')
    # #Note incorrect spelling of 'measurement' is intentional
    # backgroundBlock = buff[startBackground+80:stopBackground].decode('utf-16')

    return voltRef, voltSam, d13C_final, d18O_final, d13C_ref, d18O_ref, sampleName, firstAcq


def CIDS_parser(filePath):
    '''Reads in a .csv CIDS file, pulls out voltages, calculated bulk isotopes, and
    relevant sample information'''

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
        samples[-1].acqs[-1].d13C=float(line[line.index('Sample d13C (VPDB)')+1])
        samples[-1].acqs[-1].d46_excel=float(line[line.index('d46')+1])

      elif 'Sample d18O (SMOW)' in line[0:10]:
        samples[-1].acqs[-1].d18O_gas=float(line[line.index('Sample d18O (SMOW)')+1])
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
        del samples[-1]


    # checking for samples with not enough acqs, alerting if there are some
    temp=8
    lowAcqs= [k for k in samples if len(k.acqs)<temp]
    if not lowAcqs:
        print 'All samples have at least %d acquistions' % temp

    else:
        print 'Samples with too few acqs are:'
        print '\n'.join(['Sample '+ str(k.num)+ ': '+ k.name +' has ' + str(len(k.acqs)) + ' acquisitions' for k in lowAcqs])

    # converting the voltages and backgrounds to arrays so that they can be easily used to do calculations
    for i in range(len(samples)):
            for j in range(len(samples[i].acqs)):
                # converting voltagted to arrays, and doing 3 sig figs to match CIDS Sheet
                samples[i].acqs[j].voltSam=np.around(np.asarray(samples[i].acqs[j].voltSam),3)
                samples[i].acqs[j].voltRef=np.around(np.asarray(samples[i].acqs[j].voltRef),3)
                samples[i].acqs[j].background=np.asarray(samples[i].acqs[j].background)
                # rounding d13C and d18O of each acq to 3 sig figs to match CIDS sheet
                (samples[i].acqs[j].d13C,samples[i].acqs[j].d18O_gas) = np.around((samples[i].acqs[j].d13C,samples[i].acqs[j].d18O_gas),3)


    print 'All samples are cleaned, and voltages converted to arrays'

    return samples

def CIDS_subtractBackground(samples):
  '''Function to subtract the backgrounds off of all voltages'''

  for i in range(len(samples)):
    for j in range(len(samples[i].acqs)):
      background=samples[i].acqs[j].background
      samples[i].acqs[j].voltSam = samples[i].acqs[j].voltSam-background
      samples[i].acqs[j].voltRef = samples[i].acqs[j].voltRef-background

def D47_calculation(acq):
    '''Performs all the clumped isotope calculations for a single acq'''

    vpdb_13C=0.0112372 # values copied from CIDS spreadsheet
    vsmow_18O=0.0020052
    vsmow_17O=0.0003799
    lambda_17=0.5164

    k=acq
    R13_sa=(k.d13C/1000+1)*vpdb_13C
    R18_sa=(k.d18O_gas/1000+1)*vsmow_18O
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
    delta_measured_mean[1,:]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

    R_calculated_sample = (delta_measured_mean[0,:]/1000 + 1)*R_stoch[1,:]
    D_raw=(R_calculated_sample/R_stoch[0,:] - 1)*1000


    k.D47_raw=D_raw[2]-D_raw[0]-D_raw[1]
    k.D48_raw=D_raw[3]-D_raw[0]-D_raw[1]
    k.d47=delta_measured_mean[0,2]
    k.d48=delta_measured_mean[0,3]

    k.d45=delta_measured_mean[0,0]
    k.d46=delta_measured_mean[0,1]

    acq=k

    return acq



def D47_calculations(samples):
  '''Performs all the clumped isotope calculations on the raw voltages and bulk d13C and d18O values'''


  for i in range(len(samples)):

    for j in range(len(samples[i].acqs)):
    #   samples[i].acqs[j]=D47_calculation(samples[i].acqs[j])
      samples[i].acqs[j]=carb_gas_oxygen_fractionation(samples[i].acqs[j])

    # CI_averages(samples[i])

  return samples


def CI_averages(sample):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props1=['d45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O','d18O_min']

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev','d48','D48_raw','D48_stdev']

    acqsToUse = range(len(sample.acqs))

    if sample.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append([sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        sample.acqs[i].d48, sample.acqs[i].D47_raw,sample.acqs[i].D48_raw,sample.acqs[i].d13C,sample.acqs[i].d18O_gas, sample.acqs[i].d18O_min])
        # values[i,:]=[sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        # sample.acqs[i].d48, samle.acqs[i].D47_raw,sample.acqs[i].D48_raw,sample.acqs[i].d13C,sample.acqs[i].d18O_gas]

    if len(values) != 0:
        values = np.asarray(values)
        (sample.d45, sample.d46, sample.d47, sample.d48, sample.D47_raw, sample.D48_raw,
        sample.d13C, sample.d18O_gas, sample.d18O_min) = values.mean(axis=0)

        # (sample.d45, sample.d46, sample.d47, sample.d48, sample.D47_raw, sample.D48_raw,
        # sample.d13C, sample.d18O) = values.mean(axis=0)

        (sample.d45_stdev, sample.d46_stdev, sample.d47_stdev, sample.d48_stdev, sample.D47_stdev,
        sample.D48_stdev, sample.d13C_stdev,sample.d18O_stdev, temp) = values.std(axis=0,ddof=1)

        sample.D47_sterr=sample.D47_stdev/np.sqrt(values.shape[0])

def FlatList_exporter(samples,fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional flat list'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['User','date','Type','Sample ID','spec #\'s', 'd13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
    'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
    counter = 0
    if displayProgress:
        for item in samples:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev ])
            counter += 1
            if ((counter * 100)*100) % (len(samples)*100) == 0:
                print(str((counter*100)/len(samples)) + '% done')
    else:
        for item in samples:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev ])
    export.close()
    return

def CIDS_exporter(samples, fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional CIDS file.
    Python CIDS files are also importable, to load all important sample info
    back into the program'''

    export=open(fileName + '_python_CIDS.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for sample in samples:
        wrt.writerow(['__NewSample__'])
        for acq in sample.acqs:
            wrt.writerow(['__NewAcq__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam)[1]*['',]+acq.voltRef[0,:].tolist())
            for i in range(len(acq.voltSam)):
                wrt.writerow(['',''] + acq.voltSam[i,:].tolist() + acq.voltRef[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','ref gas d13C','ref gas d18O','name','type','d13C (vpdb)',
            'd18O_gas (vsmow)', 'd45','d46','d47', 'd48', 'D47_raw', 'D48_raw'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.d13Cref, acq.d18Oref, sample.name, sample.type, acq.d13C,
            acq.d18O_gas, acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw])

        wrt.writerow([])
        wrt.writerow(['', 'd45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O_gas'])
        for acq in sample.acqs:
            wrt.writerow(['',acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw, acq.d13C, acq.d18O_gas])
        wrt.writerow([])
        wrt.writerow(['','User','date','Type','Sample ID','spec #\'s','SkipFirstAcq','d13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
        'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
        wrt.writerow(['__SampleSummary__',sample.user, sample.date, sample.type, sample.name, sample.num, sample.skipFirstAcq, sample.d13C, sample.d13C_stdev,
        sample.d18O_gas, sample.d18O_min, sample.d18O_stdev, sample.d47, sample.d47_stdev, sample.D47_raw,
        sample.D47_stdev, sample.D47_sterr, sample.d48, sample.d48_stdev, sample.D48_raw, sample.D48_stdev])
        wrt.writerow(22*['---',])

    export.close()
    return

def CIDS_importer(filePath, displayProgress = False):
    '''Imports voltage and sample data that were exported using the CIDS_exporter function'''
    fileImport = open(filePath,'rU')
    fileReader = csv.reader(fileImport, dialect='excel')
    samples = []

    cycles_in_acqs = 7

    for line in fileReader:
        if '__NewSample__' in line:
            samples.append(CI())
            continue
        if '__NewAcq__' in line:
            samples[-1].acqs.append(ACQUISITION(0))
            continue
        if '__AcqVoltage__' in line:
            voltIndex = line.index('__AcqVoltage__')
            cycle = 0
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            samples[-1].acqs[-1].voltRef.append(voltRefs)
            continue
        if cycle < (cycles_in_acqs):
            voltSams=[float(vs) for vs in line[voltIndex+2:voltIndex+8]]
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            samples[-1].acqs[-1].voltSam.append(voltSams)
            samples[-1].acqs[-1].voltRef.append(voltRefs)
            cycle += 1
            continue
        if '__AcqInfo__' in line:
            acqIndex = line.index('__AcqInfo__')
            samples[-1].acqs[-1].acqNum = float(line[acqIndex + 1])
            samples[-1].acqs[-1].d13Cref = float(line[acqIndex + 2])
            samples[-1].acqs[-1].d18Oref = float(line[acqIndex + 3])
            samples[-1].acqs[-1].d13C = float(line[acqIndex + 6])
            samples[-1].acqs[-1].d18O_gas = float(line[acqIndex + 7])
            samples[-1].acqs[-1].voltRef = np.asarray(samples[-1].acqs[-1].voltRef)
            samples[-1].acqs[-1].voltSam = np.asarray(samples[-1].acqs[-1].voltSam)
            continue
        if '__SampleSummary__' in line:
            summaryIndex = line.index('__SampleSummary__')
            [samples[-1].user, samples[-1].date, samples[-1].type, samples[-1].name] = line[summaryIndex + 1:summaryIndex+5]
            samples[-1].num = int(line[summaryIndex + 5])
            samples[-1].skipFirstAcq = bool(line[summaryIndex + 6])

    fileImport.close()
    return samples




def Get_gases(samples):
    '''Finds which analyses are heated and equilibrated gases, and assigns them TCO2 values'''
    properNames = raw_input("Do all equilibrated gases have '25' in name? (y/n) ").lower()
    for item in samples:
        if 'BOC' in item.name.upper():
            if properNames == 'y':
                if '25' in item.name:
                    item.TCO2 = 25
                    item.type = 'eg'
                else :
                    item.TCO2 = 1000
                    item.type = 'hg'
            else :
                while True:
                    choice = raw_input('Is acq num: ' + str(item.num) + ' with name: ' + item.name + ' a (h)eated gas, an (e)quilibrated gas?, or (s)kip? ')
                    if choice.lower() == 'h':
                        item.TCO2 = 1000
                        item.type = 'hg'
                        break
                    elif choice.lower() == 'e':
                        item.TCO2 = 25
                        item.type = 'eg'
                        break
                    elif choice.lower() == 's':
                        print('Skipping this sample ')
                        break
                    else:
                        print('Not a valid entry, please try again')
    print('All heated and equilibrated gases have been found and labeled ')
    return

def Daeron_exporter(samples, fileName):
    '''Exports analyses in a csv that is formatted for Matthieu Daeron's xlcor47 script'''

    export=open(fileName +'_daeron'+ '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for item in samples:
        if np.isnan(item.TCO2):
            wrt.writerow([item.name, item.d47, item.D47_raw, item.D47_sterr, ])
        else :
            wrt.writerow([item.name, item.d47, item.D47_raw, item.D47_sterr, item.TCO2])
    export.close()
    return


def Pressure_Baseline_Processer(fileFolder):
    '''Pulls raw intensity data out of a folder of peak scans'''

    fileList = [i for i in os.listdir(fileFolder) if '.csv' in i] #copies all csv pbl files in directory
    print fileList
    day = '4-22-14' # TODO: make this I/O eventually
    fileList = [i for i in fileList if day in i] #only the csv files from certain days
    A=[]
    print fileList

    #First read in the csv files to numpy arrays
    for i in fileList:
        A.append([])
        fileName = fileFolder+i
        fileReader = csv.reader(open(fileName, 'rU'), dialect='excel')
        for line in fileReader:
            A[-1].append(line)
        A[-1].pop(0) #remove the first line with the headers

        A[-1]=np.asfarray(A[-1])
        # as a note, columns, in order, are: scanNumber, HV, mass44, mass45, mass46, mass47, mass48, mass49

    # Now, finding minimums
    int44=np.zeros((len(A),1))
    minimums=np.zeros((len(A),6))
    int49=np.zeros((len(A),1))

    for i in range(len(A)):
        max44=A[i][:,2].max()
        halfHeight=np.nonzero(A[i][:,2]>(max44/2)) #array of indices of all values with mass44 intensity is gt half of max
        peakCenter=int(np.median(halfHeight))
        halfWidth=halfHeight[0][-1]-peakCenter # currently this line is not necessary, but may be later on
        int44[i] = np.mean(A[i][peakCenter-5:peakCenter+5,2])
        int49[i] = np.mean(A[i][peakCenter-5:peakCenter+5,7])
        minimums[i,:]=A[i][peakCenter:,2:].min(axis=0) # minimum values of mass 44-49 on the right side of the peak





    return int44, int49, minimums

def D47_calculation_valued(acq, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    vpdb_13C=0.0112372 # values copied from CIDS spreadsheet
    vsmow_18O=0.0020052
    vsmow_17O=0.0003799
    lambda_17=0.5164

    R13_sa=(acq.d13C/1000+1)*vpdb_13C
    R18_sa=(acq.d18O_gas/1000+1)*vsmow_18O
    R17_sa=np.power((R18_sa/vsmow_18O),lambda_17)*vsmow_17O

    R13_ref=(acq.d13Cref/1000+1)*vpdb_13C
    R18_ref=(acq.d18Oref/1000+1)*vsmow_18O
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
    R_measured_sample=(acq.voltSam[:,1:6]/(np.tile(acq.voltSam[:,0],(5,1)).T))
    R_measured_ref=(acq.voltRef[:,1:6]/(np.tile(acq.voltRef[:,0],(5,1)).T))
    delta_measured=np.zeros(np.shape(R_measured_sample)) # Preallocating for size of delta array

    for l in range(len(R_measured_sample)):
        delta_measured[l,:]=(R_measured_sample[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
    # couches ratios in sample/std bracketing, put in delta notation

    delta_measured_mean=np.zeros((3,delta_measured.shape[1]))

    delta_measured_mean[0,:]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1,:]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

    R_calculated_sample = (delta_measured_mean[0,:]/1000 + 1)*R_stoch[1,:]
    D_raw=(R_calculated_sample/R_stoch[0,:] - 1)*1000


    d45=delta_measured_mean[0,0]
    d46=delta_measured_mean[0,1]
    d47=delta_measured_mean[0,2]
    d48=delta_measured_mean[0,3]
    D47_raw=D_raw[2]-D_raw[0]-D_raw[1]
    D48_raw=D_raw[3]-D_raw[0]-D_raw[1]


    calculatedCIValues = {'d45': d45, 'd46': d46, 'd47': d47, 'd48': d48, 'D47_raw': D47_raw, 'D48_raw': D48_raw}

    return calculatedCIValues[objName]

def CI_averages_valued(sample, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48', 'd48_stdev','D48_raw','D48_stdev']

    acqsToUse = range(len(sample.acqs))

    if sample.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append([sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        sample.acqs[i].d48, sample.acqs[i].D47_raw,sample.acqs[i].D48_raw,
        sample.acqs[i].d13C,sample.acqs[i].d18O_gas, sample.acqs[i].d18O_min])
        # values[i,:]=[sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        # sample.acqs[i].d48, samle.acqs[i].D47_raw,sample.acqs[i].D48_raw,sample.acqs[i].d13C,sample.acqs[i].d18O_gas]

    if len(values) != 0:
        values = np.asarray(values)
        (d45, d46, d47, d48, D47_raw, D48_raw,
        d13C, d18O_gas, d18O_min) = values.mean(axis=0)

        (d45_stdev, d46_stdev, d47_stdev, d48_stdev, D47_stdev,
        D48_stdev, d13C_stdev, d18O_stdev, temp) = values.std(axis=0,ddof=1)

        D47_sterr=D47_stdev/np.sqrt(values.shape[0])

        calculatedCIAverages = calculatedCIValues = {'d45': d45, 'd46': d46,
        'd47': d47, 'd48': d48, 'D47_raw': D47_raw, 'D48_raw': D48_raw,
        'd45_stdev': d45_stdev, 'd46_stdev': d46_stdev, 'd47_stdev': d47_stdev,
        'd48_stdev': d48_stdev, 'D47_stdev': D47_stdev, 'D48_stdev': D48_stdev,
        'd13C': d13C, 'd18O_gas': d18O_gas, 'd18O_min': d18O_min,
        'd13C_stdev': d13C_stdev, 'd18O_stdev': d18O_stdev, 'D47_sterr': D47_sterr}

        return calculatedCIAverages[objName]

def CI_averages_valued_individual(sample, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48','D48_raw','D48_stdev']

    valName = objName.replace('_stdev','')

    valName = valName.replace('_sterr','')

    if valName in ['D47', 'D48']:
        valName += '_raw'
    elif valName in ['d18O']:
        valName += '_gas'

    acqsToUse = range(len(sample.acqs))

    if sample.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append(getattr(sample.acqs[i],valName))
        # values.append([sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        # sample.acqs[i].d48, sample.acqs[i].D47_raw,sample.acqs[i].D48_raw,
        # sample.acqs[i].d13C,sample.acqs[i].d18O_gas, sample.acqs[i].d18O_min])
        # values[i,:]=[sample.acqs[i].d45,sample.acqs[i].d46,sample.acqs[i].d47,
        # sample.acqs[i].d48, samle.acqs[i].D47_raw,sample.acqs[i].D48_raw,sample.acqs[i].d13C,sample.acqs[i].d18O_gas]

    if len(values) != 0:
        values = np.asarray(values)
        if '_stdev' in objName:
            return values.std(axis=0,ddof=1)
        elif '_sterr' in objName:
            return values.std(axis=0,ddof=1)/np.sqrt(values.shape[0])
        else:
            return values.mean(axis=0)

def CI_comparer(samples1, samples2):
    '''compares two sets of CI data for equivalence in D47 calculations'''
    equalSoFar = True
    if not len(samples1) == len(samples2):
        print('Data sets are not the same length!')

    for i in range(len(samples1)):
        for j in range(len(samples1[i].acqs)):
            if not samples1[i].acqs[j].D47_raw == samples2[i].acqs[j].D47_raw:
                print('Acq num {1} from sample number {0} from comparison does not agree in D47'.format(i,j))
                print('analyses1 D47 = {0:.3f}, analyses2 D47 = {1:.3f}'.format(samples1[i].acqs[j].D47_raw, samples2[i].acqs[j].D47_raw))
                equalSoFar = False
        if not samples1[i].D47_stdev == samples2[i].D47_stdev:
            print('Sample number {0} from data sets do not agree in D47_stdev'.format(i))
            print('1. D47_stdev = {0:.3f}, 2. D47_stdev = {0:.3f}'.format(samples1[i].D47_stdev, samples2[i].D48_stdev))

    if equalSoFar:
        print('Data sets are equivalent')
