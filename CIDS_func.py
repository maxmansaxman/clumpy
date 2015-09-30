''' Stores all the functions and classes used by various clumped isotope scripts'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import time

global hg_slope
global hg_intercept
hg_slope = 0.0132
hg_intercept = -0.7869

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

class CI_CORRECTED_VALUE(object):
    '''subclass defining how isotopic ratios are averaged'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if -2 < instance.D47_raw < 2:
            if self.name in ['D47_CRF']:
                return np.around(CI_CRF_corrector(instance, self.name),3)

        else:
            raise ValueError('Sampe D47_raw is out of range')


    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI corretion scheme')


class CI(object):
    "A class for all the attributes of a single clumped isotope measurement"

    def __init__(self):
        self.name=''
        self.acqs=[]
        self.date = ''
        self.type=''
        self.num=np.nan
        self.user=''
        self.skipFirstAcq = False
        self.TCO2 = np.nan

    #
    # def __getattr__(self,name):
    #     if not self.acqs:
    #         print 'No acquisitions for current analysis'
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

    D47_CRF = CI_CORRECTED_VALUE('D47_CRF')



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

def carb_gas_oxygen_fractionation(analysis):
    '''calculates the d18O of a carbonate mineral from which the CO2 was digested'''
    # Done properly, this should be a function of the d18O of the gas, the reaction T,
    # and the cabonate phase of the analysis, but we're sticking to 90C calcite for now'''

    vsmow_18O=0.0020052
    vpdb_18O = 0 #don't know this right now
    d18O_vpdb = (analysis.d18O_gas-30.86)/1.03086
    analysis.d18O_min = ((d18O_vpdb+1000)/1.00821)-1000
    return analysis

def carb_gas_oxygen_fractionation_acq(acq):
    '''calculates the d18O of a carbonate mineral from which the CO2 was digested'''
    # Done properly, this should be a function of the d18O of the gas, the reaction T,
    # and the cabonate phase of the analysis, but we're sticking to 90C calcite for now'''

    vsmow_18O=0.0020052
    vpdb_18O = 0 #don't know this right now
    d18O_vpdb = (acq.d18O_gas-30.86)/1.03086
    d18O_min = ((d18O_vpdb+1000)/1.00821)-1000
    return d18O_min


def Isodat_File_Parser(fileName):
    '''Reads in a .did file (Isodat acquisition file), returns the raw voltages for
    ref gas, analysis gas, and the isodat-calculated d13C and d18O of the Acquisition'''

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
    # Rough guess of range where analysis name is, accounting for a large variation in length
    nameBlock = buff[startName+200:startName+400].decode('utf-16')
    #Exact name based on locations of unicode strings directly before and after
    analysisName = nameBlock[(nameBlock.find('Background')+18):(nameBlock.find('Identifier 1')-2)]
    # Encode as ascii for consistency
    analysisName = analysisName.encode('ascii')

    # 3.4 background values
    #find start of block with background values
    # startBackground = buff.find('CISLScriptMessageData')
    # stopBackground = buff.find('CMeasurmentErrors')
    # #Note incorrect spelling of 'measurement' is intentional
    # backgroundBlock = buff[startBackground+80:stopBackground].decode('utf-16')

    # 3.5 Date and Time
    # Find the start of Ctime block
    startTime = buff.find('CTimeObject')+49
    # Pull out the time_t time based on startTime location, seconds since the epoch (Jan 1st, 1970), GMT
    time_t_time = struct.unpack('i',buff[startTime:startTime+4])[0]
    # time_str = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(time_t_time))
    time_str = time.strftime('%m/%d/%Y', time.localtime(time_t_time))


    return voltRef, voltSam, d13C_final, d18O_final, d13C_ref, d18O_ref, analysisName, firstAcq, time_str


def CIDS_parser(filePath):
    '''Reads in a .csv CIDS file, pulls out voltages, calculated bulk isotopes, and
    relevant sample information'''

    analyses=[]
    CIDS=[]
    acqIndex=0
    analysisIndex=0
    cycle=0

   # Initializing so that the first acquisitions are placed in the first analysis
   # For now, naming the string by its number
    analyses.append(CI())



    fileReader = csv.reader(open(filePath, 'rU'), dialect='excel')
    for line in fileReader:
    #storing all lines in a list for reference. Won't need this eventually
      CIDS.append(line)
      if 'Counter' and 'Comments' in line:
        analysisIndex += 1
        analyses.append(CI())
        analyses[-1].num=analysisIndex

      #finds a new acqusition, assuming 'pre' is in the line
      if 'pre' in line:
        cycle=0
        acqIndex += 1
        acqNum=int(line[line.index('pre')-1]) # Not the same as acqIndex because acq are often skipped
        analyses[-1].acqs.append(ACQUISITION(acqNum))
        preIndex=line.index('pre') #need to save this so that following lines know where to start pulling voltages
        if not analyses[-1].name:
          try:
            analyses[-1].name = line[line.index('SAMPLE')+1]
          except ValueError:
            analyses[-1].name = line[line.index('SAMPLE:')+1]

      if cycle < 8:
        try:
          voltSams=[float(v) for v in line[preIndex+1:preIndex+7]]
          analyses[-1].acqs[-1].voltSam.append(voltSams)
        except ValueError:
          pass
        voltRefs = [float(v2) for v2 in line[preIndex+7:preIndex+13]]
        analyses[-1].acqs[-1].voltRef.append(voltRefs)
        cycle += 1

      if 'Date:' in line:
        analyses[-1].date= line[line.index('Date:')+1]
        analyses[-1].acqs[-1].D47_excel=float(line[line.index('D47=')+1])

      elif 'Time:' in line:
        analyses[-1].acqs[-1].time=line[line.index('Time:')+1]
        analyses[-1].acqs[-1].D48_excel=float(line[line.index('D48=')+1])

      elif 'Analyst' in line:
        analyses[-1].acqs[-1].d45_excel=float(line[line.index('d45')+1])

      elif 'Type' in line:
        analyses[-1].type=line[line.index('Type')+1].lower()

      elif 'Sample d13C (VPDB)' in line[0:10]:
        analyses[-1].acqs[-1].d13C=float(line[line.index('Sample d13C (VPDB)')+1])
        analyses[-1].acqs[-1].d46_excel=float(line[line.index('d46')+1])

      elif 'Sample d18O (SMOW)' in line[0:10]:
        analyses[-1].acqs[-1].d18O_gas=float(line[line.index('Sample d18O (SMOW)')+1])
        analyses[-1].acqs[-1].d47_excel=float(line[line.index('d47')+1])

      elif 'Ref Gas d13C (VPDB)' in line[0:10]:
        analyses[-1].acqs[-1].d13Cref=float(line[line.index('Ref Gas d13C (VPDB)')+1])
        try:
          analyses[-1].acqs[-1].d48_excel=float(line[line.index('d48')+1])
        except:
          print 'couldn\'t find d48 in acquistions %d' % acqIndex
          print line

      elif 'Ref Gas d18O (VSMOW)' in line[0:10]:
        analyses[-1].acqs[-1].d18Oref=float(line[line.index('Ref Gas d18O (VSMOW)')+1])

      else:
        for s in line:
          if "Background:" in s:
            backgrounds=re.findall('[0-9]{0,20}\.[0-9]{0,10}', s)
            analyses[-1].acqs[-1].background = [float(t) for t in backgrounds]

    return analyses


def CIDS_cleaner(analyses):
    '''function for cleaning up a parsed CIDS file and alerting to any corrupted analyses'''

    # Cleaning up the parsed file
    # because of the way the parser is written, an extra empty analysis is added to the end of the list
    if not analyses[-1].acqs:
        del analyses[-1]


    # checking for analyses with not enough acqs, alerting if there are some
    temp=8
    lowAcqs= [k for k in analyses if len(k.acqs)<temp]
    if not lowAcqs:
        print 'All analyses have at least %d acquistions' % temp

    else:
        print 'analyses with too few acqs are:'
        print '\n'.join(['Sample '+ str(k.num)+ ': '+ k.name +' has ' + str(len(k.acqs)) + ' acquisitions' for k in lowAcqs])

    # converting the voltages and backgrounds to arrays so that they can be easily used to do calculations
    for i in range(len(analyses)):
            for j in range(len(analyses[i].acqs)):
                # converting voltagted to arrays, and doing 3 sig figs to match CIDS Sheet
                analyses[i].acqs[j].voltSam=np.around(np.asarray(analyses[i].acqs[j].voltSam),3)
                analyses[i].acqs[j].voltRef=np.around(np.asarray(analyses[i].acqs[j].voltRef),3)
                analyses[i].acqs[j].background=np.asarray(analyses[i].acqs[j].background)
                # rounding d13C and d18O of each acq to 3 sig figs to match CIDS sheet
                (analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas) = np.around((analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas),3)


    print 'All analyses are cleaned, and voltages converted to arrays'

    return analyses

def CIDS_subtractBackground(analyses):
  '''Function to subtract the backgrounds off of all voltages'''

  for i in range(len(analyses)):
    for j in range(len(analyses[i].acqs)):
      background=analyses[i].acqs[j].background
      analyses[i].acqs[j].voltSam = analyses[i].acqs[j].voltSam-background
      analyses[i].acqs[j].voltRef = analyses[i].acqs[j].voltRef-background

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
    R_measured_analysis=(k.voltSam[:,1:6]/(np.tile(k.voltSam[:,0],(5,1)).T))
    R_measured_ref=(k.voltRef[:,1:6]/(np.tile(k.voltRef[:,0],(5,1)).T))
    delta_measured=np.zeros(np.shape(R_measured_analysis)) # Preallocating for size of delta array

    for l in range(len(R_measured_analysis)):
        delta_measured[l,:]=(R_measured_analysis[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
    # couches ratios in sample/std bracketing, put in delta notation

    delta_measured_mean=np.zeros((3,delta_measured.shape[1]))

    delta_measured_mean[0,:]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1,:]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

    R_calculated_analysis = (delta_measured_mean[0,:]/1000 + 1)*R_stoch[1,:]
    D_raw=(R_calculated_analysis/R_stoch[0,:] - 1)*1000


    k.D47_raw=D_raw[2]-D_raw[0]-D_raw[1]
    k.D48_raw=D_raw[3]-D_raw[0]-D_raw[1]
    k.d47=delta_measured_mean[0,2]
    k.d48=delta_measured_mean[0,3]

    k.d45=delta_measured_mean[0,0]
    k.d46=delta_measured_mean[0,1]

    acq=k

    return acq



def D47_calculations(analyses):
  '''Performs all the clumped isotope calculations on the raw voltages and bulk d13C and d18O values'''


  for i in range(len(analyses)):

    for j in range(len(analyses[i].acqs)):
    #   analyses[i].acqs[j]=D47_calculation(analyses[i].acqs[j])
      analyses[i].acqs[j]=carb_gas_oxygen_fractionation_acq(analyses[i].acqs[j])

    # CI_averages(analyses[i])

  return analyses


def CI_averages(analysis):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props1=['d45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O','d18O_min']

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev','d48','D48_raw','D48_stdev']

    acqsToUse = range(len(analysis.acqs))

    if analysis.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append([analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        analysis.acqs[i].d48, analysis.acqs[i].D47_raw,analysis.acqs[i].D48_raw,analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas, analysis.acqs[i].d18O_min])
        # values[i,:]=[analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        # analysis.acqs[i].d48, samle.acqs[i].D47_raw,analysis.acqs[i].D48_raw,analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas]

    if len(values) != 0:
        values = np.asarray(values)
        (analysis.d45, analysis.d46, analysis.d47, analysis.d48, analysis.D47_raw, analysis.D48_raw,
        analysis.d13C, analysis.d18O_gas, analysis.d18O_min) = values.mean(axis=0)

        # (analysis.d45, analysis.d46, analysis.d47, analysis.d48, analysis.D47_raw, analysis.D48_raw,
        # analysis.d13C, analysis.d18O) = values.mean(axis=0)

        (analysis.d45_stdev, analysis.d46_stdev, analysis.d47_stdev, analysis.d48_stdev, analysis.D47_stdev,
        analysis.D48_stdev, analysis.d13C_stdev,analysis.d18O_stdev, temp) = values.std(axis=0,ddof=1)

        analysis.D47_sterr=analysis.D47_stdev/np.sqrt(values.shape[0])

def FlatList_exporter(analyses,fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional flat list'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['User','date','Type','Sample ID','spec #\'s', 'acqs', 'd13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
    'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
    counter = 0
    if displayProgress:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, (len(item.acqs)-item.skipFirstAcq), item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev ])
            counter += 1
            if ((counter * 100)*100) % (len(analyses)*100) == 0:
                print(str((counter*100)/len(analyses)) + '% done')
    else:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, (len(item.acqs)-item.skipFirstAcq), item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev ])
    export.close()
    return

def CIDS_exporter(analyses, fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional CIDS file.
    Python CIDS files are also importable, to load all important sample info
    back into the program'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for analysis in analyses:
        wrt.writerow(['__NewSample__'])
        for acq in analysis.acqs:
            wrt.writerow(['__NewAcq__'])
            wrt.writerow(['__AcqVoltage__','',]+np.shape(acq.voltSam)[1]*['',]+acq.voltRef[0,:].tolist())
            for i in range(len(acq.voltSam)):
                wrt.writerow(['',''] + acq.voltSam[i,:].tolist() + acq.voltRef[i+1,:].tolist())
            wrt.writerow([])
            wrt.writerow(['','acq num','ref gas d13C','ref gas d18O','name','type','d13C (vpdb)',
            'd18O_gas (vsmow)', 'd45','d46','d47', 'd48', 'D47_raw', 'D48_raw'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.d13Cref, acq.d18Oref, analysis.name, analysis.type, acq.d13C,
            acq.d18O_gas, acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw])

        wrt.writerow([])
        wrt.writerow(['', 'd45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O_gas'])
        for acq in analysis.acqs:
            wrt.writerow(['',acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw, acq.d13C, acq.d18O_gas])
        wrt.writerow([])
        wrt.writerow(['','User','date','Type','Sample ID','spec #\'s','SkipFirstAcq','d13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
        'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
        wrt.writerow(['__SampleSummary__',analysis.user, analysis.date, analysis.type, analysis.name, analysis.num, analysis.skipFirstAcq, analysis.d13C, analysis.d13C_stdev,
        analysis.d18O_gas, analysis.d18O_min, analysis.d18O_stdev, analysis.d47, analysis.d47_stdev, analysis.D47_raw,
        analysis.D47_stdev, analysis.D47_sterr, analysis.d48, analysis.d48_stdev, analysis.D48_raw, analysis.D48_stdev])
        wrt.writerow(24*['---',])

    export.close()
    return

def CIDS_importer(filePath, displayProgress = False):
    '''Imports voltage and analysis data that were exported using the CIDS_exporter function'''
    fileImport = open(filePath,'rU')
    fileReader = csv.reader(fileImport, dialect='excel')
    analyses = []

    cycles_in_acqs = 7

    for line in fileReader:
        if '__NewSample__' in line:
            analyses.append(CI())
            continue
        if '__NewAcq__' in line:
            analyses[-1].acqs.append(ACQUISITION(0))
            continue
        if '__AcqVoltage__' in line:
            voltIndex = line.index('__AcqVoltage__')
            cycle = 0
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            analyses[-1].acqs[-1].voltRef.append(voltRefs)
            continue
        if cycle < (cycles_in_acqs):
            voltSams=[float(vs) for vs in line[voltIndex+2:voltIndex+8]]
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            analyses[-1].acqs[-1].voltSam.append(voltSams)
            analyses[-1].acqs[-1].voltRef.append(voltRefs)
            cycle += 1
            continue
        if '__AcqInfo__' in line:
            acqIndex = line.index('__AcqInfo__')
            analyses[-1].acqs[-1].acqNum = float(line[acqIndex + 1])
            analyses[-1].acqs[-1].d13Cref = float(line[acqIndex + 2])
            analyses[-1].acqs[-1].d18Oref = float(line[acqIndex + 3])
            analyses[-1].acqs[-1].d13C = float(line[acqIndex + 6])
            analyses[-1].acqs[-1].d18O_gas = float(line[acqIndex + 7])
            analyses[-1].acqs[-1].voltRef = np.asarray(analyses[-1].acqs[-1].voltRef)
            analyses[-1].acqs[-1].voltSam = np.asarray(analyses[-1].acqs[-1].voltSam)
            continue
        if '__SampleSummary__' in line:
            summaryIndex = line.index('__SampleSummary__')
            [analyses[-1].user, analyses[-1].date, analyses[-1].type, analyses[-1].name] = line[summaryIndex + 1:summaryIndex+5]
            analyses[-1].num = int(line[summaryIndex + 5])
            analyses[-1].skipFirstAcq = bool(line[summaryIndex + 6])

    fileImport.close()
    return analyses




def Get_gases(analyses):
    '''Finds which analyses are heated and equilibrated gases, and assigns them TCO2 values'''
    properNames = raw_input("Do all equilibrated gases have '25' in name? (y/n) ").lower()
    for item in analyses:
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
                        print('Skipping this analysis ')
                        break
                    else:
                        print('Not a valid entry, please try again')
    print('All heated and equilibrated gases have been found and labeled ')
    return

def Daeron_exporter(analyses, fileName):
    '''Exports analyses in a csv that is formatted for Matthieu Daeron's xlcor47 script'''

    export=open(fileName +'_daeron'+ '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    for item in analyses:
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
    R_measured_analysis=(acq.voltSam[:,1:6]/(np.tile(acq.voltSam[:,0],(5,1)).T))
    R_measured_ref=(acq.voltRef[:,1:6]/(np.tile(acq.voltRef[:,0],(5,1)).T))
    delta_measured=np.zeros(np.shape(R_measured_analysis)) # Preallocating for size of delta array

    for l in range(len(R_measured_analysis)):
        delta_measured[l,:]=(R_measured_analysis[l,:]/((R_measured_ref[l,:]+R_measured_ref[l+1,:])/2)-1)*1000
    # couches ratios in analysis/std bracketing, put in delta notation

    delta_measured_mean=np.zeros((3,delta_measured.shape[1]))

    delta_measured_mean[0,:]=np.mean(delta_measured, axis=0) # averaging for all cycles
    delta_measured_mean[1,:]=np.std(delta_measured, axis=0,ddof=1)  # standard deviation among all cycles
    delta_measured_mean[2,:]=delta_measured_mean[1,:]/np.sqrt(len(delta_measured)) # std error among all cycles

    R_calculated_analysis = (delta_measured_mean[0,:]/1000 + 1)*R_stoch[1,:]
    D_raw=(R_calculated_analysis/R_stoch[0,:] - 1)*1000


    d45=delta_measured_mean[0,0]
    d46=delta_measured_mean[0,1]
    d47=delta_measured_mean[0,2]
    d48=delta_measured_mean[0,3]
    D47_raw=D_raw[2]-D_raw[0]-D_raw[1]
    D48_raw=D_raw[3]-D_raw[0]-D_raw[1]


    calculatedCIValues = {'d45': d45, 'd46': d46, 'd47': d47, 'd48': d48, 'D47_raw': D47_raw, 'D48_raw': D48_raw}

    return calculatedCIValues[objName]

def CI_averages_valued(analysis, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48', 'd48_stdev','D48_raw','D48_stdev']

    acqsToUse = range(len(analysis.acqs))

    if analysis.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append([analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        analysis.acqs[i].d48, analysis.acqs[i].D47_raw,analysis.acqs[i].D48_raw,
        analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas, analysis.acqs[i].d18O_min])
        # values[i,:]=[analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        # analysis.acqs[i].d48, samle.acqs[i].D47_raw,analysis.acqs[i].D48_raw,analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas]

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

def CI_averages_valued_individual(analysis, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48','D48_raw','D48_stdev']

    valName = objName.replace('_stdev','')

    valName = valName.replace('_sterr','')

    if valName in ['D47', 'D48']:
        valName += '_raw'
    elif valName in ['d18O']:
        valName += '_gas'

    acqsToUse = range(len(analysis.acqs))

    if analysis.skipFirstAcq:
        del acqsToUse[0]

    values=[]#preallocating for value storage

    for i in acqsToUse:
        values.append(getattr(analysis.acqs[i],valName))
        # values.append([analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        # analysis.acqs[i].d48, analysis.acqs[i].D47_raw,analysis.acqs[i].D48_raw,
        # analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas, analysis.acqs[i].d18O_min])
        # values[i,:]=[analysis.acqs[i].d45,analysis.acqs[i].d46,analysis.acqs[i].d47,
        # analysis.acqs[i].d48, samle.acqs[i].D47_raw,analysis.acqs[i].D48_raw,analysis.acqs[i].d13C,analysis.acqs[i].d18O_gas]

    if len(values) != 0:
        values = np.asarray(values)
        if '_stdev' in objName:
            return values.std(axis=0,ddof=1)
        elif '_sterr' in objName:
            return values.std(axis=0,ddof=1)/np.sqrt(values.shape[0])
        else:
            return values.mean(axis=0)

def CI_comparer(analyses1, analyses2):
    '''compares two sets of CI data for equivalence in D47 calculations'''
    equalSoFar = True
    if not len(analyses1) == len(analyses2):
        print('Data sets are not the same length!')

    for i in range(len(analyses1)):
        for j in range(len(analyses1[i].acqs)):
            if not analyses1[i].acqs[j].D47_raw == analyses2[i].acqs[j].D47_raw:
                print('Acq num {1} from analysis number {0} from comparison does not agree in D47'.format(i,j))
                print('analyses1 D47 = {0:.3f}, analyses2 D47 = {1:.3f}'.format(analyses1[i].acqs[j].D47_raw, analyses2[i].acqs[j].D47_raw))
                equalSoFar = False
        if not analyses1[i].D47_stdev == analyses2[i].D47_stdev:
            print('Sample number {0} from data sets do not agree in D47_stdev'.format(i))
            print('1. D47_stdev = {0:.3f}, 2. D47_stdev = {0:.3f}'.format(analyses1[i].D47_stdev, analyses2[i].D48_stdev))

    if equalSoFar:
        print('Data sets are equivalent')

def Sample_type_checker(analyses):
    AllAnalysesHaveType = True
    for i in range(len(analyses)):
        if analyses[i].type not in ['std','sample','eg','hg']:
            AllAnalysesHaveType = False
            break

    return AllAnalysesHaveType

def Get_types_auto(analyses):
    '''Function to assign the correct type to every analysis automatically, given a few assumptions:
    stds are Carrara, NBS-19, or TV03, all gases have 'BOC' in name, and all egs have '25' in name)'''

    print('Automatically assigning analyses types ')
    choice = raw_input('(s)top process, see (n)aming guidelines, or hit any other key to continue ').lower()
    if choice == 's':
        return(analyses)
    elif choice == 'n':
        print('1. Standards must contain words "carrara", "TV03", or "NBS-19" ')
        print('2. 25 C equilibrated gases must contain "BOC" AND "25" ')
        print('3. 1000 C heated gases must contain "BOC" AND NOT "25" AND < 10 chars ')
        print('4. Analyses that already have a valid type are not modified ')
        return(analyses)
    else:
        for i in range(len(analyses)):
            if analyses[i].type in ['std', 'eg', 'hg', 'sample']:
                continue
            else:
                name = analyses[i].name.lower()
                if ('carrara' in name) or ('tv03' in name) or ('nbs-19' in name):
                    analyses[i].type = 'std';
                    continue
                elif ('boc' in name):
                    if ('25' in name):
                        analyses[i].type = 'eg'
                        continue
                    elif (len(name) < 10):
                        analyses[i].type = 'hg'
                        continue
                else:
                    analyses[i].type = 'sample'
    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Automatic assignment failed ')

    return(analyses)

def Get_types_manual(analyses):
    '''Function to assign the correct type to every analysis manually'''

    print('Manually assigning analyses types ')
    for i in range(len(analyses)):
        if analyses[i].type not in ['eg', 'hg', 'sample', 'std']:
            typeChoice = raw_input('Type for: ' + analyses[i].name + ' -> (e)g, (h)g, (s)ample, or s(t)d?').lower()
            if typeChoice == 'e':
                analyses[i].type = 'eg'
            elif typeChoice == 'h':
                analyses[i].type = 'hg'
            elif typeChoice == 't':
                analyses[i].type = 'std'
            else:
                analyses[i].type = 'sample'

    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Assignment of some types failed ')

    return(analyses)

def CI_CRF_data_corrector(analyses):
    '''Extended function to do all aspects of CRF data correction'''



def CI_CRF_corrector(analysis, objName):
    '''Function to apply the heated gas correction in the Caltech Ref Frame'''

    acid_digestion_correction = 0.081
    D47_hg_corrected = analysis.D47_raw - (analysis.d47*hg_slope + hg_intercept)
    D47_stretching = D47_hg_corrected *(-0.8453)/hg_intercept
    D47_CRF = D47_stretching + acid_digestion_correction

    return(D47_CRF)

def lsqfitma(X, Y):
    """
    Calculate a "MODEL-2" least squares fit.

    The line is fit by MINIMIZING the NORMAL deviates.

    The equation of the line is:     y = mx + b.

    This line is called the MAJOR AXIS.  All points are given EQUAL
    weight.  The units and range for X and Y must be the same.
    Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
    re-written from Kermack & Haldane (1950) Biometrika 37: 30-41;
    after a derivation by Pearson (1901) Phil. Mag. V2(6): 559-572.

    Data are input and output as follows:

    m, b, r, sm, sb = lsqfitma(X, Y)
    X    =    x data (vector)
    Y    =    y data (vector)
    m    =    slope
    b    =    y-intercept
    r    =    correlation coefficient
    sm   =    standard deviation of the slope
    sb   =    standard deviation of the y-intercept

    Note that the equation passes through the centroid:  (x-mean, y-mean)

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate sums and other re-used expressions.
    Sx = np.sum(X)
    Sy = np.sum(Y)
    xbar = Sx / n
    ybar = Sy / n
    U = X - xbar
    V = Y - ybar

    Suv = np.sum(U * V)
    Su2 = np.sum(U ** 2)
    Sv2 = np.sum(V ** 2)

    sigx = np.sqrt(Su2 / (n - 1))
    sigy = np.sqrt(Sv2 / (n - 1))

    # Calculate m, b, r, sm, and sb.
    m = (Sv2 - Su2 + np.sqrt(((Sv2 - Su2) ** 2) + (4 * Suv ** 2))) / (2 * Suv)
    b = ybar - m * xbar
    r = Suv / np.sqrt(Su2 * Sv2)

    sm = (m / r) * np.sqrt((1 - r ** 2) / n)
    sb1 = (sigy - sigx * m) ** 2
    sb2 = (2 * sigx * sigy) + ((xbar ** 2 * m * (1 + r)) / r ** 2)
    sb = np.sqrt((sb1 + ((1 - r) * m * sb2)) / n)

    return m, b, r, sm, sb


def lsqcubic(X, Y, sX, sY, tl=1e-6):
    """
    From:
    # teaching.py
#
# purpose:  Teaching module of ff_tools.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  09-Sep-2011
# modified: Sun 23 Jun 2013 04:30:45 PM BRT
#
# obs: Just some basic example function.

    Calculate a MODEL-2 least squares fit from weighted data.

    The line is fit by MINIMIZING the weighted residuals in both x & y.
    The equation of the line is:     y = mx + b,
    where m is determined by finding the roots to the cubic equation:

    m^3 + P * m^2 + Q * m + R = 0.

    Eqs for P, Q and R are from York (1966) Canad. J. Phys. 44: 1079-1086.

    Data are input and output as follows:
    m, b, r, sm, sb, xc, yc, ct = lsqcubic(X, Y, sX, sY, tl)
    X    =    x data (vector)
    Y    =    y data (vector)
    sX   =    uncertainty of x data (vector)
    sY   =    uncertainty of y data (vector)
    tl   =    test limit for difference between slope iterations

    m    =    slope
    b    =    y-intercept
    r    =    weighted correlation coefficient
    sm   =    standard deviation of the slope
    sb   =    standard deviation of the y-intercept
    xc   =    WEIGHTED mean of x values
    yc   =    WEIGHTED mean of y values
    ct   =    count: number of iterations

    Notes:  1.  (xc,yc) is the WEIGHTED centroid.
            2.  Iteration of slope continues until successive differences
                are less than the user-set limit "tl".  Smaller values of
                tl require more iterations to find the slope.
            3.  Suggested values of tl = 1e-4 to 1e-6.

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Find the number of data points and make one time calculations:
    n = len(X)
    wX = 1 / (sX ** 2)
    wY = 1 / (sY ** 2)

    # Set-up a few initial conditions:
    ct, ML = 0, 1

    # ESTIMATE the slope by calculating the major axis according
    # to Pearson's (1901) derivation, see: lsqfitma.

    MC = lsqfitma(X, Y)[0]

    test = np.abs((ML - MC) / ML)

    # Calculate the least-squares-cubic. Make iterative calculations until the
    # relative difference is less than the test conditions

    while test > tl:
        # Calculate sums and other re-used expressions:
        MC2 = MC ** 2
        W = (wX * wY) / ((MC2 * wY) + wX)
        W2 = W ** 2

        SW = np.sum(W)
        xc = (np.sum(W * X)) / SW
        yc = (np.sum(W * Y)) / SW

        U = X - xc
        V = Y - yc

        U2 = U ** 2
        V2 = V ** 2

        SW2U2wX = np.sum(W2 * U2 / wX)

        # Calculate coefficients for least-squares cubic:
        P = -2 * np.sum(W2 * U * V / wX) / SW2U2wX
        Q = (np.sum(W2 * V2 / wX) - np.sum(W * U2)) / SW2U2wX
        R = np.sum(W * U * V) / SW2U2wX
        # Find the roots to the least-squares cubic:
        LSC = [1, P, Q, R]
        MR = np.roots(LSC)

        # Find the root closest to the slope:
        DIF = np.abs(MR - MC)
        MinDif, Index = DIF.min(), DIF.argmin()

        ML = MC
        MC = MR[Index]
        test = np.abs((ML - MC) / ML)
        ct = ct + 1

    # Calculate m, b, r, sm, and sb.
    m = MC
    b = yc - m * xc
    r = np.sum(U * V) / np.sqrt(np.sum(U2) * np.sum(V2))
    sm2 = (1 / (n - 2)) * (np.sum(W * (((m * U) - V) ** 2)) / np.sum(W * U2))
    sm = np.sqrt(sm2)
    sb = np.sqrt(sm2 * (np.sum(W * (X ** 2)) / SW))

    return m, b, r, sm, sb, xc, yc, ct
