''' Stores all the functions and classes used by various clumped isotope scripts'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import time
import xlcor47_modified
from scipy.optimize import root




CIT_Carrara_CRF = 0.352
TV03_CRF = 0.650
GC_AZ_CRF = 0.654

CIT_Carrara_ARF = 0.3115 + 0.092
NBS_19_ARF = CIT_Carrara_ARF
TV03_ARF = 0.635 + 0.092
GC_AZ_ARF = 0.710





global mass47PblSlope
mass47PblSlope = 0

global CarraraCorrection
CarraraCorrection = 0.0


class CI_VALUE(object):
    '''subclass defining how important isotopic ratios are calculated, and d18O_mineral, and brand (2010)-style d13C and d18O calcs'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if len(instance.voltRef)>6:
            if self.name in ['d45', 'd46', 'd47', 'd48', 'd49','D47_raw', 'D48_raw']:
                return np.around(D47_calculation_valued(instance, self.name),7)
            elif self.name in ['d13C_brand', 'd18O_brand', 'd18O_taylor', 'd13C_taylor']:
                return np.around(bulk_comp_brand_2010(instance, self.name),7)

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
            if self.name in ['d13C','d13C_stdev','d18O_gas','d18O_stdev', 'd45', 'd46', 'd47','d47_stdev',
            'D47_raw','D47_stdev', 'D47_sterr','d48','d48_stdev', 'd49','D48_raw','D48_stdev', 'd13C_brand',
            'd18O_brand', 'd18O_taylor', 'd13C_taylor']:
                return np.around(CI_averages_valued_individual(instance, self.name),7)
            elif self.name in ['d18O_min']:
                return np.around(carb_gas_oxygen_fractionation_acq(instance),7)

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
        if self.name in ['D47_CRF']:
            return np.around(CI_CRF_corrector(instance, self.name),5)
        elif self.name in ['D47_ARF_acid']:
            return np.around(CI_ARF_acid_corrector(instance, self.name),5)
        elif self.name in ['T_D47_ARF']:
            return np.around(CI_temp_calibrations(instance, self.name),2)
        elif self.name in ['D47_ARF_stdCorr']:
            return np.around(Carrara_carbonate_correction_ARF(instance, self.name), 5)

        else:
            raise ValueError('Sample D47_raw is out of range')


    def __set__(self, obj, value):
        raise AttributeError('Cannot change CI correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete CI corretion scheme')

class CI_UNIVERSAL_VALUE(object):
    '''subclass defining how heated gases are taken'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['hg_slope', 'hg_intercept']:
            return np.around(CI_hg_values(instance, self.name),5)

        else:
            raise ValueError('Not a valid value')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change individual hg value')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete individual hg value')

class CI_VOLTAGE(object):
    '''subclass defining how raw voltages are corrected'''

    def __init__(self, name):
        self.name = name

    def __get__(self,instance,cls):
        if self.name in ['voltSam', 'voltRef']:
            return np.around(CI_background_correction(instance, self.name),3)

        else:
            raise ValueError('Not a valid value')

    def __set__(self, obj, value):
        raise AttributeError('Cannot change individual voltage correction scheme')

    def __delete__(self, instance):
        raise AttributeError('Cannot delete individual voltage correction scheme')

class CI(object):
    "A class for all the attributes of a single clumped isotope measurement"

    def __init__(self):
        self.name=''
        self.acqs=[]
        self.date = ''
        self.type=''
        self.num=np.nan
        self.user=''
        self.skipFirstAcq = True
        self.TCO2 = np.nan
        self.D48_excess = False
        self.useBrand2010 = True
        self.rxnTemp = 90
        self.mineral = 'calcite'

        self.D47_ARF = np.nan
        self.D47_error_internal = np.nan
        self.D47_error_model = np.nan
        self.D47_error_all = np.nan

    d45 = CI_AVERAGE('d45')
    d46 = CI_AVERAGE('d46')
    d47 = CI_AVERAGE('d47')
    d48 = CI_AVERAGE('d48')
    d49 = CI_AVERAGE('d49')
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
    d13C_brand = CI_AVERAGE('d13C_brand')
    d18O_brand = CI_AVERAGE('d18O_brand')
    d18O_taylor = CI_AVERAGE('d18O_taylor')
    d13C_taylor = CI_AVERAGE('d13C_taylor')



    D47_CRF = CI_CORRECTED_VALUE('D47_CRF')
    hg_slope = CI_UNIVERSAL_VALUE('hg_slope')
    hg_intercept = CI_UNIVERSAL_VALUE('hg_intercept')

    D47_ARF_acid = CI_CORRECTED_VALUE('D47_ARF_acid')
    T_D47_ARF = CI_CORRECTED_VALUE('T_D47_ARF')
    D47_ARF_stdCorr = CI_CORRECTED_VALUE('D47_ARF_stdCorr')

class ACQUISITION(object):
    "A class for all the attributes of a single clumped isotope acquision"

    def __init__(self,acqNum):
        self.acqNum=acqNum
        self.voltSam_raw=[]
        self.voltRef_raw=[]
        self.background=[]
        self.d13C = 0
        self.d18O_gas = 0
        self.d13Cref = 0
        self.d18Oref = 0
        self.date=''
        self.time=''
        self.useBrand2010 = True
        self.pressureVals = [np.nan,np.nan,np.nan]
        # self.rxnTemp = 90
        # self.mineral = 'calcite'
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
    d49 = CI_VALUE('d49')
    d13C_brand = CI_VALUE('d13C_brand')
    d18O_brand = CI_VALUE('d18O_brand')
    d18O_taylor = CI_VALUE('d18O_taylor')
    d13C_taylor = CI_VALUE('d13C_taylor')

    # d18O_min = CI_VALUE('d18O_min')
    voltSam = CI_VOLTAGE('voltSam')
    voltRef = CI_VOLTAGE('voltRef')

def defliese_acid_equation(T_c):
    ''' calculates the will defliese phosphoric acid equation'''
    T_k = T_c + 273.15
    m = 0.022434
    m_err = 0.001490
    b = -0.25424
    b_err = 0.0168
    ln_alpha = m*10**6/T_k**2 + b

    return(ln_alpha)

def murray_acid_equation(T_c):
    ''' calculates the will defliese phosphoric acid equation'''
    T_k = T_c + 273.15
    m = 0.041757
    m_err = 0.011195
    b = -0.469746
    b_err = 0.034347
    ln_alpha = m*10**6/T_k**2 + b

    return(ln_alpha)

def carb_gas_oxygen_fractionation_acq(instance):
    '''calculates the d18O of a carbonate mineral from which the CO2 was digested'''
    # Done properly, this should be a function of the d18O of the gas, the reaction T,
    # and the cabonate phase of the analysis, but we're sticking to 90C calcite for now'''

    vsmow_18O=0.0020052
    vpdb_18O = 0 #don't know this right now
    rxnFrac = {'calcite_90': 1.00821, 'calcite_50': 1.0093, 'calcite_25': 1.01025,
    'dolomite_25': 1.01178, 'dolomite_50': 1.01038, 'dolomite_90': 1.009218, 'dolomite_100':1.00913, 'gas_25': 1.0, 'gas_90': 1.0}

    # calcite fractionations from swart et al., 1991
    # dolomite fractionations from Rosenbaum and sheppard, 1986
    # dolomite 90 is an interpotalion from 2nd-order polynomial fit of R&S data

    rxnKey = instance.mineral + '_' + str(instance.rxnTemp)
    if instance.useBrand2010:
        d18O_vpdb = (instance.d18O_gas-30.92)/1.03092
    else:
        d18O_vpdb = (instance.d18O_gas-30.86)/1.03086

    d18O_min = ((d18O_vpdb+1000)/rxnFrac[rxnKey])-1000
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
    voltRef_raw=[]
    voltSam_raw=[]
    #Slightly elegant method: pulling out based on spacing after 'CIntensityData' ascii key
    # TODO: make this more flexible
    # TODO: Catch errors if wrong voltage sequence found
    startPreVolt=start+2304+168 #observed location of 'pre' cycle on ref gas side
    voltRef_raw.append(struct.unpack('6d',buff[startPreVolt:(startPreVolt+6*8)]))
    for i in range(7):
        startRefVolt=start+52+i*164 #observed location of ref gas voltage cycles
        voltRef_raw.append(struct.unpack('6d',buff[startRefVolt:(startRefVolt+6*8)]))

        startSamVolt=start+1344+i*160 #observed location of sample gas voltage cycles
        voltSam_raw.append(struct.unpack('6d',buff[startSamVolt:(startSamVolt+6*8)]))

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

    # 3.4 background values, and Pressure values
    #find start of block with background values
    startBackground = buff.find('CISLScriptMessageData')
    stopBackground = buff.find('CMeasurmentErrors')
    #Note incorrect spelling of 'measurement' is intentional
    backgroundBlock = buff[startBackground+32:stopBackground].decode('utf-16', errors = 'ignore')
    # Pulling floats out of block
    backgroundVals = re.findall('[0-9]{1,20}\.[0-9]{1,10}', backgroundBlock)
    # Only want to take P values if they exist, bc a new acq with Brett's modified code
    if '100precent' in backgroundBlock:
        # pressure vals are first, and last two in block
        pressureVals = [float(backgroundVals[0]), float(backgroundVals[-2]), float(backgroundVals[-1])]
    else:
        pressureVals = [np.nan,np.nan,np.nan]
    # Background vals are other ones, but ignoring these for now



    # 3.5 Date and Time
    # Find the start of Ctime block
    startTime = buff.find('CTimeObject')+49
    # Pull out the time_t time based on startTime location, seconds since the epoch (Jan 1st, 1970), GMT
    time_t_time = struct.unpack('i',buff[startTime:startTime+4])[0]
    # time_str = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(time_t_time))
    time_str = time.strftime('%m/%d/%Y', time.localtime(time_t_time))


    return voltRef_raw, voltSam_raw, d13C_final, d18O_final, d13C_ref, d18O_ref, analysisName, firstAcq, time_str, pressureVals

def Isodat_File_Parser_CAF(fileName):
    '''Reads in a .caf file (Classical aquisition file), returns the raw voltages for
    ref gas, analysis gas, and the isodat-calculated d13C and d18O of the Acquisition'''

    f=open(fileName,'rb')
    try:
        buff = f.read()
    finally:
        f.close()

    #1. Getting raw voltage data
    #Searching for the start of the raw voltages
    start=buff.find('CDualInletRawData')
    keys=[]
    voltRef_raw=[]
    voltSam_raw=[]
    #Slightly elegant method: pulling out based on spacing after 'CIntensityData' ascii key
    # TODO: make this more flexible
    # TODO: Catch errors if wrong voltage sequence found
    startPreVolt=start+2304+168 #observed location of 'pre' cycle on ref gas side
    voltRef_raw.append(struct.unpack('6d',buff[startPreVolt:(startPreVolt+6*8)]))
    for i in range(7):
        startRefVolt=start+52+i*164 #observed location of ref gas voltage cycles
        voltRef_raw.append(struct.unpack('6d',buff[startRefVolt:(startRefVolt+6*8)]))

        startSamVolt=start+1344+i*160 #observed location of sample gas voltage cycles
        voltSam_raw.append(struct.unpack('6d',buff[startSamVolt:(startSamVolt+6*8)]))

    # test of where voltages are
    voltTest = []
    for i in range(len(buff)):
        thisStart = start + 2 + i
        theseVolts = np.array(struct.unpack('6d', buff[thisStart:(thisStart+6*8)]))
        if (theseVolts > 1e4).all() and (theseVolts < 1e6).all():
            print('Acceptable sequence at: {0}'.format(thisStart))
            voltTest.append(theseVolts)

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

    # 3.4 background values, and Pressure values
    #find start of block with background values
    startBackground = buff.find('CISLScriptMessageData')
    stopBackground = buff.find('CMeasurmentErrors')
    #Note incorrect spelling of 'measurement' is intentional
    backgroundBlock = buff[startBackground+32:stopBackground].decode('utf-16', errors = 'ignore')
    # Pulling floats out of block
    backgroundVals = re.findall('[0-9]{1,20}\.[0-9]{1,10}', backgroundBlock)
    # Only want to take P values if they exist, bc a new acq with Brett's modified code
    if '100precent' in backgroundBlock:
        # pressure vals are first, and last two in block
        pressureVals = [float(backgroundVals[0]), float(backgroundVals[-2]), float(backgroundVals[-1])]
    else:
        pressureVals = [np.nan,np.nan,np.nan]
    # Background vals are other ones, but ignoring these for now



    # 3.5 Date and Time
    # Find the start of Ctime block
    startTime = buff.find('CTimeObject')+49
    # Pull out the time_t time based on startTime location, seconds since the epoch (Jan 1st, 1970), GMT
    time_t_time = struct.unpack('i',buff[startTime:startTime+4])[0]
    # time_str = time.strftime('%m/%d/%Y %H:%M:%S', time.localtime(time_t_time))
    time_str = time.strftime('%m/%d/%Y', time.localtime(time_t_time))


    return voltRef_raw, voltSam_raw, d13C_final, d18O_final, d13C_ref, d18O_ref, analysisName, firstAcq, time_str, pressureVals

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
        try:
            acqNum=int(line[line.index('pre')-1])     # Not the same as acqIndex because acq are often skipped
        except(ValueError):
            acqNum = lastAcqNum + 1
        lastAcqNum = acqNum
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
          analyses[-1].acqs[-1].voltSam_raw.append(voltSams)
        except ValueError:
          pass
        voltRefs = [float(v2) for v2 in line[preIndex+7:preIndex+13]]
        analyses[-1].acqs[-1].voltRef_raw.append(voltRefs)
        cycle += 1

      if 'Date:' in line:
        analyses[-1].date= line[line.index('Date:')+1]
        # analyses[-1].acqs[-1].D47_excel=float(line[line.index('D47=')+1])

      elif 'Time:' in line:
        analyses[-1].acqs[-1].time=line[line.index('Time:')+1]
        # analyses[-1].acqs[-1].D48_excel=float(line[line.index('D48=')+1])

    #   elif 'Analyst' in line:
        # analyses[-1].acqs[-1].d45_excel=float(line[line.index('d45')+1])

    #   elif 'Type' in line:
        # analyses[-1].type=line[line.index('Type')+1].lower()

      elif 'Sample d13C (VPDB)' in line[0:10]:
        analyses[-1].acqs[-1].d13C=float(line[line.index('Sample d13C (VPDB)')+1])
        # analyses[-1].acqs[-1].d46_excel=float(line[line.index('d46')+1])

      elif 'Sample d18O (SMOW)' in line[0:10]:
        analyses[-1].acqs[-1].d18O_gas=float(line[line.index('Sample d18O (SMOW)')+1])
        # analyses[-1].acqs[-1].d47_excel=float(line[line.index('d47')+1])

      elif 'Ref Gas d13C (VPDB)' in line[0:10]:
        analyses[-1].acqs[-1].d13Cref=float(line[line.index('Ref Gas d13C (VPDB)')+1])
        # try:
        #   analyses[-1].acqs[-1].d48_excel=float(line[line.index('d48')+1])
        # except:
        #   print 'couldn\'t find d48 in acquistions %d' % acqIndex
        #   print line

      elif 'Ref Gas d18O (VSMOW)' in line[0:10]:
        analyses[-1].acqs[-1].d18Oref=float(line[line.index('Ref Gas d18O (VSMOW)')+1])

      else:
        for s in line:
          if "Background:" in s:
            backgrounds=re.findall('[0-9]{0,20}\.[0-9]{0,10}', s)
            analyses[-1].acqs[-1].background = [float(t) for t in backgrounds]

    return analyses

def CIDS_cleaner(analyses, checkForOutliers = False):
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
        print 'deleting those with just one acq...'
        for i in lowAcqs:
            if len(i.acqs) < 3:
                analyses.remove(i)



    # converting the voltages and backgrounds to arrays so that they can be easily used to do calculations
    for i in range(len(analyses)):
            for j in range(len(analyses[i].acqs)):
                # converting voltagted to arrays, and doing 3 sig figs to match CIDS Sheet
                analyses[i].acqs[j].voltSam_raw=np.around(np.asarray(analyses[i].acqs[j].voltSam_raw),3)
                analyses[i].acqs[j].voltRef_raw=np.around(np.asarray(analyses[i].acqs[j].voltRef_raw),3)
                analyses[i].acqs[j].background=np.asarray(analyses[i].acqs[j].background)
                # rounding d13C and d18O of each acq to 3 sig figs to match CIDS sheet
                (analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas) = np.around((analyses[i].acqs[j].d13C,analyses[i].acqs[j].d18O_gas),3)


    if checkForOutliers:
        Check_for_wrong_acqs(analyses)
    print 'All analyses are cleaned, and voltages converted to arrays'

    return analyses

def Check_for_wrong_acqs(analyses, showPlots = True):
    ''' Checks all analyses for individual acquisions that are clearly incorrect, using Grubbs' (1969) test '''
    print('Checking for outlier acqs...')
    sigLimit = 4
    alpha = 0.01
    # crit_vals_2sided = np.array([np.nan, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201])
    crit_vals_2sided = np.array([np.nan, 63.657, 9.965, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169, 3.106])
    for i in range(len(analyses)):
        # acqs_temp = analyses[i].acqs
        acqs_temp = analyses[i].acqs[analyses[i].skipFirstAcq:-1]
        d46_temp = np.asarray([j.d46 for j in acqs_temp])
        d46_res = np.abs(d46_temp - d46_temp.mean())
        G_2sided = d46_res.max()/d46_temp.std()

        if G_2sided > crit_vals_2sided[len(d46_temp)]:
            print('Sample {0} fails Grubbs\' outlier test at acq {1} with {2}-alpha significance'.format(i, np.argmax(d46_res), alpha))
            if showPlots:
                fig2, ax2 = plt.subplots(4)
                acq_temp = acqs_temp[np.argmax(d46_res)]
                voltLabels = ['mass44','mass45', 'mass46', 'mass47']
                ax2[0].set_title('Voltages for outlier acq number {0}, sample {1}'.format(acq_temp.acqNum, i))
                for i in range(4):
                    ax2[i].plot(np.arange(0,8),acq_temp.voltRef[:,i], '-', label = 'ref gas')
                    ax2[i].plot(np.arange(1,8), acq_temp.voltSam[:,i], '--', label = 'sample gas')

                    ax2[i].set_ylabel('{0} voltage (mV)'.format(voltLabels[i]))
                ax2[4].set_xlabel('cycle')
                ax2[0].legend(loc = 'best')

                fig0, ax0 = plt.subplots()
                # Plot d46
                ax0.plot(np.asarray(range(len(acqs_temp)))+analyses[i].skipFirstAcq,d46_temp, 'bo')
                # change axis color to blue
                ax0.set_ylabel(ur'$\delta 46$', {'color': 'b'})
                for tl in ax0.get_yticklabels():
                    tl.set_color('b')
                # Plot d47
                ax1 = ax0.twinx()
                ax1.plot(np.asarray(range(len(acqs_temp)))+analyses[i].skipFirstAcq,[j.d45 for j in acqs_temp], 'rd')
                # change axis color to red
                ax1.set_ylabel(ur'$\delta 45$', {'color': 'r'})
                for tl2 in ax1.get_yticklabels():
                    tl2.set_color('r')
                ax0.set_xlabel('cycle')

            delChoice = raw_input('Delete cycle {0}? (y/n) '.format(np.argmax(d46_res))).lower()
            if delChoice == 'y':
                del analyses[i].acqs[np.argmax(d46_res)+analyses[i].skipFirstAcq]
            plt.close('all')

    return

def FlatList_exporter(analyses,fileName, displayProgress = False):
    '''Exports a CSV file that is the same format as a traditional flat list'''

    export=open(fileName + '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['User','date','Type','Sample ID','spec #\'s', 'acqs', '100% bellow P (mbar)','d13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
    'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev', 'hg_slope', 'hg_intercept','D47_CRF', 'D47_ARF',
    'D47_ARF std error', 'mineral', 'rxnTemp', 'D47_ARF_acid', 'T(C)', 'D47_ARF_stdCorr', 'D48_excess', 'd13C_brand', 'd18O_brand'])
    counter = 0
    if displayProgress:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, (len(item.acqs)-item.skipFirstAcq), item.acqs[0].pressureVals[0],item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev, item.hg_slope, item.hg_intercept,
            item.D47_CRF, np.around(item.D47_ARF, 5), np.around(item.D47_error_all, 5), item.mineral, item.rxnTemp, item.D47_ARF_acid, item.T_D47_ARF, item.D47_ARF_stdCorr, item.D48_excess, item.d13C_brand, item.d18O_brand])
            counter += 1
            if ((counter * 100)*100) % (len(analyses)*100) == 0:
                print(str((counter*100)/len(analyses)) + '% done')
    else:
        for item in analyses:
            wrt.writerow([item.user, item.date, item.type, item.name, item.num, (len(item.acqs)-item.skipFirstAcq), item.acqs[0].pressureVals[0],item.d13C, item.d13C_stdev, item.d18O_gas, item.d18O_min,
            item.d18O_stdev,item.d47,item.d47_stdev,item.D47_raw, item.D47_stdev,item.D47_sterr,item.d48,item.d48_stdev,item.D48_raw,item.D48_stdev, item.hg_slope, item.hg_intercept,
            item.D47_CRF, np.around(item.D47_ARF, 5), np.around(item.D47_error_all,5), item.mineral, item.rxnTemp, item.D47_ARF_acid, item.T_D47_ARF, item.D47_ARF_stdCorr, item.D48_excess, item.d13C_brand, item.d18O_brand])
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
            wrt.writerow(['','acq num','date','ref gas d13C','ref gas d18O','name','type','d13C (vpdb)',
            'd18O_gas (vsmow)', 'd45','d46','d47', 'd48', 'D47_raw', 'D48_raw','100% bellow P','16V sample P','16V wg P'])
            wrt.writerow(['__AcqInfo__',acq.acqNum, acq.date, acq.d13Cref, acq.d18Oref, analysis.name, analysis.type, acq.d13C,
            acq.d18O_gas, acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw,acq.pressureVals[0],acq.pressureVals[1], acq.pressureVals[2]])

        wrt.writerow([])
        wrt.writerow(['', 'd45','d46','d47','d48','D47_raw','D48_raw','d13C','d18O_gas'])
        for acq in analysis.acqs:
            wrt.writerow(['',acq.d45, acq.d46, acq.d47, acq.d48, acq.D47_raw, acq.D48_raw, acq.d13C, acq.d18O_gas])
        wrt.writerow([])
        wrt.writerow(['','User','date','Type','Sample ID','spec #\'s','SkipFirstAcq','mineral','rxnTemp','d13C (vpdb)','d13C_stdev','d18O_gas (vsmow)','d18O_mineral (vpdb)',
        'd18O_stdev','d47','d47_stdev','D47 (v. Oz)','D47_stdev','D47_sterr','d48', 'd48_stdev','D48','D48_stdev'])
        wrt.writerow(['__SampleSummary__',analysis.user, analysis.date, analysis.type, analysis.name, analysis.num, analysis.skipFirstAcq, analysis.mineral, analysis.rxnTemp,analysis.d13C, analysis.d13C_stdev,
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
            analyses[-1].acqs[-1].voltRef_raw.append(voltRefs)
            continue
        if cycle < (cycles_in_acqs):
            voltSams=[float(vs) for vs in line[voltIndex+2:voltIndex+8]]
            voltRefs = [float(vr) for vr in line[voltIndex+8:voltIndex+14]]
            analyses[-1].acqs[-1].voltSam_raw.append(voltSams)
            analyses[-1].acqs[-1].voltRef_raw.append(voltRefs)
            cycle += 1
            continue
        if '__AcqInfo__' in line:
            acqIndex = line.index('__AcqInfo__')
            analyses[-1].acqs[-1].acqNum = float(line[acqIndex + 1])
            analyses[-1].acqs[-1].date = line[acqIndex + 2]
            analyses[-1].acqs[-1].d13Cref = float(line[acqIndex + 3])
            analyses[-1].acqs[-1].d18Oref = float(line[acqIndex + 4])
            analyses[-1].acqs[-1].d13C = float(line[acqIndex + 7])
            analyses[-1].acqs[-1].d18O_gas = float(line[acqIndex + 8])
            analyses[-1].acqs[-1].voltRef_raw = np.asarray(analyses[-1].acqs[-1].voltRef_raw)
            analyses[-1].acqs[-1].voltSam_raw = np.asarray(analyses[-1].acqs[-1].voltSam_raw)
            try:
                analyses[-1].acqs[-1].pressureVals[0] = float(line[acqIndex + 15])
                analyses[-1].acqs[-1].pressureVals[1] = float(line[acqIndex + 16])
                analyses[-1].acqs[-1].pressureVals[2] = float(line[acqIndex + 17])
            except(IndexError, ValueError):
                pass
            continue
        if '__SampleSummary__' in line:
            summaryIndex = line.index('__SampleSummary__')
            [analyses[-1].user, analyses[-1].date, analyses[-1].type, analyses[-1].name] = line[summaryIndex + 1:summaryIndex+5]
            analyses[-1].num = int(line[summaryIndex + 5])
            analyses[-1].skipFirstAcq = bool(line[summaryIndex + 6])
            if line[summaryIndex + 7] in ['calcite', 'dolomite', 'gas']:
                analyses[-1].mineral = line[summaryIndex + 7]
            if line[summaryIndex + 8] in ['90', '25', '50', '100']:
                analyses[-1].rxnTemp = int(line[summaryIndex + 8])



    fileImport.close()
    return analyses

def Get_carbonate_stds(analyses):
    '''Finds which analyses are carbonate standards, and assigns them accepted D47nominal values'''
    acid_correction_dict = {90: 0.092, 50: 0.040, 25: 0.0}
    for item in analyses:
        if item.type == 'std':
                if 'carrara' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = CIT_Carrara_ARF -acid_correction_dict[item.rxnTemp]
                elif 'tv03' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = TV03_ARF -acid_correction_dict[item.rxnTemp]
                elif 'nbs-19' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = NBS_19_ARF -acid_correction_dict[item.rxnTemp]
                elif 'gc-az' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = GC_AZ_ARF -acid_correction_dict[item.rxnTemp]
                elif 'gc_az' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = GC_AZ_ARF -acid_correction_dict[item.rxnTemp]
                elif 'tv04' in item.name.lower():
                    item.TCO2 = ''
                    item.D47nominal = TV03_ARF -acid_correction_dict[item.rxnTemp]
                else:
                    item.TCO2 = ''
                    item.D47nominal = ''

    print('All carbonate standards have been found and had accepted D47s assigned ')
    return

def Get_gases(analyses):
    '''Finds which analyses are heated and equilibrated gases, and assigns them TCO2 values'''
    # properNames = raw_input("Do all equilibrated gases have '25' in name? (y/n) ").lower()
    properNames = 'y'
    for item in analyses:
        if item.type in ['hg', 'eg']:
            if properNames == 'y':
                if '25' in item.name:
                    item.TCO2 = 25
                    item.D47nominal = ''
                    item.type = 'eg'
                else :
                    item.TCO2 = 1000
                    item.D47nominal = ''
                    item.type = 'hg'
            else :
                while True:
                    choice = raw_input('Is acq num: ' + str(item.num) + ' with name: ' + item.name + ' a (h)eated gas, an (e)quilibrated gas?, or (s)kip? ')
                    if choice.lower() == 'h':
                        item.D47nominal = ''
                        item.TCO2 = 1000
                        item.type = 'hg'
                        break
                    elif choice.lower() == 'e':
                        item.D47nominal = ''
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

def Daeron_exporter_crunch(analyses, fileName):
    '''Exports analyses in a csv that is formatted for Matthieu Daeron's online ClumpyCrunch 1.0 script'''

    export=open(fileName +'_daeron'+ '.csv','wb')
    wrt=csv.writer(export,dialect='excel')
    wrt.writerow(['type','ID', 'd45', 'd46', 'd47', 'd48', 'd49', 'sd47', 'D17O', 'd13Cwg_pdb', 'd18Owg_pdbco2', 'D47raw', 'TeqCO2', 'D47nominal'])
    for item in analyses:
        if item.type == 'sample':
            wrt.writerow([item.type, item.name, item.d45, item.d46, item.d47, item.d48, item.d49,
            item.d47_stdev/np.sqrt(len(item.acqs)-item.skipFirstAcq), 0.0, item.acqs[0].d13Cref, (item.acqs[0].d18Oref - 41.48693)/1.04148693,
            item.D47_raw])
        else :
            wrt.writerow([item.type, item.name, item.d45, item.d46, item.d47, item.d48, item.d49,
            item.d47_stdev/np.sqrt(len(item.acqs)-item.skipFirstAcq), 0.0, item.acqs[0].d13Cref, (item.acqs[0].d18Oref - 41.48693)/1.04148693,
            item.D47_raw, item.TCO2, item.D47nominal])
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
def bulk_comp_brand_2010(acq, objName):
    ''' Brand 2010 style calculation of bulk composition'''
    lambda_17 = 0.528
    K = 0.01022451
    vsmow_R18 = 0.0020052
    vsmow_R17 = 0.00038475
    vpdb_R13 = 0.011180

    # Derived quantities
    vpdb_R18 = vsmow_R18*1.03092*1.01025
    # vpdb_R17 = vsmow_R17*(1.03092*1.01025)**lambda_17

    # vpdb_R18 = 0.002088389
    vpdb_R17 = 0.00039310

    # Dummy holder eventually need to add in D17 of sample
    D17O_anomaly = 0
    K_general = np.exp(D17O_anomaly)*vpdb_R17*vpdb_R18**-lambda_17


    # convert wg d18O to vpdb, using the recommended value of Coplen et al., 1983 (once NBS-19 is corrected to +2.2 permille)
    d18Oref_vpdb = (acq.d18Oref-41.48693)/1.04148693

    # Calculate working gas ratios
    R13_ref=(acq.d13Cref/1000+1)*vpdb_R13
    R18_ref=(d18Oref_vpdb/1000+1)*vpdb_R18
    R17_ref=np.power((R18_ref/vpdb_R18),lambda_17)*vpdb_R17

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

    # Calculate ref gas stochastic ratios
    R45_ref_stoch = R13_ref + 2*R17_ref
    R46_ref_stoch=2*R13_ref*R17_ref+R17_ref*R17_ref+2*R18_ref
    R47_ref_stoch=2*R13_ref*R18_ref+2*R17_ref*R18_ref+R13_ref*R17_ref*R17_ref
    R48_ref_stoch=2*R17_ref*R18_ref*R13_ref+R18_ref*R18_ref
    R49_ref_stoch=R13_ref*R18_ref*R18_ref

    # assemble into array
    R_ref_stoch = np.array([R45_ref_stoch, R46_ref_stoch, R47_ref_stoch, R48_ref_stoch, R49_ref_stoch])
    # multiply by deltas to get 'true' ratios in sample gas. This assumes that wg is stochastic in composition
    R_measured_mean = (delta_measured_mean[0,:]/1000+1)*R_ref_stoch
    # # assemble extra arguments tuple to pass to solver
    # extraArgs = (R_measured_mean[0:3], lambda_17, K,)
    # # inital guess from true ratios of wg. Order is 13C, 17O, 18O
    # R_guess_init = [R13_ref, R17_ref, R18_ref]
    # # Minimize residuals for bulk solver
    # bulk_comp_result = root(bulk_comp_solver, R_guess_init, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})
    # r13C_brand, r17O_brand, r18O_brand = bulk_comp_result.x

    # TAYLOR POLYNOMIAL VERSION OF THE SAME
    # From Daeron et al., 2016, Appendix B
    # Calculate K specifically for the ref gas composition
    K = np.exp(D17O_anomaly) * vpdb_R17 * vpdb_R18 ** -lambda_17
    #taylor polynomials of the d46-d18O equations
    A_taylor = -3 * K**2 * (vpdb_R18**(2*lambda_17))
    B_taylor = 2 * K * R_measured_mean[0]*(vpdb_R18**lambda_17)
    C_taylor = 2 * vpdb_R18
    D_taylor = -R_measured_mean[1]

    a_taylor = A_taylor*lambda_17*(2*lambda_17-1) + B_taylor*lambda_17*(lambda_17-1)/2
    b_taylor = 2*A_taylor*lambda_17+B_taylor*lambda_17+C_taylor
    c_taylor = A_taylor + B_taylor + C_taylor + D_taylor
    # solve using quadratic eqn
    # in vpdb_co2 space
    d18O_taylor = 1000*(-b_taylor + (b_taylor**2-4*a_taylor*c_taylor)**0.5)/(2*a_taylor)
    R18_taylor = (d18O_taylor/1000+1)*vpdb_R18
    R17_taylor = K * R18_taylor**lambda_17
    R13_taylor = R_measured_mean[0] - 2*R17_taylor
    d13C_taylor = (R13_taylor/vpdb_R13-1)*1000

    d13C_vpdb_brand = d13C_taylor
    d18O_vpdb_brand = d18O_taylor
    d18O_vsmow_brand = d18O_vpdb_brand*1.04148693 + 41.48693


    #assemble dict
    bulk_comps = {'d13C_brand': d13C_vpdb_brand, 'd18O_brand': d18O_vsmow_brand, 'd18O_taylor': d18O_taylor, 'd13C_taylor': d13C_taylor}
    return(bulk_comps[objName])

def Brand_2010_bulk_comp_setter(analyses):
    'sets all acqs bulk comps to the Brand 2010-derived values'

    for i in analyses:
        i.useBrand2010 = True
        for j in i.acqs:
            j.d13C = j.d13C_brand
            j.d18O_gas = j.d18O_brand
            j.useBrand2010 = True

    return

def bulk_comp_solver(calcRatios, *extraArgs):
    ''' function to find roots of bulk compositions'''
    # unpacking extra args tuple
    measRatios, a, K = extraArgs
    # Santrock et al. (1985) eqn 10:
    R13_residual = calcRatios[0] + 2*calcRatios[1] - measRatios[0]
    # eqn 11:
    R18_residual = 2*calcRatios[2] + 2*calcRatios[0]*calcRatios[1] + calcRatios[1]**2 - measRatios[1]
    # eqn 19:
    R17_residual = K*calcRatios[2]**a - calcRatios[1]
    return(R13_residual, R18_residual, R17_residual)

def D47_calculation_valued(acq, objName):
    '''Performs all the clumped isotope calculations for a single acq'''

    vpdb_13C=0.0112372 # values copied from CIDS spreadsheet
    vsmow_18O=0.0020052
    vsmow_17O=0.0003799
    lambda_17=0.5164

    # if using Brand et al. (2010) correction
    if acq.useBrand2010:
        lambda_17 = 0.528
        K = 0.01022451
        vsmow_18O = 0.0020052
        vsmow_17O = 0.00038475
        vpdb_13C = 0.011180

        # Derived quantities
        vpdb_18O = vsmow_18O*1.03092*1.01025
        # vpdb_17O = vsmow_17O*(1.03092*1.01025)**lambda_17
        # vpdb_17O = K*vpdb_R18**a
        vpdb_17O = 0.0003931

        # convert wg d18O to vpdb, using the recommended value of Coplen et al., 1983 (once NBS-19 is corrected to +2.2 permille)
        d18Oref_vpdb = (acq.d18Oref-41.48693)/1.04148693

        # Calculate working gas ratios
        R13_ref=(acq.d13Cref/1000+1)*vpdb_13C
        R18_ref=(d18Oref_vpdb/1000+1)*vpdb_18O
        R17_ref=np.power((R18_ref/vpdb_18O),lambda_17)*vpdb_17O

        d18O_gas_vpdb = (acq.d18O_gas-41.48693)/1.04148693
        #
        R13_sa=(acq.d13C/1000+1)*vpdb_13C
        R18_sa=(d18O_gas_vpdb/1000+1)*vpdb_18O
        R17_sa=np.power((R18_sa/vpdb_18O),lambda_17)*vpdb_17O

    else:
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
    d49=delta_measured_mean[0,4]
    D47_raw=D_raw[2]-D_raw[0]-D_raw[1]
    D48_raw=D_raw[3]-D_raw[0]-D_raw[1]


    calculatedCIValues = {'d45': d45, 'd46': d46, 'd47': d47, 'd48': d48, 'd49': d49,'D47_raw': D47_raw, 'D48_raw': D48_raw}

    return calculatedCIValues[objName]

def CI_averages_valued_individual(analysis, objName):
    '''Computes the mean, std dev, and std error for every attribute useful for a CI measurement'''

    props2=['d13C','d13C_stdev','d18O_gas','d18O_min','d18O_stdev',
        'd47','d47_stdev','D47_raw','D47_stdev', 'D47_sterr','d48','D48_raw','D48_stdev', 'd13C_brand', 'd18O_brand']

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
        print('1. Standards must contain words "carrara", "TV03", or "NBS-19", "or GC-AZ" ')
        print('2. 25 C equilibrated gases must contain "BOC" AND "25" ')
        print('3. 1000 C heated gases must contain "BOC" AND NOT "25" AND < 15 chars ')
        print('4. Analyses that already have a valid type are not modified ')
        return(analyses)
    else:
        for i in range(len(analyses)):
            if analyses[i].type in ['std', 'eg', 'hg', 'sample']:
                continue
            else:
                name = analyses[i].name.lower()
                if ('carrara' in name) or ('tv03' in name) or ('nbs-19' in name) or ('gc-az' in name) or ('gc_az' in name):
                    analyses[i].type = 'std';
                    continue
                elif ('boc' in name):
                    if ('25' in name):
                        analyses[i].type = 'eg'
                        analyses[i].rxnTemp = 25
                        analyses[i].mineral = 'gas'
                        continue
                    elif (len(name) < 15):
                        analyses[i].type = 'hg'
                        analyses[i].rxnTemp = 25
                        analyses[i].mineral = 'gas'
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
                analyses[i].rxnTemp = 25
                analyses[i].mineral = 'gas'
            elif typeChoice == 'h':
                analyses[i].type = 'hg'
                analyses[i].rxnTemp = 25
                analyses[i].mineral = 'gas'
            elif typeChoice == 't':
                analyses[i].type = 'std'
            else:
                analyses[i].type = 'sample'

    if Sample_type_checker(analyses):
        print('Types successfully assigned ')
    else:
        print('Assignment of some types failed ')

    return(analyses)

def CI_CRF_data_corrector(analyses, showFigures = True):
    '''Extended function to do all aspects of CRF data correction'''

    # First, make a list of all hgs

    # Next, make a hg D48 line, using a York regression
    # hg_slope_48, hg_intercept_48, r_48, sm, sb = lsqfitma([i.d48 for i in hgs], [i.D48_raw for i in hgs])
    CI_48_excess_checker(analyses)
    # Make a new hg collection with any 48 excesses removed
    hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
    d47_hgs = np.asarray([i.d47 for i in hgs])
    D47_raw_hgs = np.asarray([i.D47_raw for i in hgs])
    d47_stdev_hgs = np.asarray([i.d47_stdev for i in hgs])
    D47_sterr_hgs = np.asarray([i.D47_sterr for i in hgs])
    # Now, make hg D47 line, using a York regression
    global hg_slope, hg_intercept
    hg_slope, hg_intercept, r_47, sm, sb, xc, yc, ct = lsqcubic(d47_hgs, D47_raw_hgs,d47_stdev_hgs, D47_sterr_hgs)

    D47_raw_hgs_model = d47_hgs*hg_slope+hg_intercept
    if showFigures:
        plt.figure(0)
        plt.figure(0).hold(True)
        plt.errorbar(d47_hgs, D47_raw_hgs, xerr = d47_stdev_hgs, yerr = D47_sterr_hgs, fmt = 'bo')
        plt.plot(d47_hgs, D47_raw_hgs_model, '-')
        plt.ylabel(ur'$\Delta_{47} \/ (\u2030)$')
        plt.xlabel(ur'$\delta^{47} \/ (\u2030)$')
        plt.show()

        # Calculation of the D47_CRF values are done automatically given its definition, based on the global hg slope and int
        # Plot the stds
        stds_Carrara = [i for i in analyses if (i.type == 'std' and 'carrara' in i.name.lower() and not i.D48_excess)]
        stds_TV03 = [i for i in analyses if (i.type == 'std' and 'tv03' in i.name.lower() and not i.D48_excess)]
        stds_GC_AZ = [i for i in analyses if (i.type == 'std' and 'gc-az' in i.name.lower() and not i.D48_excess)]
        plt.figure(1)
        plt.figure(1).hold(True)
        plt.errorbar([CIT_Carrara_CRF for i in range(len(stds_Carrara))],[j.D47_CRF-CIT_Carrara_CRF for j in stds_Carrara], yerr = [k.D47_sterr for k in stds_Carrara], fmt = 'o')
        plt.errorbar([TV03_CRF for i in range(len(stds_TV03))],[j.D47_CRF-TV03_CRF for j in stds_TV03], yerr = [k.D47_sterr for k in stds_TV03], fmt = 'o')
        plt.errorbar([GC_AZ_CRF for i in range(len(stds_GC_AZ))],[j.D47_CRF-GC_AZ_CRF for j in stds_GC_AZ], yerr = [k.D47_sterr for k in stds_GC_AZ], fmt = 'o')

        plt.xlim(0.3, 0.7)
        plt.xlabel(ur'$\mathrm{}\Delta_{47, CRF} \/ (\u2030)}$')
        plt.ylabel(ur'$\mathrm{\Delta_{47, measured}-\Delta_{47, expected} \/ (\u2030)}$')
        plt.show()

    return

def CI_hg_slope_finder(hgs):
    # hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
    d47_hgs = np.asarray([i.d47 for i in hgs])
    D47_raw_hgs = np.asarray([i.D47_raw for i in hgs])
    d47_stdev_hgs = np.asarray([i.d47_stdev for i in hgs])
    D47_sterr_hgs = np.asarray([i.D47_sterr for i in hgs])
    # Now, make hg D47 line, using a York regression
    global hg_slope, hg_intercept
    hg_slope, hg_intercept, r_47, sm, sb, xc, yc, ct = lsqcubic(d47_hgs, D47_raw_hgs,d47_stdev_hgs, D47_sterr_hgs)

    return hg_slope

def CI_48_excess_checker(analyses, showFigures = False):
    '''Checks for 48 excess using hg slope and int, based on a certain pre-specified tolerance'''
    D48_excess_tolerance = 4.0
    hgs = [i for i in analyses if i.type == 'hg']
    hg_slope_48, hg_intercept_48, r_48, sm_48, sb_48, xc_48, yc_, ct_ = lsqcubic(np.asarray([i.d48 for i in hgs]), np.asarray([i.D48_raw for i in hgs]),
    np.asarray([i.d48_stdev for i in hgs]), np.asarray([i.D48_stdev for i in hgs]))
    for i in analyses:
        D48_predicted = i.d48*hg_slope_48 + hg_intercept_48
        D48_excess_value = i.D48_raw - D48_predicted
        if np.abs(D48_excess_value) > D48_excess_tolerance:
            i.D48_excess = True
            print('D48 excess found for sample: ' + i.name + ', ' + str(i.num))


    if showFigures:
        plt.figure(1)
        plt.subplot(2,1,1)
        d48s = np.asarray([i.d48 for i in hgs])
        D48s_model = d48s * hg_slope_48 + hg_intercept_48
        # First, plotting the D48 line, TODO: include regression
        plt.figure(1).hold(True)
        plt.errorbar(np.asarray([i.d48 for i in hgs]), np.asarray([i.D48_raw for i in hgs]),
        xerr = np.asarray([i.d48_stdev for i in hgs]), yerr = np.asarray([i.D48_stdev for i in hgs]),
        fmt = 'bo')
        plt.plot(d48s, D48s_model,'r-')
        plt.xlabel(ur'$\delta^{48} \/ ( \u2030 $)')
        plt.ylabel(ur'$\Delta_{48} \/ ( \u2030 $)')

        # Plotting all samples in D47 vs d47 space,
        plt.subplot(2,1,2)
        plt.errorbar(np.asarray([i.d47 for i in analyses]), np.asarray([i.D47_raw for i in analyses]),
        xerr = np.asarray([i.d47_stdev for i in analyses]), yerr = np.asarray([i.D47_sterr for i in analyses]),
        fmt = 'o')
        plt.xlabel(ur'$\delta^{47} \/ ( \u2030 $)')
        plt.ylabel(ur'$\Delta_{47} \/ ( \u2030 $)')
        # plt.savefig('D48_line.pdf', format = 'pdf')
        plt.show()


    return

def Carrara_carbonate_correction_ARF(analysis, objName):
    ''' Function to apply a linear correction to carbonates based on CIT Carrara'''
    if analysis.type in ['hg', 'eg']:
        return(analysis.D47_ARF_acid)
    else:
        return(analysis.D47_ARF_acid - CarraraCorrection)

def Carrara_carbonate_correction_CRF(analysis, objName):
    ''' Function to apply a linear correction to carbonates based on CIT Carrara'''
    if analysis.type in ['hg', 'eg']:
        return(analysis.D47_CRF)
    else:
        return(analysis.D47_CRF - CarraraCorrection_CRF)

def CI_CRF_corrector(analysis, objName):
    '''Function to apply the heated gas correction in the Caltech Ref Frame'''

    acid_correction_dict = {'90': 0.081, '50': 0.040, '25': 0.0}
    acid_digestion_correction = acid_correction_dict[str(analysis.rxnTemp)]
    try:
        D47_hg_corrected = analysis.D47_raw - (analysis.d47*hg_slope + hg_intercept)
        D47_stretching = D47_hg_corrected *(-0.8453)/hg_intercept
        D47_CRF = D47_stretching + acid_digestion_correction
    except NameError:
        return(np.NaN)

    return(D47_CRF)

def CI_ARF_acid_corrector(analysis, objName):
    '''Function to apply the acid correction to an arf value'''
    acid_correction_dict = {'90': 0.092, '50': 0.040, '25': 0.0}
    acid_digestion_correction = acid_correction_dict[str(analysis.rxnTemp)]
    try:
        D47_ARF_acid = analysis.D47_ARF + acid_digestion_correction
    except NameError:
        return(np.NaN)

    return(D47_ARF_acid)

def CI_temp_calibrations(analysis, objName):
    '''Function to apply a temperature calibration'''
    #for now, just doing the boniface + Henkes ARF calibration
    # T = np.sqrt(0.0421e6/(D47_ARF_acid - 0.211)) - 273.15

    # New ARF function, based on Stolper 2015, Bonifacie 2011, Guo 2009, and Ghosh 2006 data
    # projected into ARF using absolute slope
    a = 0.00108331
    b = 0.0285392
    c = 0.258652-analysis.D47_ARF_acid
    TKe6 = (-b +np.sqrt(b**2-4*a*c))/(2*a)
    T_C = np.sqrt(1e6/TKe6)-273.15
    return(T_C)

def CI_D47_to_temp_ARF(D47_ARF_acid):
    '''Function to apply a temperature calibration'''
    #for now, just doing the boniface + Henkes ARF calibration
    # T = np.sqrt(0.0421e6/(D47_ARF_acid - 0.211)) - 273.15

    # NEW Boniface 2016 calibration
    TKe6 = (D47_ARF_acid-0.1262)/0.0422

    # New ARF function, based on Stolper 2015, Bonifacie 2011, Guo 2009, and Ghosh 2006 data
    # projected into ARF using absolute slope
    # a = 0.00108331
    # b = 0.0285392
    # c = 0.258652-D47_ARF_acid
    # a = 0.001083
    # b = 0.02854
    # c = 0.25865-D47_ARF_acid
    # TKe6 = (-b +np.sqrt(b**2-4*a*c))/(2*a)
    if TKe6 > 0:
        T_C = np.sqrt(1e6/TKe6)-273.15
    else:
        T_C = np.nan
    return(T_C)

def CI_temp_to_D47_ARF(T_C):
    '''Function to apply a temperature calibration'''
    #for now, just doing the boniface + Henkes ARF calibration
    # T = np.sqrt(0.0421e6/(D47_ARF_acid - 0.211)) - 273.15

    # New ARF function, based on Stolper 2015, Bonifacie 2011, Guo 2009, and Ghosh 2006 data
    # projected into ARF using absolute slope
    # a = 0.00108331
    # b = 0.0285392
    # c = 0.258652
    TKe6 = 1e6/(T_C+273.15)**2
    # D47_ARF_acid = a*TKe6**2+b*TKe6+c
    D47_ARF_acid= 0.0422*TKe6 + 0.1262

    return(D47_ARF_acid)

def CI_D47_to_temp_CRF(D47_CRF_acid):
    '''Function to apply a temperature calibration'''
    # stolper 2015 CRF calibration
    a = 1.006e-3
    b = 2.620e-2
    c = 0.2185-D47_CRF_acid
    TKe6 = (-b +np.sqrt(b**2-4*a*c))/(2*a)
    if TKe6 > 0:
        T_C = np.sqrt(1e6/TKe6)-273.15
    else:
        T_C = np.nan
    return(T_C)

def CI_water_calibration(d18O_min, T_C):
    '''Function to determine d18O of the water from which calcite crystallized
    Using the Kim and O'Neil 1997 calibration line
    Note, only valid from 10-40 C'''

    #for now, just doing the boniface + Henkes ARF calibration
    d18O_min_vsmow = d18O_min*1.03092 + 30.92
    R18_vsmow = 0.0020052
    R18_min = (d18O_min_vsmow/1000+1)*R18_vsmow
    T_K = T_C + 273.15
    # The Kim and Oneil 1997 equation:
    alpha_calcite_H2O = np.exp(((18.03*(1e3/T_K)-32.42)/1000))
    R18_H2O = R18_min/alpha_calcite_H2O
    d18O_H2O_vsmow = (R18_H2O/R18_vsmow-1)*1000

    return(d18O_H2O_vsmow)


def CI_carb_water_calibration(d18O_min, T_C, mineral):
    ''' function to determine fluid composition, in dolomite or calcite'''
    d18O_min_vsmow = d18O_min*1.03092 + 30.92
    R18_vsmow = 0.0020052
    R18_min = (d18O_min_vsmow/1000+1)*R18_vsmow
    T_K = T_C + 273.15
    if mineral == 'dolomite':
        # Horita 2014 calibration
        alpha_carb_H2O = np.exp((3.140*(1e6/T_K**2)-3.14)/1000)
    else:
        #O'neil 1969, Friedman and o'neil 1977 calibration (via horita)
        alpha_carb_H2O = np.exp((2.789*(1e6/T_K**2)-2.89)/1000)

    R18_H2O = R18_min/alpha_carb_H2O
    d18O_H2O_vsmow = (R18_H2O/R18_vsmow-1)*1000
    return(d18O_H2O_vsmow)

def CI_hg_values(instance, objName):
    '''Placeholder for hg slope and int when used in CI classes'''
    try:
        values = {'hg_slope': hg_slope, 'hg_intercept': hg_intercept}
        return(values[objName])
    except(NameError, AttributeError):
        return(np.NaN)

def Daeron_data_creator(analyses, useCarbStandards = False):
    '''Creates a list of dictionaries in the format needed for M. Daeron's data
    processing script '''
    acid_correction_dict = {90: 0.092, 50: 0.040, 25: 0.0}
    daeronData = []
    if useCarbStandards:
        for i in analyses:
            daeronData.append({'label': i.name, 'd45': i.d45, 'd46':i.d46,'d47': i.d47, 'd48': i.d48, 'd49': i.d49,
            'sd47': i.d47_stdev/np.sqrt(len(i.acqs)-i.skipFirstAcq), 'D17O': 0.0, 'd13Cwg_pdb': i.acqs[0].d13Cref, 'd18Owg_pdbco2': (i.acqs[0].d18Oref - 41.48693)/1.04148693,
            'D47raw': i.D47_raw, 'D47_raw_sterr': i.D47_sterr} )
            if i.type == 'hg':
                daeronData[-1]['TCO2eq'] = 1000.0
                daeronData[-1]['D47nominal'] = xlcor47_modified.CO2eqD47(1000.0)
            elif i.type == 'eg':
                daeronData[-1]['TCO2eq'] = 25.0
                daeronData[-1]['D47nominal'] = xlcor47_modified.CO2eqD47(25.0)
            elif i.type == 'std':
                if 'carrara' in i.name.lower():
                    daeronData[-1]['D47nominal'] = CIT_Carrara_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan
                elif 'tv03' in i.name.lower():
                    daeronData[-1]['D47nominal'] = TV03_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan
                elif 'nbs-19' in i.name.lower():
                    daeronData[-1]['D47nominal'] = NBS_19_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan
                elif 'gc-az' in i.name.lower():
                    daeronData[-1]['D47nominal'] = GC_AZ_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan
                elif 'gc_az' in i.name.lower():
                    daeronData[-1]['D47nominal'] = GC_AZ_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan
                elif 'tv04' in i.name.lower():
                    daeronData[-1]['D47nominal'] = TV03_ARF - acid_correction_dict[i.rxnTemp]
                    daeronData[-1]['TCO2eq'] = np.nan

    else:
        for i in analyses:
            daeronData.append({'label': i.name, 'd45': i.d45, 'd46':i.d46,'d47': i.d47, 'd48': i.d48, 'd49': i.d49,
            'sd47': i.d47_stdev/np.sqrt(len(i.acqs)-i.skipFirstAcq), 'D17O': 0.0, 'd13Cwg_pdb': i.acqs[0].d13Cref, 'd18Owg_pdbco2': (i.acqs[0].d18Oref - 41.48693)/1.04148693,
            'D47raw': i.D47_raw, 'D47_raw_sterr': i.D47_sterr} )
            if i.type == 'hg':
                daeronData[-1]['TCO2eq'] = 1000.0
                daeronData[-1]['D47nominal'] = xlcor47_modified.CO2eqD47(1000.0)
            elif i.type == 'eg':
                daeronData[-1]['TCO2eq'] = 25.0
                daeronData[-1]['D47nominal'] = xlcor47_modified.CO2eqD47(25.0)
    return(daeronData)

def Daeron_data_processer(analyses, showFigures = False):
    '''Performs Daeron-style correction to put clumped isotope data into the ARF'''
    useCarbStandards = False
    askCarbStandards = raw_input('Use carbonate standards for ARF correction? (y/n) ').lower()
    if askCarbStandards == 'y':
        useCarbStandards = True

    daeronData = Daeron_data_creator(analyses, useCarbStandards)

    global daeronBestFitParams, CorrelationMatrix
    daeronBestFitParams, CorrelationMatrix = xlcor47_modified.process_data(daeronData)

    for i in range(len(daeronData)):
        analyses[i].D47_ARF = daeronData[i]['corD47']
        analyses[i].D47_error_internal = daeronData[i]['scorD47_internal']
        analyses[i].D47_error_model = daeronData[i]['scorD47_model']
        analyses[i].D47_error_all = daeronData[i]['scorD47_all']

    if showFigures:
        xlcor47_modified.plot_data( daeronData, daeronBestFitParams, CorrelationMatrix ,'filename')

    # Calculating CIT Carrara offset
    stds_Carrara = [i for i in analyses if (i.type == 'std' and 'carrara' in i.name.lower() and not i.D48_excess)]
    global CarraraCorrection
    CarraraCorrection = np.mean([i.D47_ARF_acid for i in stds_Carrara]) - CIT_Carrara_ARF
    print('Carrara Correction is: '+ str(CarraraCorrection))

    return
def ExportSequence(analyses, pbl = False):
    '''Most common export sequence'''
    print('Exporting to temporary CIDS and FlatList files ')
    print('Exporting full acqs to a CIDS sheet...')
    exportNameCIDS = 'autoCIDS_Export_'
    if pbl:
        exportNameCIDS += '_pblCorr'
    CIDS_exporter(analyses, exportNameCIDS)
    print('Exporting analyses to a flatlist...')
    exportNameFlatlist = 'autoFlatListExport'
    if pbl:
        exportNameFlatlist += '_pblCorr'
    FlatList_exporter(analyses,exportNameFlatlist)
    print('Analyses successfully exported')
    doDaeron = raw_input('Export analyses for a Daeron-style ARF reduction (y/n)? ')
    if doDaeron.lower() == 'y':
        exportNameDaeron = 'autoDaeronExport'
        Get_gases(analyses)
        Get_carbonate_stds(analyses)
        Daeron_exporter_crunch(analyses,exportNameDaeron)

    return
def ExportSequence_named(analyses, fileName, pbl = False):
    '''Most common export sequence'''
    # split filename to get part before file type
    fileName = fileName.rpartition('.')[0]
    print('Exporting to temporary CIDS and FlatList files ')
    print('Exporting full acqs to a CIDS sheet...')
    exportNameCIDS = 'autoCIDS_Export_' + fileName
    if pbl:
        exportNameCIDS += '_pblCorr'
    CIDS_exporter(analyses, exportNameCIDS)
    print('Exporting analyses to a flatlist...')
    exportNameFlatlist = 'autoFlatListExport' + fileName
    if pbl:
        exportNameFlatlist += '_pblCorr'
    FlatList_exporter(analyses,exportNameFlatlist)
    print('Exporting daeron-formatted file for rambaldi.pythonanywhere.com...')
    exportNameDaeron = 'autoDaeronExport' + fileName
    Get_gases(analyses)
    Get_carbonate_stds(analyses)
    Daeron_exporter_crunch(analyses,exportNameDaeron)
    print('Analyses successfully exported')

    return

def plot_T_scale(axis, in_ARF = False):
    ''' Plots a secondary y-axis with temperature scale to match the D47 on ARF '''
    # 1. Get D47 labels
    tlabels = axis.get_yticklabels(which = 'both')
    # 2. For some reason, have to save fig first
    fig = plt.gcf()
    fig.savefig('temp.pdf')
    y2_temps = []
    for i in range(len(tlabels)):
        tl = tlabels[i]
        D47_str = tl.get_text()
        if len(D47_str) > 0:
            if in_ARF:
                y2_temps.append('{0:.0f}'.format(np.around(CI_D47_to_temp_ARF(float(D47_str)),0)))
            else:
                y2_temps.append('{0:.0f}'.format(np.around(CI_D47_to_temp_CRF(float(D47_str)),0)))
        else:
            y2_temps.append('')
    axis2 = axis.twinx()
    axis2.grid(False)
    axis2.set_ylabel(ur'Temperature ($^{\circ}$C)')
    axis2.set_yticklabels(y2_temps)
    return

def CI_background_correction(instance, objName):
    # voltSamTemp = np.copy(instance.voltSam_raw)
    # voltRefTemp = np.copy(instance.voltRef_raw)

    slopeArray = np.array([0, 0 , 0, mass47PblSlope, 0, 0])
    # interceptArray  = np.array([0, 0, 0, mass47PblIntercept, 0, 0])

    # mass 47 correction
    try:
        # print('Correcting voltages now')
        # print('using this mass47 slope: '+ str(mass47PblSlope))
        voltSamTemp = instance.voltSam_raw - np.transpose(np.tile(instance.voltSam_raw[:,0],(len(slopeArray),1)))*slopeArray
        voltRefTemp = instance.voltRef_raw - np.transpose(np.tile(instance.voltRef_raw[:,0],(len(slopeArray),1)))*slopeArray
    except(IndexError):
        return 0

    voltages = {'voltSam': voltSamTemp, 'voltRef': voltRefTemp}

    return(voltages[objName])

def Set_mass_47_pbl(pblSlope):
    global mass47PblSlope
    mass47PblSlope = pblSlope
    # mass47PblIntercept = pblIntercept

def D47_mixing(x_1, d13C_1, d18O_1, D47_1, x_2, d13C_2, d18O_2, D47_2):
    ''' Mixing between two endmembers, in concentrations, assuming d13C in vpdb and d18O in vsmow'''

    vpdb_13C=0.01118 # values copied from CIDS spreadsheet
    vsmow_18O=0.0020052
    vsmow_17O=0.00038475
    lambda_17=0.528

    R13_1=(d13C_1/1000+1)*vpdb_13C
    R18_1=(d18O_1/1000+1)*vsmow_18O
    R17_1=np.power((R18_1/vsmow_18O),lambda_17)*vsmow_17O

    R13_2=(d13C_2/1000+1)*vpdb_13C
    R18_2=(d18O_2/1000+1)*vsmow_18O
    R17_2=np.power((R18_2/vsmow_18O),lambda_17)*vsmow_17O


    # calculating stochastic ratios
    # to keep things organized and avoid repetetive code lines,
    # organizing calculations in arrays where item 1 is sample gas, item 2 is ref gas
    R13=np.array([R13_1, R13_2])
    R17=np.array([R17_1, R17_2])
    R18=np.array([R18_1, R18_2])

    x = np.array([x_1, x_2])


    R45_stoch=R13+2*R17
    R46_stoch=2*R13*R17+R17*R17+2*R18
    R47_stoch=2*R13*R18+2*R17*R18+R13*R17*R17
    R48_stoch=2*R17*R18*R13+R18*R18
    R49_stoch=R13*R18*R18

    R47 = (np.array([D47_1, D47_2])/1000 + 1)*R47_stoch

    c12 = 1/(1+R13)
    c13 = c12*R13
    c16 = 1/(1+R17+R18)
    c17 = R17*c16
    c18 = R18*c16

    c44 = c12*c16*c16
    c47 = R47*c44

    c12_m = np.sum(c12*x)
    c13_m = np.sum(c13*x)
    c16_m = np.sum(c16*x)
    c17_m = np.sum(c17*x)
    c18_m = np.sum(c18*x)

    c44_m = c12_m*c16_m*c16_m
    c47_m = np.sum(c47*x)

    R13_m = c13_m/c12_m
    R17_m = c17_m/c16_m
    R18_m = c18_m/c16_m

    R47_m_stoch = 2*R13_m*R18_m+2*R17_m*R18_m+R13_m*R17_m*R17_m

    D47_m = (c47_m/c44_m/R47_m_stoch-1)*1000
    d13C_m = (R13_m/vpdb_13C-1)*1000
    d18O_m = (R18_m/vsmow_18O-1)*1000

    return(d13C_m,d18O_m,D47_m)

def lsqfitma(X, Y):
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
