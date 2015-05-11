''' Script to experiment with parsing .did files for clumped isotope parsing'''
import csv
import struct
import numpy as np



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

# Pulling out other auxiliary info
# Ref gas isotope composition
startRefGas=buff.find('CEvalDataSecStdTransferPart')

d13C_ref=struct.unpack('d',buff[startRefGas+203:startRefGas+203+8])
d18O_ref=struct.unpack('d',buff[startRefGas+423:startRefGas+423+8])
