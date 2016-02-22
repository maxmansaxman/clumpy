'''Program to process raw isodat intensities for clumped CO2 measurements, and turn them into meaningful isotope ratios'''


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import CIDS_func
import os
from scipy.optimize import minimize


print('Welcome to the carbonate clumped isotope importer/exporter')
filePath = raw_input('autoCIDS file to process? ').rstrip()
filePath = os.path.abspath(filePath)
analyses = CIDS_func.CIDS_importer(filePath)
print('Analyses successfully imported')
print('Processing data in all relevant reference frames')
print('Processing data in caltech ref frame')
CIDS_func.CI_CRF_data_corrector(analyses)
print('Processing data in absolute ref frame, with the Daeron method')
CIDS_func.Daeron_data_processer(analyses,showFigures = True)

analyses_uncorrected = list(analyses)
all_D47_ARF_before = np.asarray([i.D47_ARF for i in analyses_uncorrected])
all_D47_ARF_errors_before = np.asarray([i.D47_error_all for i in analyses_uncorrected])


hgs = [i for i in analyses if (i.type == 'hg' and not i.D48_excess)]
# d47_hgs = np.asarray([i.d47 for i in hgs])
# D47_raw_hgs = np.asarray([i.D47_raw for i in hgs])
# d47_stdev_hgs = np.asarray([i.d47_stdev for i in hgs])
# D47_sterr_hgs = np.asarray([i.D47_sterr for i in hgs])

def hg_slope_calculator(pblGuess):
    ''' function to calculate the hg slope given a background correction for mass 47'''

    mass47slope = pblGuess[0]
    # mass47slope = pblGuess
    # mass47intercept = pblGuess[1]
    # mass47intercept = -3

    # Step 1, set new mass 47 pbl function
    CIDS_func.Set_mass_47_pbl(mass47slope)
    # Step 2, calculate new hg line with this
    hg_slope = abs(CIDS_func.CI_hg_slope_finder(hgs))
    print(str(mass47slope) + ','+ str(hg_slope))

    return hg_slope


pblGuessInit = [0]
# pblGuessInit = [-0.3]
# bnds = ((-1,0),(-5,1))
bnds = [(-1,0),]

res = minimize(hg_slope_calculator, pblGuessInit, method = 'SLSQP', bounds = bnds, options = {'ftol': 1e-9, 'eps': 1e-3, 'disp' : True, 'maxiter': 1000})
# res = minimize(hg_slope_calculator, pblGuessInit, method = 'nelder-mead', options = {'ftol': 1e-10, 'xtol':1e-10, 'disp' = True, 'maxfev': 1000})
# res = minimize(hg_slope_calculator, pblGuessInit, method = 'L-BFGS-B', bounds = bnds, options = {'ftol': 1e-2, 'gtol': 1e-12, 'disp' : True, 'maxiter': 1000})
# res = minimize(hg_slope_calculator, pblGuessInit, method = 'BFGS', options = {'ftol': 1e-12, 'xtol':1e-10, 'disp' : True, 'maxfev': 1000})



print('Re-processing data in caltech ref frame')
CIDS_func.CI_CRF_data_corrector(analyses)
print('Processing data in absolute ref frame, with the Daeron method')
CIDS_func.Daeron_data_processer(analyses, showFigures = True)

all_D47_ARF_after = np.asarray([i.D47_ARF for i in analyses])
all_D47_ARF_errors_after = np.asarray([i.D47_error_all for i in analyses])

plt.figure()
plt.plot(all_D47_ARF_before, (all_D47_ARF_after-all_D47_ARF_before),'o')
plt.xlabel(ur'$\mathrm{\Delta_{47, ARF}, \/before\/ PBL\/ correction \/(\u2030)}$')
plt.ylabel(ur'$\mathrm{(\Delta_{47,after} - \Delta_{47,before}) (\u2030)}$')
plt.show()

plt.figure()
plt.errorbar(all_D47_ARF_before, (all_D47_ARF_after-all_D47_ARF_before),yerr = all_D47_ARF_errors_before, fmt='o')
plt.xlabel(ur'$\mathrm{\Delta_{47, ARF}, \/before\/ PBL\/ correction \/(\u2030)}$')
plt.ylabel(ur'$\mathrm{(\Delta_{47,after} - \Delta_{47,before}) (\u2030)}$')
plt.show()


print('Exporting data to new sheets')
CIDS_func.ExportSequence(analyses, pbl = True)
