'''plots all permutations of enrichments'''

import numpy as np
import matplotlib.pyplot as plt

variables=np.shape(enrichments)[1]

varNames=['d45','d46','d47','d48','D47','D48']

for i in range(variables):
  for j in range(1,variables):
    if i==j:
      break
    else:
      plt.figure()
      plt.plot(enrichments[:,i],enrichments[:,j],'.')
      plt.xlabel(varNames[i])
      plt.ylabel(varNames[j])


for i in range(variables):
    plt.figure()
    plt.plot(deltas[:,i],excels[:,i],'.')
    plt.xlabel(varNames[i]+' python')
    plt.ylabel(varNames[i]+ ' excel')

for i in range(variables):
  for j in range(1,variables):
    if i==j:
      break
    else:
      plt.figure()
      plt.plot(deltas[:,i],deltas[:,j],'.')
      plt.xlabel(varNames[i])
      plt.ylabel(varNames[j])
