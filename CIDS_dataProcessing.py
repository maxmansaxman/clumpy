'''Program to process raw isdat intensities, and turn them into meaningful isotope ratios'''

import csv

class CI:
  "A class for all the attributes of a single clumped isotope measurement"

  def __init__(self, name):
        self.name=name



class ACQ:
  "A class for all the attributes of a single clumped isotope acquision"










fileReader = csv.reader(open('~/Users/Max/Dropbox/CarbonateClumping/PressureBaselineCorrection/CIDS_04April_2014.csv', 'rU'), dialect='excel')
for row in fileReader:
    print row, ", "
