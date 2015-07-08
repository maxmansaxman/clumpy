# clumpy
Repository for all things clumped isotopes.
The purpose of this repository is to offer streamlined, open-source programs for processing and analyzing carbonate clumped isotope data.
This package is designed to standardize the data reduction processs, ensuring that all operations 
(e.g., projection into the absolute reference frame/CDES, standard correction, 48 excess handling) are performed in a unified way by all users in the community.
Programs are written in python 2.x for maximum readability, clarity, and flexibility in transferring between systems.

# Features
1. Isodat Parser: .did file processing: The heart of this package is a parser that reads Isodat's proprietary data filed (.did files) at the 
binary level to extract useful data such as voltages, sample ID, and the Isodat-calculated bulk d13C and d18O values in standard reference frames.
This parser can run in manual mode, where the user imports single acquisition files in turn, specifying to which samples they belong,
or in automatic mode, where a block of acquisitions in a directiory are automatically parsed and divided into their likely samples.

2. Excel Parser: Reads in raw acquisiton data and sample names from a CIDS sheet with a Caltech-style format.
Useful for importing old data that has been processed in a special way.

2. The Clumped Isotope class (CI): class with associated methods for handling data obtained with either of the above parsers. Built with the class are functions 
for properly calculating more complicated clumped isotope parameters, such as ∆47 and ∆48. 

# Coming soon
* Exports of summaries for each acquisition for data vetting and storage
* A unified program for projecting datasets into the absolute reference frame
* A GUI to make program operation more intuitive, without requiring command-line kung-fu



