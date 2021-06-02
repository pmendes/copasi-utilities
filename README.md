# copasi-utilities
This repository contains assorted files that are used to work with [COPASI](http://copasi.org). I wrote these files to aid in my own work and make them available for anyone else to use/modify.

## Descriptions
* **copasiplot2gnuplot.pl** - PERL script to convert a file saved from a COPASI plot window with several runs (eg from parameter scan) to a data file that [gnuplot](http://gnuplot.info/) can index the various runs (this basically just substitutes lines that have no numbers into two blank lines).
* **normalize_sequences.pl** - PERL script to convert a report file from COPASI parameter estimation or optimization to a regularly spaced file (in terms of function evaluations). More detailed information in the comments at the top of the file. 
