# copasi-utilities
This repository contains assorted files that are used to work with [COPASI](http://copasi.org). I wrote these files to aid in my own work and make them available for anyone else to use/modify.

A utility previously called *MEC.py* is now called *model_replicator* and has been moved to its [own repository](https://github.com/copasi/model_replicator).

## Descriptions
* **copasiplot2gnuplot.pl** - PERL script to convert a file saved from a COPASI plot window with several runs (eg from parameter scan) to a data file that [gnuplot](http://gnuplot.info/) can index the various runs (this basically just substitutes lines that have no numbers into two blank lines).

* **normalize_sequences.pl** - PERL script to convert a report file from COPASI parameter estimation or optimization to a regularly spaced file (in terms of function evaluations). More detailed information in the comments at the top of the file. 

* **progress_of_fit_stats.pl** - PERL script to calculate statistics on files containing several parameter estimation or optimization results which have been normalized to contain regular intervals (eg processed with *normalize_sequences.pl*). More detailed information in the comments at the top of the file.

* **run_pe_cps.sh** - BASH script to run a series of repeated COPASI parameter estimations (or optimizations) calculating their statistics. Useful for comparing algorithms. Uses PERL scripts *normalize_sequences.pl* and *progress_of_fit_stats.pl*. More details in the comments at the top of the file. 
* **poptbench** - BASH script to benchmark a parallel optimization algorithm against the serial version. Requires full path to serial and parallel versions of CopasiSE.

* **profilecopasi** - BASH script to profile CPU, memory, and I/O features of a CopasiSE run. Requires setting path to CopasiSE executable and uses the [audria](https://github.com/scaidermern/audria) package. Writes out a CSV file with the data of one or several runs.

* **model_report.py** - Python script to write a text summary of the model file (either in COPASI or SBML formats). Requires the *pandas* and *html2text* packages.
