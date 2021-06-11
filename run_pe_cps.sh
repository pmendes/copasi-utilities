#! /usr/bin/bash
#
# usage: run_pe_cps.sh basefile iterations label {interval {final}}
#
# Written by Pedro Mendes <mendes@copasi.org>, May 2021. 
# Published under MIT License.
# https://raw.githubusercontent.com/pmendes/copasi-utilities/main/LICENSE
#
# This is designed to run a COPASI file with the parameter estimation
# task repeatedly for a fixed number of times and then calculate 
# statistics of its output in the form of mean, standard deviation, N, 
# min, and max values per each value of function evaluations. This 
# uses a couple of PERL scripts to normalize the data and to calculate
# the statistics; it will stop if it does not find those two scripts in
# the same folder (and marked as executable).
#
# Requires input of a COPASI file with parameter estimation task marked 
# for execution, an appropriately defined report file named 'pe_out.tsv'
# (which will later be renamed).

# Input arguments:
#  basefile:   basefile.cps is the name of the COPASI file to run; pass
#              only the filename without the extension '.cps'
#  iterations: how many times to run the parameter estimation
#  label:      label to describe this run, will be used to name files 
#              (together with basefile)
#  interval:   spacing of function evaluations for which we want 
#              statistics (defaults to 500)
#  final:      final value of function evaluations (defaults to 20001);
#              note that interval and final are related with the 
#              characteristics of the algorithm used (eg population and
#              generations, or swarm size and iterations, etc)

# ADJUSTABLE PARAMETERS

# cpsexec is the command line to call CopasiSE, adjust paths and options
cpsexec="/usr/local/copasi-dev/bin/CopasiSE -c /usr/share/copasi --nologo --report-file "

# default values for interval and final
interval=500
final=20001

# SANITY CHECKS

# are there sufficient arguments?
if [ $# -lt 3 ]
  then
   echo "usage: run_pe_cps.sh basefile iterations label {interval {final}}";
   exit 1
fi

# iterations must be a number
re='^[1-9][0-9]*$'
if ! [[ $2 =~ $re ]] ; then
   echo "ERROR: argument 2 (iterations) must be a positive integer";
   echo "usage: run_pe_cps.sh basefile iterations label {interval {final}}";
   exit 2
fi

# does the file exist?
if ! [ -f "$1.cps" ] ; then
   echo "ERROR: file $1.cps does not exist";
   exit 3
fi

# read optional arguments
if [ $# -gt 3 ]
  then
   interval=$4;
fi
if [ $# -gt 4 ]
  then
   final=$5;
fi

# do we have the other scripts here?
if ! [ -x normalize_sequences.pl ] ; then
   echo "ERROR: required script \"normalize_sequences.pl\" missing or not executable";
   exit 4
fi
if ! [ -x progress_of_fit_stats.pl ] ; then
   echo "ERROR: required script \"progress_of_fit_stats.pl\" missing or not executable";
   exit 4
fi

# DO THE WORK

# remove previous report file if one was there
rm "$1_$3.tsv" > /dev/null 2>&1

# run all the iterations
echo "Running $1.cps $2 iterations"
for (( i=0; i<$2; i++ ))
do
 $cpsexec "$1_$3.tsv" "$1.cps"
done

# normalize the output to a regular grid of fevals
echo "Normalizing sequences to $1_$3.norm.tsv"
./normalize_sequences.pl $1_$3.tsv $interval $final

# calculate the statistics
echo "Calculating statistics to $1_$3.norm.stats.tsv"
./progress_of_fit_stats.pl $1_$3.norm.tsv

# done
echo "Done."