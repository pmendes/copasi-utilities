#!/usr/bin/perl -w
#
# progress_of_fit_stats.pl inputfile
#
# Written by Pedro Mendes <mendes@copasi.org>, May 2021. 
# Published under MIT License.
# https://raw.githubusercontent.com/pmendes/copasi-utilities/main/LICENSE
#
# Calculates statistics for objective function values in normalized COPASI 
# progress of fit or progress of optimization data files. This is intended
# for files that have been already normalized to have fixed interval size,
# as can be done with normalize_sequences.pl. This is not intended to
# process the files directly produced by COPASI (as they don't have fixed
# intervals).
#
# Parameters:
#  inputfile: file to be processed
#
# The file is expected to contain several sequences with the same spacing.
# Statistics are calculated for the values of objective function value 
# (second column of input file) for each value of function evaluations
# (first column of input file). The statistics calculated are:
#  Mean value
#  Standard deviation
#  Number of samples
#  Minimum value
#  Maximum value
#
# Output file is named based on the inputfile; examples:
#  input: filename           output: filename.stats.tsv
#  input: basename.tsv       output: basename.stats.tsv   
#  input: basename.ext       output: basename.stats.tsv   
#  input: basename.norm.tsv  output: basename.norm.stats.tsv   

use POSIX;

if( $#ARGV < 0 ) { die "usage: progress_of_fit_stats.pl inputfile\n"; }
$in = $ARGV[0];

# extract filename without extension if possible
$basename = $in;
$basename =~ s/(.+)\..*$/$1/;

$fendvals = "$basename.stats.tsv";

open( IFILE, "$in" )  || die "can't read $in";

# file for end values
open( O1FILE, ">$fendvals" )  || die "can't open $fendvals";

# counter for sequence number
$n = -1;
# max counter of all sequences
$top = 0;
# counter for sequence replicates
$replicate = -1;
# feval read in last line
$lastline = -1;

# arrays to calculate the statistics
# the number of iterations for each position
@iter = ();
# the mean ssq values
@mean = ();
# the count of values for this position
@count = ();
# the cumulative square distance from the mean of ssq
@m2 = ();
# the minimum value of ssq
@minval = ();
# the maximum value of ssq
@maxval = ();

# read input file
while( <IFILE> ) 
{
 # read first and second numbers (feval and ssq), and ignore the rest (parameter values)
 if(/^([0-9]+)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/)
 {
  if( $1 < $lastline || $replicate == -1 )
  {
   # keep the largest counter
   if( $n > $top ) {$top = $n;}
   # this is new replicate, reset sequence counter
   $n = 0;
   # increment replicate count
   $replicate++;
  }
  if( $replicate == 0 )
  {
   # this only executes on the first replicate sequence
   # initialize counters and accumulators
   $count[$n] = 0;
   $mean[$n] = 0.0;
   $m2[$n] = 0.0;
   $minval[$n] = POSIX::FLT_MAX;
   $maxval[$n] = 0.0; # ssq is always positive...
  }
  # store the iteration number 
  $iter[$n] = $1;
  # remember its value until next line
  $lastline = $1;
  $count[$n] ++;
  $delta = $2 - $mean[$n];
  $mean[$n] += $delta / $count[$n];
  $delta2 = $2 - $mean[$n];
  $m2[$n] += $delta * $delta2;
  if( $minval[$n] > $2 ) {$minval[$n] = $2};
  if( $maxval[$n] < $2 ) {$maxval[$n] = $2};
  $n++
 }
}

# write out all the stats
print( O1FILE "#fevals\tMean\tStd dev\tN\tMin\tMax\n");
for( $n=0; $n < $top; $n++ )
{
 $stdev = sqrt($m2[$n] / $count[$n]);
 print( O1FILE "$iter[$n]\t$mean[$n]\t$stdev\t$count[$n]\t$minval[$n]\t$maxval[$n]\n");
}

close( O1FILE );
close( IFILE );
