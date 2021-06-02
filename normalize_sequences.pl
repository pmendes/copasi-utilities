#!/usr/bin/perl
#
# normalize_sequences.pl inputfile {-v} {interval {final}}
#
# Written by Pedro Mendes, May 2021. Published under MIT License.
# https://raw.githubusercontent.com/pmendes/copasi-utilities/main/LICENSE
#
# Normalizes irregularly spaced sequences of data to regularly spaced 
# sequences with a specified interval. This is primarily intended for 
# files written by COPASI in Parameter Estimation or Optimization 
# tasks where the software only outputs a data point where there is
# an improvement of the optimization task. This results in files 
# that have irregular spacings in terms of function evaluations.
# This script processes those files and outputs a corresponding file
# where the data is regularly spaced.
#
# Parameters:
#  inputfile: file to be processed
#  interval: regular spacing for output; default is 10
#  final: value requested for end of the sequence; if none is passed
#         data will only go up to the last value present in inputfile
#
# Input file is expected to contain at least two numbers per line: the
# first is the counter for function evaluations, the second for the
# objective function. Lines can contain any number of extra information
# which will be passed on to the output too. Lines without numbers are
# passed verbatim to the output file.
#
# Output file is named based on inputfile; examples:
#  input: filename           output: filename.norm.tsv
#  input: basename.tsv       output: basename.norm.tsv   
#  input: basename.ext       output: basename.norm.tsv   

if( $#ARGV < 0 ) { die "usage: normalize_sequences.pl inputfile {-v} {interval {final}}\n"; }
$in = $ARGV[0];

$verbose=0;

$interval = 10;
if( $#ARGV > 0) 
{
 $arg = $ARGV[1];
 if( $arg =~ /-vv/ ) 
 {
  $verbose = 2;
 }
 else 
 {
  if( $arg =~ /-v/ ) 
  {
   $verbose = 1;
  }
  else
  {
   if( $arg =~ /\d+/ ) 
   {
    $interval = $arg;
   }
   else { die "usage: normalize_sequences.pl inputfile {-v} {interval {final}}\n"; }
  }
 }
}

$final = 0;
if( $#ARGV > 1) 
{
 $arg = $ARGV[2];
 if( $verbose )
 {
  # last line was -v, so this must be interval
  if( $arg =~ /\d+/ ) 
  {
   $interval = $arg;
  }
  else { die "usage: normalize_sequences.pl inputfile {-v} {interval {final}}\n"; }
 }
 else
 {
  # since there was no -v this argument must be $final
  if( $arg =~ /\d+/ )
  { 
   $final = $arg;
  }
  else { die "usage: normalize_sequences.pl inputfile {-v} {interval {final}}\n"; }
 }
}

if( $verbose && $#ARGV > 2) 
{
 # we had -v so we can still read $final
 $arg = $ARGV[3];
 if( $arg =~ /\d+/ )
 { 
  $final = $arg;
 }
 else { die "usage: normalize_sequences.pl inputfile {-v} {interval {final}}\n"; }
}

# extract filename without extension if possible
$basename = $in;
$basename =~ s/(.+)\..*$/$1/;

$normf = "$basename.norm.tsv";

open( IFILE, "$in" )  || die "can't read $in";

# file for normalized sequence
open( O1FILE, ">$normf" )  || die "can't open $normf";

$ordinal = 0; # counter for current sequence
$lastprinted = 0; # last sequence counter seen
$lastrest = ""; # last rest of line seen
$lastread = 0; #last line read

while( <IFILE> ) 
{
 # read first number and rest of line
 if(/^([0-9]+)(.*)$/)
 {
  $ordinal = $1;
  if( $verbose) {print "reading $ordinal\n";}
  $rest = $2;

  # first we check if the new value is smaller than the last one
  # and if the last printed is smaller than the last read; because
  # in that case we still need to print the last value 
  if( $ordinal<$lastprinted && $lastprinted < $lastread)
  {
   # we print the last one read in the corresponding interval
   $i = $lastprinted + $interval;
   print( O1FILE "$i$lastrest\n");
   if( $verbose) {print "printing $i with value $lastread\n";}
   $lastprinted = $i;
  }  
  
  # check if we are in a new sequence and if we still need to print
  # the last sequence until the final line requested
  if( $ordinal<$lastread && $final>$lastprinted )
  {
   for($i = $lastprinted+$interval; $i <= $final; $i+=$interval)
   {
    print( O1FILE "$i$lastrest\n");
    if( $verbose) {print "printing $i with value $ordinal\n";}
    $lastprinted = $i;
   }
  }
 
  # if $ordinal is the start of the sequence let's process it now
  if( $ordinal==1 || $ordinal<$lastread )
  {
   # print as many intervals as needed until reaching $ordinal
   for($i = 1; $i <= $ordinal; $i+=$interval)
   {
    print( O1FILE "$i$rest\n");
    if( $verbose) {print "printing $i with value $ordinal\n";}
    $lastprinted = $i;
   }
   $lastrest = $rest;
   $lastread = $ordinal;
   next;
  }
  $steps = $ordinal - $lastprinted;
  if( $steps < $interval )
  {
    # if step less than interval, then we just keep the rest of the line but do not output now
    $lastrest = $rest;
    $lastread = $ordinal;
  }
  else
  {
   if( $steps == $interval )
   {
    # exact step size as requested, write out this line as is
    print( O1FILE "$_");
    if( $verbose) {print "printing $ordinal with value $ordinal\n";}
    $lastprinted = $ordinal;
    $lastread = $ordinal;
    $lastrest = $rest;
   } 
   else
   {
    # larger step taken than required, keep writing lines with $oldrest until needed
    for($i = $lastprinted + $interval; $i < $ordinal; $i+=$interval)
    {
     print( O1FILE "$i$lastrest\n");
     if( $verbose) {print "printing $i with value $lastread\n";} 
     $lastprinted = $i;
    }
    $lastrest = $rest;
    $lastread = $ordinal;
    # now check if we have an exact step, in which case we print this
    if( $ordinal - $lastprinted == $interval )
    {
     # exact step size as requested, write out this line as is
     print( O1FILE "$_");
     if( $verbose) {print "printing $ordinal with value $ordinal\n";}
     $lastprinted = $ordinal;
    }
   }
  }
 }
 else
 {
  # this line doen't start with a number, so we first check
  # if last printed is smaller than the last read; because 
  # in that case we still need to print the last value 
  if( $lastprinted < $lastread)
  {
   # we print the last one read in the corresponding interval
   $i = $lastprinted + $interval;
   print( O1FILE "$i$lastrest\n");
   if( $verbose) {print "printing $i with value $lastread\n";}
   $lastprinted = $i;
  }
  # now it this line is empty it signals the end of a sequence 
  # so if needeed, let's print it until the end
  if( !/\S/ && $final>$lastprinted )
  {
   for($i = $lastprinted+$interval; $i <= $final; $i+=$interval)
   {
    print( O1FILE "$i$lastrest\n");
    if( $verbose) {print "printing $i with value $lastread\n";}
    $lastprinted = $i;
   }
  }
  # finally we just print the line as it came in
  print( O1FILE "$_");
  if( $verbose>1) {print "printing non-numeric line\n";}
 }
}

# before we leave we check if last printed is smaller than 
# the last read; because in that case we still need to 
# print the last value 
if( $lastprinted < $lastread)
{
 # we print the last one read in the corresponding interval
 $i = $lastprinted + $interval;
 print( O1FILE "$i$lastrest\n");
 if( $verbose) {print "printing $i with value $lastread\n";}
 $lastprinted = $i;
}

# now check if we need print extra lines until the end
if( $final>$lastprinted )
{
 for($i = $lastprinted+$interval; $i <= $final; $i+=$interval)
 {
  print( O1FILE "$i$lastrest\n");
  if( $verbose) {print "printing $i with value $lastread\n";}
  $lastprinted = $i;
 }
}  

# close the files, we're done
close( O1FILE );
close( IFILE );
