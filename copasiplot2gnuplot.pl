#!/usr/bin/perl -w
#
# adds two blank lines in place of any line that has't got any numbers
# this is useful to convert COPASI plot output files into gnuplot 
# files with indexed groups (changes lines of nan to two blank lines)

if( $#ARGV < 1 ) { die "usage: copasiplot2gnuplot.pl inputfile outputfile\n"; }
$in = $ARGV[0];
$out = $ARGV[1];

open( IFILE, "$in" )  || die "can't read $in";
open( OFILE, ">$out" )  || die "can't write $out";

while( <IFILE> ) 
{
 if(/[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?/)
 {
  print OFILE $_;
 }
 else
 {
  print OFILE "\n\n";
 }
}

close $out;