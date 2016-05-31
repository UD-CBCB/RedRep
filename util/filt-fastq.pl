#!/usr/bin/perl
#RedRep Utility filt-fastq.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser
#reads from stdin
#ARGS min_len min_qual

while(<STDIN>)
{	$_=~/LENGTH=(\d+)\sMEAN_QUAL=(\S+)/;
	print $_ if($1>=$ARGV[0] && $2>=$ARGV[1]);
}
