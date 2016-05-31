#!/usr/bin/perl
#RedRep Utility fastq2fasta.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser

$i=0;
while(<>)
{	if(/^\@/&&$i==0)
	{	s/^\@/\>/;
		print;
	}
	elsif($i==1)
	{	print;
		$i=-3;
	}
	$i++;
}
