#!/usr/bin/perl
#RedRep Utility fasta-abbrev.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser
#ARGS inputfile-fasta outputfile-fasta outputfile-dictionary

open(DAT,$ARGV[0]);
open(NEW,"> ".$ARGV[1]);
open(DICT,"> ".$ARGV[2]);
my $i=0;
while(<DAT>)
{
	if($_=~/^>/)
	{	$_=~/^>(.+SAMPLE=(\S+)\s.+)$/;
		print DICT "$2-$i\t$1\n";
		print NEW ">$2-$i\n";
		$i++;
	}
	else
	{	print NEW $_;
	}
}
close(DAT);
close(NEW);
close(DICT);
