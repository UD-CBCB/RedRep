#!/usr/bin/perl
#RedRep Utility fasta-abbrev.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser
#ARGS inputfile-fasta outputfile-fasta outputfile-dictionary
#Reduces fasta header and chomps seqeunce to 1 line to get smallest possible filesize
#Sequence is renamed using the SAMPLE= field plus an iterator and adds a size tag for usearch input (SAMPLE-#;size=1;)
#Produces dictionary of headers, so that original sequence name can be restored if desired
#Typical use is by redrep-cluster.pl as input for usearch clustering

open(DAT,$ARGV[0]);
open(NEW,"> ".$ARGV[1]);
open(DICT,"> ".$ARGV[2]);
my $i=0;
my $first=1;
while(<DAT>)
{	chomp;				#chomp ensures 1 line fasta sequence -- required for redrep-cluster if split is needed
	if($_=~/^>/)
	{	print NEW "\n" if($first==0);
		$first=0;
		$_=~/^>(.+SAMPLE=(\S+)\s.+)$/;
		print DICT "$2-$i\t$1\n";
		print NEW ">$2-${i};size=1;\n";
		$i++;
	}
	else
	{	print NEW $_;
	}
}
close(DAT);
close(NEW);
close(DICT);
