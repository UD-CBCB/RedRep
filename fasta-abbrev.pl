#!/usr/bin/perl

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
