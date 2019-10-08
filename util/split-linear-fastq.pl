#!/usr/bin/perl -w
#RedRep Utility split-linear-fastq.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser
#reads from STDIN
#outputs to STDOUT

use strict;

my $count = 0;
while(my $line = <STDIN>){
    chomp($line);
    my @vals = split(/\t/, $line);

    print $vals[0]."\n";
    print $vals[1]."\n";
    print $vals[2]."\n";
    print $vals[3]."\n";
}
