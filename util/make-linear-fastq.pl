#!/usr/bin/perl
#RedRep Utility make-linera-fastq.pl
#Copyright 2016 Shawn Polson, Keith Hopper, Randall Wisser
#reads from STDIN
#outputs to STDOUT


use strict;

my $count = 0;
while(my $line = <STDIN>){
    chomp($line);
    print $line;
    $count = ($count + 1) % 4;

    if($count == 0){
        print "\n";
    }else{
        print "\t";
    }
}
