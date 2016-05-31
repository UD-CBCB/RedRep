#!/usr/bin/perl -w
use strict;

sub promptUser;

#	Accepts a cd-hit or psi-cd-hit .clstr file and determines library	
#		composition per cluster based on LibName (below).				
#																		
#	INPUT: .clstr file with seq descriptions as follows:				
#			[LibName-Identifier] or [LibName_Identifier]				
#			LibName cannot include spaces, hyphens, or underscores		
#																		
#	OUTPUT: tab delimited table											
#																		
#	ARGUMENTS: none; (user is prompted for input/ouput file names)		
#																		
#	Shawn Polson, 5-2-2008

my($inputFile, $outputFile, $clustSize, $line, $j, $key);
my @clustComp;
my @clustSeed;
my %keyVals;
my $pass=0;
my $clustNum=-1;
#my $delim="mge:[^:]+:[^:]+";		# aclame clusters
my $delim='\w\w\w';					# metagenome sequences
#my $delim="[^\-]+";					# ncbi viral genomes

$inputFile = $ARGV[0]; # &promptUser("Enter CD-Hit .clstr file");
open(DAT, $inputFile) or die "Cannot open input file $inputFile\n";

$outputFile = $ARGV[1]; # &promptUser("Enter output file","$inputFile.freq");
open(OUT, "> $outputFile") or die "Cannot open output file $outputFile\n";

while(<DAT>)
{	$line=$_;
	if($line =~ "^>")
	{	$clustNum++;
	} 
	
	# READ CLUSTER COMPONENT ENTRIES
#	if ($line =~ /^\d+\D+\d+aa.+>([^\_^\-]+)[\-\_]/)
#	if ($line =~ /\d+[na][ta].+>(\w\w\w)/)
	if ($line =~ /\d+[na][ta], >($delim)/)
	{	
		$clustComp[$clustNum]{$1}++;
		$keyVals{$1}=1;
		($clustSeed[$clustNum])=($line=~/\d+[na][ta], >(\S+)/) if($line=~/\*$/);
	}
}  #END WHILE(<DATA>) LOOP
close(DAT);

#PRINT LIBRARY FREQUENCY FILE
#PRINT HEADERS
print OUT "\t";
foreach $key (sort keys %keyVals) 
{	print OUT "$key\t";
}
print OUT "\n";
#PRINT DATA
for($j=0; $j<=$clustNum; $j++)
{	$pass=1;
	print OUT "Cluster_$j:$clustSeed[$j]\t";
	foreach $key (sort keys %keyVals) 
	{	if($pass==2)
		{	print OUT "\t";
		}
		if(defined($clustComp[$j]{$key}))
		{	print OUT $clustComp[$j]{$key};
		}
		else
		{	print OUT "0";
		}
		$pass=2;
	}
	print OUT "\n";
}
close(OUT);



#----------------------------(	promptUser	)-----------------------------#
#																		  #
#  FUNCTION:	promptUser												  #
#																		  #
#  PURPOSE: Prompt the user for some type of input, and return the		  #
#		input back to the calling program.								  #
#																		  #
#  ARGS:	$promptString - what you want to prompt the user with		  #
#		$defaultValue - (optional) a default value for the prompt		  #
#																		  #
#-------------------------------------------------------------------------#

sub promptUser 
{

	my ($promptString, $defaultValue) = @_;
	if ($defaultValue) 
	{	print $promptString, "[", $defaultValue, "]: ";
	} 
	else 
	{	print $promptString, ": ";
	}
	$| = 1;				  # force a flush after our print
	$_ = <STDIN>;		  # get the input from STDIN (presumably the keyboard)
	chomp;
	if (defined($defaultValue)) 
	{	return $_ ? $_ : $defaultValue;	   # return $_ if it has a value
	}
	else 
	{	return $_;
	}
}
