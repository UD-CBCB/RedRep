#!/usr/bin/perl

# MANUAL FOR redrep-refmap.pl

=pod

=head1 NAME

redrep-refmap.pl -- De novo clustering of reduced representation libraries

=head1 SYNOPSIS

 redrep-refmap.pl --in FILE_OR_DIR_NAME --out DIRNAME [PARAMETERS]
                     [--help] [--manual]

=head1 DESCRIPTION

Accepts fastq output from redrep-qc.pl and fasta-formatted reference sequence(s), then maps fastq reads to fasta reference.

=head1 OPTIONS

=over 3

=item B<-1, -i, --in, --in1>=FILENAME

Input file in fastq format or directory of fastq files. (Required)

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA. (Required)

=item B<-l, --log>=FILENAME

Log file output path. [ Default output-dir/log.map.txt ]

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep_temp>

Retain temporary intermediate files.

=item B<-n, --ndiff>

bwa '-n' option: max #diff (int) or missing prob under 0.02 err rate (float) [0.04]

=item B<-s, --stop>

bwa '-R' option: stop searching when there are >INT equally best hits [100]

=item B<-g, --gap>

bwa '-o' option: maximum number or fraction of gap opens [1]

=item B<-t, --threads, --ncpu>=integer

Number of cpu's to use for threadable operations.

=item B<-d, --debug>

Produce detailed log.

=item B<-v, --ver, --version>

Displays the current version.

=item B<-h, --help>

Displays the usage message.

=item B<-m, --man, --manual>

Displays full manual.

=back

=head1 VERSION HISTORY

=over 3

=item 0.0 - 1/16/2013: Draft1

=item 1.0 - 1/23/2013: Release

=item 1.1 - 1/31/2013: Add parameter for number of gaps allowed

=item 1.2 - 3/27/2013: Fixed bug when ref fasta was in a different directory

=item 1.3 - 7/1/2013: Updated bwa path

=item 1.4 - 3/28/2014: Change bwa method from aln to mem; add version output to log; use path verisons of samtools, java, bwa; check ENV for $PICARDJARS and $REDREPBIN

-item 2.0 - 11/18/2016: Version 2 of picard tools now supported.  Major Release.

=back

=head1 DEPENDENCIES

=item bwa (tested with version 0.7.13-r1126)

=item samtools (tested with version 1.3.1)

=item picard-tools (tested with version 2.4.1)

=item java (tested with version 1.8.0_111-b16 OpenJDK Runtime Environment)

=head2 Requires the following Perl libraries:

=over 3

=item strict

=item Getopt::Long

=item File::Basename

=item Pod::Usage

=item POSIX

=back

=head2 Requires the following external programs be in the system PATH:

=over 3

=item

=back

=head1 AUTHOR

Written by Shawn Polson, University of Delaware

=head1 REPORTING BUGS

Report bugs to polson@udel.edu

=head1 COPYRIGHT

Copyright 2012-2016 Shawn Polson, Randall Wisser, Keith Hopper.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Please acknowledge author and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

my $script=join(' ',@ARGV);

use strict;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use POSIX;

sub binClass;
sub cmd;
sub logentry;


### ARGUMENTS WITH NO DEFAULT
my($in,$outDir,$help,$manual,$force,$metaFile,$keep,$debug,$version,$diff,$refFasta);


### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu		=	1;
my $diff=0.04;
my $stop=100;
my $gap=1;



GetOptions (	"i|in=s"					=>	\$in,
				"o|out=s"					=>	\$outDir,
				"r|ref=s"					=>	\$refFasta,
				"l|log=s"					=>	\$logOut,

				"n|ndiff=s"					=>	\$diff,
				"s|stop=i"					=>	\$stop,
				"g|gap=s"					=>	\$gap,

				"f|force"					=>	\$force,
				"d|debug"					=>	\$debug,
				"k|keep_temp"				=>	\$keep,

				"t|threads|ncpu=i"			=>	\$ncpu,

				"v|ver|version"				=>	\$version,
				"h|help"					=>	\$help,
				"m|man|manual"				=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
my $ver="redrep-refmap.pl Ver. 2.0 (11/18/2016 rev)";
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);
#pod2usage( -msg  => "ERROR!  Required argument -m (metadata file) not found.\n", -exitval => 2)  if (! $metaFile);

if($debug)
{	require warnings; import warnings;
	require Data::Dumper; import Data::Dumper;
	$keep=1;
}


### DECLARE OTHER GLOBALS
my $sys;												# system call variable
#my $stub=fileparse($inFile, qr/\.[^.]*(\.gz)?$/);


### THROW ERROR IF OUTPUT DIRECTORY ALREADY EXISTS (unless $force is set)
if(-d $outDir)
{	if(! $force)
	{	pod2usage( -msg  => "ERROR!  Output directory $outDir already exists.  Use --force flag to overwrite.", -exitval => 2);
	}
	else
	{	$sys=`rm -R $outDir`;
	}
}


### OUTPUT FILE LOCATIONS
my $intermed="$outDir/intermed";

### CREATE OUTPUT DIR
mkdir($outDir);
mkdir($outDir."/intermed");
my(@files);
if(-d $in)
{	opendir(DIR,$in);
	@files=grep /\.fastq$/, readdir(DIR);
	close(DIR);
	s/^/$in\// for @files;  #prepend $in (dirpath) to each element
}
elsif (-f $in)
{	push(@files,$in);
}
else
{	pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) file not found.\n", -exitval => 2) if (! $in);
}


### BIN/SCRIPT LOCATIONS
my %ver;
my %path;
my $execDir=$ENV{'REDREPBIN'};
my $java = "java -Xmx2g -d64 -jar ";
$ver{'java'}=`java -version 2>&1|grep version`;
$path{'java'}=`which java`;
my $bwa="bwa";
$ver{'bwa'}=`$bwa 2>&1|grep Version`;
$path{'bwa'}=`which $bwa`;
my $samtools="samtools";
$ver{'samtools'}=`$samtools 2>&1|grep Version`;
$path{'samtools'}=`which $samtools`;
my $picard = $ENV{'PICARDJARS'}."/picard.jar";
$ver{'picard'}=`$java $picard AddCommentsToBam --version 2>&1`;
chomp %ver;
chomp %path;


### CREATE LOG FILES
$logOut="$outDir/log.map.txt" if (! $logOut);
open(our $LOG, "> $logOut");
print $LOG "$0 $script\n";
print $LOG "RedRep Scripts: $execDir\n";
print $LOG "Dependency: $path{'java'} ($ver{'java'})\n";
print $LOG "Dependency: $path{'bwa'} ($ver{'bwa'})\n";
print $LOG "Dependency: $path{'samtools'} ($ver{'samtools'})\n";
print $LOG "Dependency: $picard ($ver{'picard'})\n";
logentry("SCRIPT STARTED ($ver)");


### CHECK FOR EXTERNAL DEPENDENCIES
# Not yet implemented



############
### MAIN

	# ref Fasta stub
	my($stub2,$refPath)=fileparse($refFasta, qr/\.[^.]*$/);

foreach my $inFile(@files)
{
	my $stub=fileparse($inFile, qr/\.[^.]*$/);

	### OUTPUT FILE LOCATIONS
	my $bwa_out="$intermed/$stub.bwa";
	my $bwa_bam="$intermed/$stub.bam";
	my $rg_bam="$intermed/$stub.rg.bam";
	my $sort_bam="$outDir/$stub.rg.sort.bam";

	if(! -e $refFasta.".fai")
	{	logentry("REFERENCE FASTA INDEX NOT FOUND: BUILDING");
		$sys=cmd("$samtools faidx $refFasta","Building reference index");
	}
	if(! -e "$refPath/$stub2.dict")
	{	logentry("REFERENCE DICTIONARY NOT FOUND: BUILDING");
		$sys=cmd("$java $picard CreateSequenceDictionary R=$refFasta O=$refPath/$stub2.dict","Building reference dictionary");
	}

	if(! -e $refFasta.".amb" || ! -e $refFasta.".ann" || ! -e $refFasta.".bwt" || ! -e $refFasta.".pac" || ! -e $refFasta.".sa")
	{	logentry("REFERENCE FASTA BWT INDEXES NOT FOUND: BUILDING");
		$sys=cmd("$bwa index $refFasta","Build reference BWT indexes");
	}

	logentry("REFERENCE MAPPING");
#	$sys=cmd("$bwa aln -t $ncpu -n $diff -R $stop -o $gap -f $bwa_out $refFasta $inFile","BWA Reference Mapping");
	$sys=cmd2("$bwa mem -t $ncpu $refFasta $inFile 1> $bwa_out ","BWA Reference Mapping");

	logentry("MAKE BAM FILE");
#	$sys=cmd("$bwa samse $refFasta $bwa_out $inFile | $samtools view -bS -F 4 - -o $bwa_bam","Convert BWA to BAM");
	$sys=cmd("cat $bwa_out | $samtools view -bS -F 4 - -o $bwa_bam","Convert BWA to BAM");

	logentry("ADD READ GROUPS TO ALIGNMENT BAM");
	$sys=cmd("$java $picard AddOrReplaceReadGroups I=$bwa_bam O=$rg_bam SORT_ORDER=coordinate RGID=$stub RGLB=$stub RGPL=illumina RGPU=$stub RGSM=$stub VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=500000","Add read groups to BAM");

	logentry("SORT ALIGNMENT BAM FILE");
	$sys=cmd("$java $picard ReorderSam I=$rg_bam O=$sort_bam REFERENCE=$refFasta","Sorting BAM file");

	logentry("BUILDING BAM INDEX");
	$sys=cmd("$samtools index $sort_bam","Building BAM index");
}

	# FILE CLEAN UP
	logentry("FILE CLEAN UP");
	if(! $keep || ! $debug)
	{	$sys=cmd("rm -R $intermed","Remove intermediate file directory");
	}



logentry("SCRIPT COMPLETE");
close($LOG);

exit 0;


#######################################
############### SUBS ##################


#######################################
### binClass
# run system command and collect output and error states
sub binClass

{	my $number=shift;
	my $size=shift;
	my $start=shift;

	my $bottom=$start+($number*$size);
	my $top=($start+(($number+1)*$size))-1;
	return $bottom."-".$top;
}




#######################################
### cmd
# run system command and collect output and error states
sub cmd
{	my $cmd=shift;
	my $message=shift;

	logentry("System call: $cmd") if($debug);

	my $sys=`$cmd 2>&1`;
	my $err=$?;
	if ($err)
	{	print $LOG "$message\n$cmd\nERROR $err\n$sys\n" if($LOG);
		pod2usage( -msg  => "ERROR $err!  $message\n", -exitval => 2);
		return 1;
	}
	else
	{	print $LOG $sys;
		return $sys;
	}
}

#######################################
### cmd2
# run system command and collect output and error states (when std out must be used separately)
sub cmd2
{	my $cmd=shift;
	my $message=shift;

	logentry("System call: $cmd") if($debug);

	my $sys=`$cmd 2> err`;
	my $err=$?;
	if ($err)
	{	$sys=`cat err`;
		print $LOG "$message\n$cmd\nERROR $err\n$sys\n" if($LOG);
		pod2usage( -msg  => "ERROR $err!  $message\n", -exitval => 2);
		return 1;
	}
	else
	{	print $LOG $sys;
		return $sys;
	}
	$sys=`rm err`;
}


#######################################
### logentry
# Enter time stamped log entry
sub logentry
{	my $message=shift;
	print $LOG POSIX::strftime("%m/%d/%Y %H:%M:%S > $message\n", localtime);
}



__END__
