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

=item B<-p, --pe>

bwa '-p' option: indicates paired end sequence files

=item B<-n, --ndiff>

bwa '-n' option: max #diff (int) or missing prob under 0.02 err rate (float) [0.04]

=item B<-s, --stop>

bwa '-R' option: stop searching when there are >INT equally best hits [100]

=item B<-g, --gap>

bwa '-o' option: maximum number or fraction of gap opens [1]

=item B<--mem>=integer

Max RAM usage for GATK/Picard Java Virtual Machine in GB (default 2)

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

=item 2.0 - 11/18/2016: Version 2 of picard tools now supported.  Major Release.

=item 2.1 - 2/2/2017: Code cleanup including minimizing OS dependencies

=back

=head1 DEPENDENCIES

=head2 Requires the following Perl Modules:

=head3 Non-Core (Not installed by default in some installations of Perl)

=over 3

=item File::Which

=item Parallel::ForkManager

=back

=head3 Core Modules (Installed by default in most Perl installations)

=over 3

=item strict

=item Exporter

=item File::Basename

=item File::Copy

=item File::Path

=item Getopt::Long

=item Pod::Usage

=item POSIX

=item Sys::Hostname

=back

=head2 Requires the following external programs be in the system PATH:

=over 3

=item bwa (tested with version 0.7.16a r-1181)

=item samtools (tested with version 1.4.1)

=item picard-tools (tested with version 2.19.0 SNAPSHOT)

=item java (tested with version 1.8.0_151-b12 OpenJDK Runtime Environment)

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

# standard libs
use strict;
use Getopt::Long;
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);
use Pod::Usage;
use Sys::Hostname qw(hostname);

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils qw(check_dependency check_jar cmd cmd_STDOUT find_job_info get_execDir logentry logentry_then_die);


### ARGUMENTS WITH NO DEFAULT
my($in,$outDir,$help,$manual,$force,$metaFile,$keep,$version,$diff,$refFasta,$pe,$javaarg);
our($debug);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu		=	1;
my $diff=0.04;
my $stop=100;
my $gap=1;
my $mem=2;						#in GB


GetOptions (	"i|in=s"					=>	\$in,
				"o|out=s"					=>	\$outDir,
				"r|ref=s"					=>	\$refFasta,
				"l|log=s"					=>	\$logOut,

				"p|pe"						=>	\$pe,
				"n|ndiff=s"					=>	\$diff,
				"s|stop=i"					=>	\$stop,
				"g|gap=s"					=>	\$gap,

				"java=s"													=>	\$javaarg,

				"f|force"					=>	\$force,
				"d|debug"					=>	\$debug,
				"k|keep_temp"				=>	\$keep,

				"t|threads|ncpu=i"			=>	\$ncpu,
				"mem=i"									=>	\$mem,

				"v|ver|version"				=>	\$version,
				"h|help"					=>	\$help,
				"m|man|manual"				=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
my $ver="redrep-refmap.pl Ver. 2.1 (2/1/2017 rev)";
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);

if($debug)
{	require warnings; import warnings;
	require Data::Dumper; import Data::Dumper;
	$keep=1;
}


### SET DEFAULT METHOD OF FILE PROPAGATION
my $mv = \&move;
if($keep)
{	$mv = \&copy;
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
	{	$sys=remove_tree($outDir);
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


### CREATE LOG FILES
$logOut="$outDir/log.map.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)");
print $LOG "Command: $0 $script\n";
print $LOG "Executing on ".hostname."\n";
print $LOG find_job_info()."\n";


### REDREP BIN/SCRIPT LOCATIONS
my $execDir=$ENV{'REDREPBIN'};
my $utilDir=$ENV{'REDREPUTIL'};
my $libDir=$ENV{'REDREPLIB'};
print $LOG "RedRep Bin: $execDir\n";
print $LOG "RedRep Utilities: $utilDir\n";
print $LOG "RedRep Libraries: $libDir\n";


### CHECK FOR AND SETUP EXTERNAL DEPENDENCIES
logentry("Checking External Dependencies");

	# java
	our $java=check_dependency("java","-version","s/\r?\n/ | /g");
	$java = "${java} -Xmx${mem}g $javaarg -d64 -jar ";

	# bwa
	my $bwa=check_dependency("bwa",' 2>&1 |grep "Version"',"s/\r?\n/ | /g");

	# samtools
	my $samtools=check_dependency("samtools","--version","s/\r?\n/ | /g");

	#picard
	# Picard tools (as of v 2.19.0) gives an error code of 1 when version number is checked.  To get around the $? check in &cmd(), added " 2>&1;v=0" to the version check flag to make error code 0, but still get version redirected to SDTOUT
	my $picard=check_jar("Picard Tools",$ENV{'PICARDJAR'},"AddCommentsToBam --version 2>&1;v=0","s/\r?\n/ | /g","Environment Variable PICARDJAR is not defined or is a not pointing to a valid file!  Please define valid picard tools location with command: export PICARDJAR='PATH_TO_PICARD_JAR");


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
	my $flags="";
	$flags="-p" if($pe);
	$sys=cmd_STDOUT("$bwa mem $flags -t $ncpu $refFasta $inFile 1> $bwa_out ","BWA Reference Mapping");

	logentry("MAKE BAM FILE");
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
	{	$sys=remove_tree($intermed);
	}



logentry("SCRIPT COMPLETE");
close($LOG);

exit 0;


#######################################
############### SUBS ##################



__END__
