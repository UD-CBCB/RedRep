#!/usr/bin/perl

my $ver="redrep-refmap.pl Ver. 2.3 [08/07/2020 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils;
use RedRep::Utils qw(cmd_STDOUT split_in_files);
use Getopt::Long qw(:config no_ignore_case);
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dircopy);
use File::Path qw(make_path remove_tree);
use Filesys::Df qw(df);
use Pod::Usage;
use Sys::Hostname qw(hostname);


### ARGUMENTS WITH NO DEFAULT
my($debug,$in,$outDir,$help,$manual,$force,$metaFile,$keep,$version,$diff,$refFasta,$pe,$stage,$tmpdir,$tmp_in_outdir,$no_stage_intermed);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu			=	1;
my $diff			= 0.04;
my $stop			= 100;
my $gap				= 1;
my $mem				= 5;						#in GB
our $verbose	= 4;
my $javaarg="";

GetOptions (
	"i|in=s"							=>	\$in,
	"o|out=s"							=>	\$outDir,
	"r|ref=s"							=>	\$refFasta,
	"l|log=s"							=>	\$logOut,

	"p|pe"								=>	\$pe,
	"n|ndiff=s"						=>	\$diff,
	"s|stop=i"						=>	\$stop,
	"g|gap=s"							=>	\$gap,

	"java=s"							=>	\$javaarg,

	"f|force"							=>	\$force,
	"d|debug:+"						=>	\$debug,
	"k|keep|keep_temp"		=>	\$keep,
	"V|verbose:+"					=>	\$verbose,
	"T|tmpdir=s"					=>	\$tmpdir,
	"S|stage"							=>	\$stage,
	"tmp_in_outdir"				=>	\$tmp_in_outdir,

	"t|threads|ncpu=i"		=>	\$ncpu,
	"mem=i"								=>	\$mem,

	"v|ver|version"				=>	\$version,
	"h|help"							=>	\$help,
	"m|man|manual"				=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);

if($debug) {
	require warnings; import warnings;
	$keep=1;
	$verbose=10;
	if($debug>1) {
		$verbose=100;
		require Data::Dumper; import Data::Dumper;
		require diagnostics; import diagnostics;
	}
	$tmp_in_outdir=1 if($debug>2);
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


### CREATE OUTPUT DIR AND OPEN LOG
mkdir($outDir) || die("ERROR: Can't create output directory $outDir");
$logOut="$outDir/log.refmap.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)",2);


### OUTPUT FILE LOCATIONS
$tmpdir								=	get_tmpdir($tmpdir);
my $intermed					=	"$tmpdir/intermed";
$intermed							=	"$outDir/intermed" if($tmp_in_outdir);
my $bam_dir						=	$outDir."/bam";
my $stage_dir					=	$tmpdir."input/";

mkdir($intermed);
mkdir($bam_dir);
mkdir($stage_dir);


### INITIATE LOG
logentry("Command: $0 $script\n",2);
logentry("Executing on ".hostname."\n",2);
logentry(find_job_info(),2);
logentry("Output directory $outDir (" . sprintf("%.2f",df("$outDir",1073741824)->{'bfree'}) . " GB free)\n",2);
logentry("Temporary directory $tmpdir (" . sprintf("%.2f",df("$tmpdir",1073741824)->{'bfree'}) . " GB free)\n",2);
logentry("Log verbosity: $verbose\n",2) if($verbose && $verbose>0);
logentry("Running in Debug Mode $debug\n",2) if($debug && $debug>0);
logentry("Keeping intermediate files.  WARNING: Can consume significant extra disk space\n",2) if($keep && $keep>0);


### REDREP BIN/SCRIPT LOCATIONS
my $execDir=$ENV{'REDREPBIN'};
my $utilDir=$ENV{'REDREPUTIL'};
my $libDir=$ENV{'REDREPLIB'};
logentry("RedRep Bin: $execDir\n",2);
logentry("RedRep Utilities: $utilDir\n",2);
logentry("RedRep Libraries: $libDir\n",2);


### CHECK FOR AND SETUP EXTERNAL DEPENDENCIES
logentry("Checking External Dependencies",3);

	# java
	#our $java=check_dependency("java","-version","s/\r?\n/ | /g",1);
	my $java_opts = "-Xmx${mem}g $javaarg -d64";

	# bwa
	my $bwa=get_path_bwa(1);

	# samtools
	my $samtools = get_path_samtools(1);

	# GATK
	my $gatk  = get_path_gatk(1);
	$gatk .= qq( --java-options "$java_opts");


### FIND FILES
logentry("COMPILING INPUT FILES",3);
my @files=split_in_files($in, qr/\.fastq$/);
logentry("FASTQ FILES DETECTED: ".scalar(@files),4);


### STAGE Files
if($stage) {
	logentry("STAGING FILES TO $tmpdir",3);
	@files=stage_files($stage_dir,\@files);
	logentry(scalar(@files)." input fastq files staged to $stage_dir",4);
}


############
### MAIN

	# ref Fasta stub
	my($stub2,$refPath)=fileparse($refFasta, qr/\.[^.]*$/);

foreach my $inFile(@files)
{
	logentry("BEGIN PROCESSING FILE $inFile",3);
	my $stub=fileparse($inFile, qr/\.[^.]*$/);

	### OUTPUT FILE LOCATIONS
	my $bwa_out="$intermed/$stub.bwa";
	my $bwa_bam="$intermed/$stub.bam";
	my $rg_bam="$intermed/$stub.rg.bam";
	my $sort_bam="$bam_dir/$stub.rg.sort.bam";

	# check Reference Fasta build indexes if needed
	check_ref_fasta($refFasta);

	if(! -e $refFasta.".amb" || ! -e $refFasta.".ann" || ! -e $refFasta.".bwt" || ! -e $refFasta.".pac" || ! -e $refFasta.".sa")
	{	logentry("REFERENCE FASTA BWT INDEXES NOT FOUND: BUILDING",4);
		$sys=cmd("$bwa index $refFasta","Build reference BWT indexes");
	}

	logentry("REFERENCE MAPPING",4);
	my $flags="";
	$flags="-p" if($pe);
	$sys=cmd_STDOUT("$bwa mem $flags -t $ncpu $refFasta $inFile 1> $bwa_out ","BWA Reference Mapping");

	logentry("MAKE BAM FILE",4);
	$sys=cmd("cat $bwa_out | $samtools view -bS -F 4 - -o $bwa_bam","Convert BWA to BAM");

	logentry("ADD READ GROUPS TO ALIGNMENT BAM",5);
	$sys=cmd("$gatk AddOrReplaceReadGroups --INPUT $bwa_bam --OUTPUT $rg_bam --SORT_ORDER coordinate --RGID $stub --RGLB $stub --RGPL illumina --RGPU $stub --RGSM $stub --VALIDATION_STRINGENCY LENIENT --MAX_RECORDS_IN_RAM 500000","Add read groups to BAM");

	logentry("SORT ALIGNMENT BAM FILE",5);
	$sys=cmd("$gatk ReorderSam --INPUT $rg_bam --OUTPUT $sort_bam --SEQUENCE_DICTIONARY $refFasta","Sorting BAM file");

	logentry("BUILDING BAM INDEX",5);
	$sys=cmd("$samtools index $sort_bam","Building BAM index");

	logentry("FINISH PROCESSING FILE $inFile",3);
}


### FILE CLEAN UP
logentry("FILE CLEAN UP",3);
if($keep && ! $tmp_in_outdir) {
	logentry("Saving intermediate tmp directory",4);
	$sys=dircopy($intermed, "$outDir/intermed");
}
logentry("Removing tmp files",4);
$sys=remove_tree($tmpdir);


### WRAP UP
logentry("SCRIPT COMPLETE",2);
logentry("TOTAL EXECUTION TIME: ".script_time(),2);
close($LOG);

exit 0;


#######################################
############### SUBS ##################



__END__


#######################################
########### DOCUMENTATION #############
=pod

=head1 NAME

redrep-refmap.pl -- De novo clustering of reduced representation libraries

=head1 SYNOPSIS

 redrep-refmap.pl --in FILE_OR_DIR_NAME --out DIRNAME --ref REFERENCE_FASTA [PARAMETERS]
                     [--help] [--manual]

=head1 DESCRIPTION

Accepts fastq output from redrep-qc.pl and fasta-formatted reference sequence(s), then maps fastq reads to fasta reference.

=head1 OPTIONS

=head2 REQUIRED PARAMETERS

=over 3

=item B<-i, --in>=FILENAME

Input fastq files.  (Required)

Can be one of the folllowing formats:

1. Path to one or more fastq files (comma separated; must have ".fastq" extension)
2. Path to one or more FOFN files (file with complete paths to one or more fastq files -- 1 per line; must have ".fofn" or ".fofn.list" extension)
3. Path to one or more directories of fastq files (comma separated)
4. Path to one or more FODN files (file with complete paths to directories containing one or more fastq files -- 1 directory per line; must have ".fodn" or ".fodn.list" extension)
5. Comma separated list of any combination of 1-4.

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA. Having a pre-built bwa index in the same directory will speed analysis.  If no bwa index is provided, script will build index and save to location of reference fasta (needs write permissions) (Required)

=back

=head2 OPTIONAL PARAMETERS

=head3 Program Specific Parameters

=over 3

=item B<-p, --pe>

bwa '-p' option: indicates paired end sequence files

=item B<-n, --ndiff>

bwa '-n' option: max #diff (int) or missing prob under 0.02 err rate (float) [0.04]

=item B<-s, --stop>

bwa '-R' option: stop searching when there are >INT equally best hits [100]

=item B<-g, --gap>

bwa '-o' option: maximum number or fraction of gap opens [1]

=back

=head3 Program Behavior and Resource Control

=over 3

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep>

Retain temporary intermediate files.

=item B<--mem>=integer

Max RAM usage for GATK Java Virtual Machine in GB (default 2)

=item B<-S, --stage>

When specified, input files will be staged to the tmpdir.  Can increase performance on clusters and other situations where the output directory is on network attached or cloud storage.

=item B<-t, --threads>=integer

Number of cpu's to use for threadable operations.

=item B<-T, --tmpdir>

Set directory (local directory recommended) for fast temporary and staged file operations.  Defaults to system environment variable REDREP_TMPDIR, TMPDIR, TEMP, or TMP in that order (if they exist).  Otherwise assumes /tmp.

=back

=head3 Logging Parameters

=over 3

=item B<-l, --log>=FILENAME

Log file output path. [ Default output-dir/log.map.txt ]

=item B<-V, --verbose>[=integer]

Produce detailed log.  Can be involed multiple times for additional detail levels or level can be specified. 0: ERROR, 1: WARNINGS, 2: INFO, 3-6: STATUS LEVELS, 7+:DEBUGGING INFO (default=4)

=back

=head3 Help

=over 3

=item B<-h, --help>

Displays the usage message.

=item B<-m, --man, --manual>

Displays full manual.

=item B<-v, --ver, --version>

Displays the current version.

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

=item 2.11 = 10/9/2019: Minor code cleanup.  Verbosity settings added.  Last version with GATK 3.x compatibility.

=item 2.2 - 10/28/2019: Added GATK v4 compatibility.  Many other enhancements including parallelization, documentation, and logging.

=item 2.3 - 8/1/2020: Fixed disk space reporting bugs.

=back

=head1 DEPENDENCIES

=head2 Requires the following external programs be in the system PATH:

=over 3

=item bwa (tested with version 0.7.16a r-1181)

=item samtools (tested with version 1.4.1)

=item java (tested with version 1.8.0_151-b12 OpenJDK Runtime Environment)

=item GATK (tested with version 4.1.4.0)

=back

=head1 AUTHOR

Written by Shawn Polson, University of Delaware

=head1 REPORTING BUGS

Report bugs to polson@udel.edu

=head1 COPYRIGHT

Copyright 2012-2020 Shawn Polson, Randall Wisser, Keith Hopper.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Please acknowledge author and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut
