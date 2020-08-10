#!/usr/bin/perl

my $ver="redrep-SNPcall.pl Ver. 2.3 [08/07/2020 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils;
use RedRep::Utils qw(build_argument_list build_contig_list read_fofn read_listfile split_in_files);
use Getopt::Long qw(:config no_ignore_case);
use Parallel::ForkManager;
use Pod::Usage;
use POSIX qw(ceil floor);
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dircopy);
use File::Path qw(make_path remove_tree);
use Filesys::Df qw(df);
use Sys::Hostname qw(hostname);


### ARGUMENTS WITH NO DEFAULT
my($debug,$in,$outDir,$help,$manual,$force,$metaFile,$keep,$version,$refFasta,$dcov,$dfrac,$dt,
    $intervals,$stage,$tmpdir,$tmp_in_outdir,$no_stage_intermed);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu			=	1;
my $mem				= 50;						#in GB
our $verbose	= 4;
my $javaarg="";
my $gatkarg_HaplotypeCaller="";

GetOptions (
	"i|in=s"													=>	\$in,
	"o|out=s"													=>	\$outDir,
	"r|ref=s"													=>	\$refFasta,
	"l|log=s"													=>	\$logOut,

	"t|threads|ncpu=i"								=>	\$ncpu,
	"mem=i"														=>	\$mem,
	"T|tmpdir=s"											=>	\$tmpdir,
	"S|stage"													=>	\$stage,
	"tmp_in_outdir"										=>	\$tmp_in_outdir,

	"java=s"													=>	\$javaarg,
	"gatk_HaplotypeCaller=s"					=>	\$gatkarg_HaplotypeCaller,
	"L|intervals=s"										=>	\$intervals,

	"f|force"													=>	\$force,
	"d|debug:+"												=>	\$debug,
	"k|keep|keep_temp"								=>	\$keep,
	"V|verbose:+"											=>	\$verbose,

	"v|ver|version"										=>	\$version,
	"h|help"													=>	\$help,
	"m|man|manual"										=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);
pod2usage( -msg  => "ERROR!  Arguments -dcov and -dfrac are incompatible, chose one.\n", -exitval => 2)  if ($dcov && $dfrac);


### DEBUG MODE
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
my $hmm_threads=1;							# 1 hmm thread with 8 splits is marginally faster than 2hmm threads with 4 splits (within 5%)
my $manager = new Parallel::ForkManager( floor($ncpu/$hmm_threads) );

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
$logOut="$outDir/log.SNPcall.txt" if (! $logOut);
open(our $LOG,">",$logOut);
logentry("SCRIPT STARTED ($ver)",2);


### OUTPUT FILE LOCATIONS
$tmpdir								=	get_tmpdir($tmpdir);
my $intermed					=	"$tmpdir/intermed";
$intermed							=	"$outDir/intermed" if($tmp_in_outdir);
my $gvcf_dir					=	$outDir."/gvcf";
my $stage_dir					=	$tmpdir."input/";
my $contig_list				= $intermed."/contig.list";

mkdir($intermed);
mkdir($gvcf_dir);
mkdir($stage_dir);


### FIND FILES
my @files=split_in_files($in, qr/\.bam$/);
logentry("BAM FILES DETECTED: ".scalar(@files),3);
my @bai_files;


### STAGE Files
if($stage) {
	logentry("STAGING FILES TO $tmpdir",3);

	# bai not found are built later on tmpdir (in case of write permission issues?)

	foreach my $file (@files) {
		if(-e "$file.bai" && -f "$file.bai") {
			push(@bai_files,"$file.bai");
		}
	}

	@files=stage_files($stage_dir,\@files);
	logentry(scalar(@files)." input bam files staged to $stage_dir",4);

	@bai_files=stage_files($stage_dir,\@bai_files);
	logentry(scalar(@bai_files)." input bam index files staged to $stage_dir",4);

}

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
	#our $java=check_dependency("java","-version","s/\r?\n/ | /g");    #not currently needed; retain for future
	my $java_opts = "-Xmx${mem}g $javaarg -d64";

	# samtools
	my $samtools = get_path_samtools(1);

	#GATK
	my $gatk  = get_path_gatk(1);
	$gatk .= qq( --java-options "$java_opts");


### OUTPUT FILE LOCATIONS
	my $gvcf_fofn_out			=	"$outDir/gvcf.fofn.list";
	my $contig_list				= $tmpdir."/contig.list";


############
### MAIN

	# ref Fasta stub
	my $stub2=fileparse($refFasta, qr/\.[^.]*$/);

	## CHECK AND BUILD FASTA INDEX AND DICTIONARY
	check_ref_fasta($refFasta);
	build_contig_list($refFasta,$contig_list);

	## BAM CHECK AND BUILD INDEXES
	foreach my $file (@files) {
		unless(-e "$file.bai" && -f "$file.bai") {
			logentry("BAM Index not found for $file",3);
			logentry("BUILDING BAM INDEX $file.bai",4);
			$sys=cmd("$samtools index $file","Calling samtools to building BAM index");
			push(@bai_files,"$file.bai");
		}
	}


	#my $bam_files=build_argument_list(@files,"--input");
	my $interval_args="";
	my @intervals;
	if($intervals) {
		if($intervals=~/.list$/ && -e $intervals && -f $intervals) {
			@intervals=read_listfile($intervals);
		}
		else {
			@intervals=split(/,/,$intervals);
		}
		$interval_args=build_argument_list(\@intervals,"--intervals");
	}

  logentry("GATK VARIANT CALLING: ERC GVCF SINGLE-SAMPLE DISCOVERY MODE",3);

	open(FOFN,">",$gvcf_fofn_out);
	foreach my $file (@files)
	{	$manager->start and next;
		logentry("BEGIN PROCESSING FILE $file",4);
		my $stub=fileparse($file, qr/\.[^.]*(\.gz)?$/);
		my $outFile="$gvcf_dir/$stub.g.vcf";
		$sys=cmd("$gatk HaplotypeCaller --input $file --output $outFile --reference $refFasta -ERC GVCF --native-pair-hmm-threads $hmm_threads -RF GoodCigarReadFilter $interval_args $gatkarg_HaplotypeCaller","Run HaplotypeCaller on ${file}");
		print FOFN "$outFile\n";
		if(-e $outFile && -f $outFile) {
			logentry("FINISH PROCESSING FILE $file",4);
		}
		else {
			logentry("No g.vcf file created for input $file.  Continuing analysis.",1);
		}
		$manager->finish;
	}
	$manager->wait_all_children;
	close(FOFN);
	my @gvcfs=read_fofn($gvcf_fofn_out);
	my $gvcf_count=scalar(@gvcfs);
	logentry("COMPLETED SNP CALLING: ".$gvcf_count." g.vcf files created.",3);


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

redrep-SNPcall.pl -- De novo clustering of reduced representation libraries

=head1 SYNOPSIS

 redrep-SNPcall.pl --in FILE_OR_DIR_NAME --out DIRNAME --ref REFERENCE_FASTA [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts bam output from redrep-refmap.pl and fasta-formatted reference sequence(s) with indexes, then identifies variants.

=head1 OPTIONS

=head2 REQUIRED PARAMETERS

=over 3

=item B<-i, --in>=FILENAME

Input bam file(s).  (Required)

Can be one of the folllowing formats:

1. Path to one or more bam files (comma separated; must have ".bam" extension)
2. Path to one or more FOFN files (file with complete paths to one or more bam files -- 1 per line; must have ".fofn.list" extension)
3. Path to one or more directories of bam files (comma separated)
4. Path to one or more FODN files (file with complete paths to directories containing one or more bam files -- 1 directory per line; must have ".fodn.list" extension)
5. Comma separated list of any combination of 1-4.

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA file. (Required)

=back

=head2 OPTIONAL PARAMETERS

=head3 Program Specific Parameters

=over 3

=item B<-L, --intervals>=string

Genomic intervals to SNP call.  May consist of one or more ranges separated by a comma (e.g. -L chr1:1-100,chr2:34-500,chr3) or a list file contining 1 interval per line (must have .list extension).  (default: range of the entire reference sequence is anlayzed)

=back

=head3 Program Behavior and Resource Control

=over 3

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep>

Retain temporary intermediate files.

=item B<--mem>=integer

Max RAM usage for GATK Java Virtual Machine and some other software components in GB.  Does not necessarily guarantee maximum RAM usage (default 50)

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

Log file output path. [ Default output-dir/log.snp.txt ]

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

=item 0.0 - 1/18/2013: Draft1

=item 1.0 - 1/29/2013: Release

=item 1.2 - 7/13/2013: Dependencies Updated

=item 1.4 - 10/14/2013: GATK -dt, -dcov, -dfrac parameters added.

=item 1.5 - 3/28/2014: add version output to log; use path verisons of samtools, java; check ENV for $GATKJAR and $REDREPBIN

=item 1.6 - 9/2/2014: added --mem option to change amount of memory allocated to GATK JVM

=item 2.0 - 11/18/2016: added support for Picard Tools version 2.  Major Release.

=item 2.01 - 12/13/2016: added GATK -L parameter pass through

=item 2.1 - 1/26/2017: added GATK haplotyper w/ ERC GVCF support, g.vcf workflow, RedRep::Utils library support, code cleanup

=item 2.11 = 10/9/2019: Minor code cleanup.  Verbosity settings added.  Last version suppporting GATK v3.x.

=item 2.2 = 10/28/2019: Added GATK v4 compatibility.  Many other enhancements including parallelization, documentation, and logging.

=item 2.3 - 8/1/2020: Fixed disk space reporting bugs.

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

=item Cwd

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

=item samtools (tested with version 1.4.1)

=item java (tested with version 1.8.0_151-b12 OpenJDK Runtime Environment)

=back

=head2 Requires the following external programs be installed in the locations defined by environment variables

=over 3

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

Please acknowledge authors and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut
