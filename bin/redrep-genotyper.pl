#!/usr/bin/perl

my $ver="redrep-genotyper.pl Ver. 2.22 [01/07/2020 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils;
use RedRep::Utils qw(build_argument_list build_fofn build_listfile build_vcf_index check_genomicsdb get_timestamp read_listfile retrieve_contigs split_in_files);
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
use Time::Seconds;

use Data::Dumper;

sub split_in_files;

### ARGUMENTS WITH NO DEFAULT
my($debug,$in,$outDir,$help,$manual,$force,$keep,$version,$refFasta,$intervals,$stage,$tmpdir,$tmp_in_outdir,$no_stage_intermed,$genomics_db);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu					=	1;
my $mem						= 50;						#in GB
our $verbose			= 4;
my $gdb_batchsize	=	192;
my $javaarg="";
my $gatkarg_GenotypeGVCFs="";
my $gatkarg_GenomicsDBImport="";
my $gatkarg_CombineGVCFs="";


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
	"gatk_CombineGVCFs=s"							=>	\$gatkarg_CombineGVCFs,
	"gatk_GenotypeGVCFs=s"						=>	\$gatkarg_GenotypeGVCFs,
	"gatk_GenomicsDBImport=s"					=>	\$gatkarg_GenomicsDBImport,
	"L|intervals=s"										=>	\$intervals,
	"D|genomicsdb|db:s"								=>	\$genomics_db,
	"gdb_batchsize"										=>	\$gdb_batchsize,

	"f|force"													=>	\$force,
	"d|debug:+"												=>	\$debug,
	"k|keep|keep_temp"								=>	\$keep,
	"V|verbose:+"											=>	\$verbose,

	"v|ver|version"										=>	\$version,
	"h|help"													=>	\$help,
	"m|man|manual"										=>	\$manual
);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);


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
if($keep) {
	$mv = \&copy;
}


### DECLARE OTHER GLOBALS
my $sys;												# system call variable
my $manager = new Parallel::ForkManager($ncpu);


### THROW ERROR IF OUTPUT DIRECTORY ALREADY EXISTS (unless $force is set)
if(-d $outDir) {
	if(! $force) {
		pod2usage( -msg  => "ERROR!  Output directory $outDir already exists.  Use --force flag to overwrite.", -exitval => 2);
	}
	else {
		$sys=remove_tree($outDir);
	}
}


### CREATE OUTPUT DIR AND OPEN LOG
mkdir($outDir) || die("ERROR: Can't create output directory $outDir");
$logOut="$outDir/log.genotyper.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)",2);


### OUTPUT FILE LOCATIONS
$tmpdir								=	get_tmpdir($tmpdir);
my $intermed					=	"$tmpdir/intermed";
$intermed							=	"$outDir/intermed" if($tmp_in_outdir);
my $gvcf_fofn_out					=	"$intermed/gvcf.fofn.list";
my $interval_gvcf_dir	=	"$intermed/gvcf_intervals/";
my $interval_vcf_dir	=	"$intermed/vcf_intervals/";
my $combined_vcf_out	=	"$outDir/combined.vcf";
my $stage_dir					= $tmpdir."input/";
my $contig_list				= $intermed."/contig.list";
my $interval_list			=	$intermed."/interval.list";
my $interval_vcfs_fofn=	$intermed."/interval_vcfs.fofn.list";

mkdir($intermed) || logentry_then_die("Can't create temporary directory $intermed");
mkdir($interval_gvcf_dir) || logentry_then_die("Can't create temporary directory $interval_gvcf_dir");
mkdir($interval_vcf_dir) || logentry_then_die("Can't create temporary directory $interval_vcf_dir");
mkdir($stage_dir) || logentry_then_die("Can't create temporary directory $stage_dir");


### INITIATE LOG
logentry("Command: $0 $script\n",2);
logentry("Executing on ".hostname."\n",2);
logentry(find_job_info(),2);
logentry("Output directory $outDir (" . (df("$outDir")->{'bfree'}/1,073,741,824) . " GB free)\n",2);
logentry("Temporary directory $tmpdir (" . (df("$tmpdir")->{'bfree'}/1,073,741,824) . " GB free)\n",2);
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


############
### MAIN

	# SPLIT DIRS, FILES, FOFNs, FODNs, etc
	logentry("PARSING INPUT",3);
	my @files;
	my @vcf_indexes;

	@files=split_in_files($in, qr/\.g\.vcf$/, $gvcf_fofn_out);
	my $gvcf_count=scalar(@files);
	logentry("G.VCF FILES DETECTED: ".$gvcf_count,4);

	logentry("CHECKING/BUILDING VCF INDEX FILES",3);
	foreach my $file (@files) {
		my $new_file=build_vcf_index($file);
		push(@vcf_indexes,$new_file);
	}

	if($stage) {
		logentry("STAGING FILES TO $tmpdir",3);

		@files=stage_files($stage_dir,\@files);
		logentry(scalar(@files)." input vcf files staged to $stage_dir",4);

		@vcf_indexes=stage_files($stage_dir,\@vcf_indexes);
		logentry(scalar(@vcf_indexes)." vcf index files staged to $stage_dir",4);

	}


	# Check Reference Fasta and build index and dictionary if necessary
	check_ref_fasta($refFasta);
	my @contigs=retrieve_contigs($refFasta);


	# Define intervals.  Each interval will be a separate thread.
	my @intervals;
	if($intervals) {
		if(-e $intervals && -f $intervals) {
			@intervals=read_listfile($intervals);
		}
		else {
			@intervals=split(/,/,$intervals);
		}
	}
	else {
		@intervals=@contigs;
	}
	build_listfile(\@intervals,$interval_list);


	if($genomics_db) {
		logentry("BEGIN EXPORT TO GATK GENOMICSDB $genomics_db",3);
		my $db_out_method="--genomicsdb-workspace-path";

		# Calculate Batch Size for GenomicsDBImport
		if($gvcf_count>$gdb_batchsize) {
			$gdb_batchsize=ceil($gvcf_count/ceil($gvcf_count/$gdb_batchsize));
		}

		# If the database already exists.  Validate and set to add to database.
		if (-e $genomics_db) {
			$genomics_db=check_genomicsdb($genomics_db);
			$db_out_method="--genomicsdb-update-workspace-path";
			my $genomics_db_bak=$genomics_db."-".get_timestamp().".bak";
			logentry("GenomicsDB specified already exists.  A backup copy will be made ($genomics_db_bak).  Input gvcf files will be appended to existing database.",1);
			dircopy($genomics_db, $genomics_db_bak) || logentry_then_die("GenomicsDB $genomics_db could not be backed up to $genomics_db_bak.  File could not be created.  Please ensure sufficient disk space and write permissions.");
			symlink($genomics_db,"$outDir/$genomics_db");
		}
#		elsif($genomics_db eq "") {
#			logentry("GenomicsDB name not specified.  Using default ($outDir/redrep.gatk.db).",1);
#			$genomics_db="$outDir/redrep.gatk.db";
#		}
		else {
			$genomics_db="$outDir/$genomics_db";
		}
		$sys=cmd("$gatk GenomicsDBImport $db_out_method $genomics_db --variant $gvcf_fofn_out --reference $refFasta --intervals $interval_list --batch-size $gdb_batchsize $gatkarg_GenomicsDBImport","Building GenomicsDB: $genomics_db");
		logentry("FINISH EXPORT TO GATK GENOMICSDB $genomics_db",3);
		$genomics_db="gendb://".$genomics_db;
	}


	my $gvcf_files = build_argument_list(\@files,"--variant");

	logentry("BEGIN GENOTYPING",3);
	foreach my $interval(@intervals) {
		$manager->start and next;
		my $geno_in;
		if($genomics_db) {
			$geno_in=$genomics_db;
		}
		else {
			logentry("Beginning CombineGVCFs on interval $interval",4);
			my $interval_gvcf_out="$interval_gvcf_dir/interval-$interval.g.vcf";
			$sys=cmd("$gatk CombineGVCFs --reference $refFasta --intervals $interval $gvcf_files --output $interval_gvcf_out --tmp-dir $tmpdir $gatkarg_CombineGVCFs","Calling GATK CombineGCVFs");
			logentry("Completed CombineGVCFs on interval $interval",4);
			$geno_in=$interval_gvcf_out;
		}

		my $interval_vcf_out="$interval_vcf_dir/interval-$interval.vcf";
		logentry("Beginning GenotypeGVCFs on interval $interval",4);
		$sys=cmd("$gatk GenotypeGVCFs --reference $refFasta --intervals $interval --variant $geno_in --output $interval_vcf_out --tmp-dir $tmpdir $gatkarg_GenotypeGVCFs","Calling GATK GenotypeGCVFs");
		logentry("Completed GenotypeGVCFs on interval $interval",4);
		$manager->finish;
	}
	$manager->wait_all_children;
	logentry("COMPLETE GENOTYPING",3);


	# MAKE INTERVAL VCF OUTPUT FOFN
	my @interval_vcf_outs=get_dir_filelist($interval_vcf_dir,qr/\.vcf$/);
	build_fofn(\@interval_vcf_outs,$interval_vcfs_fofn);


	# MERGE VCFs
	logentry("MERGING VCFs",3);
	$sys=cmd("$gatk MergeVcfs --INPUT $interval_vcfs_fofn --OUTPUT $combined_vcf_out --TMP_DIR $tmpdir","Calling GATK MergeVcfs");
	logentry("VCF MERGE COMPLETE",3);


	# FILE CLEAN UP
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

redrep-genotyper.pl -- Combine gvcf files and call genotypes

=head1 SYNOPSIS

 redrep-genotyper.pl --in FILE_OR_DIR_NAME --out DIRNAME [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts bam output from redrep-refmap.pl and fasta-formatted reference sequence(s) with indexes, then identifies variants.

=head1 OPTIONS

=head2 REQUIRED PARAMETERS

=over 3

=item B<-i, --in>=FILENAME

Input GVCF(s).  (Required)

Can be one of the folllowing formats:

1. Path to one or more g.vcf files (comma separated; must have ".g.vcf" extension)
2. Path to one or more FOFN files (file with complete paths to one or more g.vcf files -- 1 per line; must have ".fofn" or ".fofn.list" extension)
3. Path to one or more directories of g.vcf files (comma separated)
4. Path to one or more FODN files (file with complete paths to directories containing one or more g.vcf files -- 1 directory per line; must have ".fodn" or ".fodn.list" extension)
5. Comma separated list of any combination of 1-4.

NOTE: Additional (or alternative) inputs can be included in an existing GATK genomicsDB (see --genomicdb)

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA. (Required)

=back

=head2 OPTIONAL PARAMETERS

=head3 Program Specific Parameters

=over 3

=item B<-D,--genomicsdb,--db>=DIRECTORY_NAME

If specified inputs will be imported into a GATK GenomicsDB before genotyping (typically speeds analysis significantly, but may fail with large numbers of inputs or intervals/contigs.).  If the directory does not exist, a new GenomicsDB will be created.  If the directory exists and is a valid genomicsDB, then inputs will be added to the existing database.

=item B<--gdb_batchsize>=integer

If more than this number of samples are being imported into GATK GenomicsDB at once, then it will be split into batches of not more than this number.  For instance, if set at 192: 194 samples would be processed in batches of 97; 384 would be processed in batches of 192.  Reduces memory usage, at the cost of performance.  Only has effect in combination with --genomicsdb.  (default=192)

=item B<-L, --intervals>=string

Genomic intervals to genotype.  May consist of one or more ranges separated by a comma (e.g. -L chr1:1-100,chr2:34-500,chr3) or a list file contining 1 interval per line (must have .list extension).  (default: range of the entire reference sequence is anlayzed)

=back

=head3 Program Behavior and Resource Control

=over 3

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep>

Retain temporary intermediate files.

=item B<--mem>=integer

Max RAM usage for GATK Java Virtual Machine in GB (default 50).  Each thread will be allocated this much memory.

=item B<-S, --stage>

When specified, input files will be staged to the tmpdir.  Can increase performance on clusters and other situations where the output directory is on network attached or cloud storage.

=item B<-t, --threads>=integer

Number of cpu's to use for threadable operations.  WARNING: EACH THREAD will require the amount of RAM specified in --mem (default 50GB).  Using too many threads may lead to crashes due to excessive RAM usage or simultaneous files open.

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

=item 2.11 - 10/9/2019: Minor code cleanup.  Verbosity settings added.  Last version supporting GATK v3.x.

=item 2.2 - 10/28/2019: Added GATK v4 compatibility.  Added genomicsDB compatibility and many other enhancements including parallelization, documentation, and logging.

=back

=head1 DEPENDENCIES

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

Copyright 2012-2019 Shawn Polson, Randall Wisser, Keith Hopper.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Please acknowledge authors and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut
