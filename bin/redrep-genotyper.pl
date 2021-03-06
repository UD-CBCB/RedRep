#!/usr/bin/perl

my $ver="redrep-genotyper.pl Ver. 2.3 [08/13/2020 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils;
use RedRep::Utils qw(archive_gdb build_argument_list build_fofn build_listfile build_vcf_index check_genomicsdb gdb_history get_timestamp read_listfile retrieve_contigs split_in_files);
use Getopt::Long qw(:config no_ignore_case);
use Carp qw(cluck confess longmess);
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
our($outDir, $intermed, $tmpdir, $debug, $no_stage_intermed, $tmp_in_outdir); 	# globals to properly handle file unstaging in case of fatal errors
my($in,$help,$manual,$force,$version,$refFasta,$intervals,$stage,$no_stage_intermed,$savetmp,$genomics_db,$genomics_db_final,$no_geno,$gvcf_files,$db_load_args);
our($keep);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
our $ncpu					=	1;
my $mem						= 50;						#in GB
our $verbose			= 4;
my $gdb_batchsize	=	192;
my $javaarg="";
my $gatkarg_GenotypeGVCFs="";
my $gatkarg_GenomicsDBImport="";
my $gatkarg_CombineGVCFs="";
my $gdb_append=0;


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
	"D|genomicsdb|db=s"								=>	\$genomics_db,
	"gdb_batchsize=i"									=>	\$gdb_batchsize,
	"db_load_only"										=>	\$no_geno,
	"db_load_args=s"									=>	\$db_load_args,

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
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory/db) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);
pod2usage( -msg  => "ERROR!  Option --genomicsdb must be specified when --db_load_only is used.\n", -exitval => 2)  if ($no_geno && ! $genomics_db);


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
mkdir($outDir) || confess("ERROR: Can't create output directory $outDir");
$logOut="$outDir/log.genotyper.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)",2);


### OUTPUT FILE LOCATIONS
$tmpdir								=	get_tmpdir($tmpdir);
$intermed							=	"$tmpdir/intermed";
$intermed							=	"$outDir/intermed" if($tmp_in_outdir);
my $gvcf_fofn_out			=	"$intermed/gvcf.fofn.list";
my $interval_gvcf_dir	=	"$intermed/gvcf_intervals/";
my $interval_vcf_dir	=	"$intermed/vcf_intervals/";
my $combined_vcf_out	=	"$outDir/combined.vcf";
my $stage_dir					= $tmpdir."input/";
my $contig_list				= $intermed."/contig.list";
my $interval_list			=	$intermed."/interval.list";
my $interval_vcfs_fofn=	$intermed."/interval_vcfs.fofn.list";

if($genomics_db) {
	$genomics_db =~ s!^gendb://!!;					# remove gendb prefix if present
	$genomics_db =~ s!/$!!;						     	# remove trailing slash if present
	# If doesn't contain a path and if it is not referring to an existing subdir in working dir
	if($genomics_db !~ m!/! && ! -d $genomics_db) {
		$genomics_db = "${outDir}/${genomics_db}";    # place in outDir
	}
}

mkdir($intermed) || logentry_then_die("Can't create temporary directory $intermed");
mkdir($interval_gvcf_dir) || logentry_then_die("Can't create temporary directory $interval_gvcf_dir");
mkdir($interval_vcf_dir) || logentry_then_die("Can't create temporary directory $interval_vcf_dir");
mkdir($stage_dir) || logentry_then_die("Can't create temporary directory $stage_dir");


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
	my $java_opts = "-Xmx${mem}g $javaarg ";

	# samtools
	my $samtools = get_path_samtools(1);

	#GATK
	my $gatk  = get_path_gatk(1);
	$gatk .= qq( --java-options "$java_opts");

	#DETECT IF TILEDB FILE LOCKING IS DISABLED
	my $tile_init_state;
	$tile_init_state=$ENV{'TILEDB_DISABLE_FILE_LOCKING'} if ($ENV{'TILEDB_DISABLE_FILE_LOCKING'});

############
### MAIN

	# SPLIT DIRS, FILES, FOFNs, FODNs, etc
	logentry("PARSING INPUT",3);

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

	my @files;
	my @files_orig;
	my @vcf_indexes;
	my $genomics_db_in;
	my $gvcf_count;

	# GENOMICSDB INPUT
	if($in =~ m{^gendb://(.+)$}) {
		my $in_dir=$1;
		#push(@files, $in);

		logentry_then_die("GenomicsDB $genomics_db specified as input does not exist.") if (! -e $in_dir);
		if($stage) {
			logentry("STAGING FILES TO $tmpdir",3);
			$in_dir=stage_gdb($in_dir,$tmpdir);
			logentry("GenomicsDB $in staged to $in_dir (for read operations only).",4);
			logentry("Working copy of GenomicsDB does not appear to be on a shared file system.  Setting environment variable TILEDB_DISABLE_FILE_LOCKING=0",4);
			$ENV{'TILEDB_DISABLE_FILE_LOCKING'}=0;
		}
		#$genomics_db="gendb://";
		$genomics_db=check_genomicsdb($in_dir,1);
		$genomics_db_in=1;
	}
	#FILE INPUT
	else {
		@files=split_in_files($in, qr/\.g\.vcf$/, $gvcf_fofn_out);
		@files_orig=@files;
		$gvcf_count=scalar(@files);
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

		$gvcf_files = build_argument_list(\@files,"--variant");

		if($genomics_db) {
			logentry("SETTING WORKING COPY OF GENOMICS_DB",3);
			$genomics_db_final=$genomics_db;		# keep this version as final output location
			$genomics_db =~ m{/?([^/]+)/?$};
			$genomics_db=$intermed."/".$1;

			#If already existing database, copy to working location
			if(-e $genomics_db_final) {
				logentry("Existing GenomicsDB ($genomics_db_final) being staged to ($genomics_db) as working copy",4);
				dircopy($genomics_db_final,$genomics_db);
			} else {
				logentry("No GenomicsDB found at $genomics_db_final, new genomicsDB being staged to ($genomics_db) as working copy",4);
			}

			if($tmp_in_outdir) {
				logentry("WARNING: When using --tmp_in_outdir working copy of the genomicsDB will be stored in $outDir, this may cause problems if $outDir is on a shared file system (e.g. NFS).",2);
				logentry("Working copy of GenomicsDB may be on a shared file system.  Setting environment variable TILEDB_DISABLE_FILE_LOCKING=1.",4);
				$ENV{'TILEDB_DISABLE_FILE_LOCKING'}=1;
			} else {
				logentry("Working copy of GenomicsDB does not appear to be on a shared file system.  Setting environment variable TILEDB_DISABLE_FILE_LOCKING=0.",4);
				$ENV{'TILEDB_DISABLE_FILE_LOCKING'}=0;
			}
		}
	}

	#LOAD NEW DATA TO GENOMICSDB
	if($genomics_db && ! $genomics_db_in) {
		logentry("BEGIN EXPORT TO GATK GENOMICSDB $genomics_db",3);
		my $db_out_method="--genomicsdb-workspace-path";

		# Calculate Batch Size for GenomicsDBImport
		if($gvcf_count>$gdb_batchsize) {
			$gdb_batchsize=ceil($gvcf_count/ceil($gvcf_count/$gdb_batchsize));
		}

		# If the database already exists.  Validate and set to add to database.
		if (-e $genomics_db) {
			$genomics_db=check_genomicsdb(${genomics_db});
			$db_out_method="--genomicsdb-update-workspace-path";
		}
		$sys=cmd("$gatk GenomicsDBImport $db_out_method $genomics_db $db_load_args --variant $gvcf_fofn_out --reference $refFasta --intervals $interval_list --batch-size $gdb_batchsize $gatkarg_GenomicsDBImport","Building GenomicsDB: $genomics_db");
		logentry("FINISH EXPORT TO GATK GENOMICSDB $genomics_db",3);
#		$genomics_db="gendb://".$genomics_db;
	}



	unless($no_geno) {
		logentry("BEGIN GENOTYPING",3);
		foreach my $interval(@intervals) {
			$manager->start and next;
			my $geno_in;
			if($genomics_db) {
				$geno_in="gendb://".$genomics_db;
			}	else {
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
	}	else {
		logentry("SKIPPING GENOTYPING: option --no_geno specified",3);
	}

	# FILE CLEAN UP
	logentry("FILE CLEAN UP",3);
	if($keep && ! $tmp_in_outdir) {
		logentry("Saving intermediate tmp directory",4);
		eval { $sys=dircopy($intermed, "$outDir/intermed") } or do {
			$savetmp=1;
			logentry("Could not copy intermediate file directory ($intermed) to final output directory ($outDir) as requested with --keep option.  Retaining temporary files after execution.",1);
		};
	}
	my $gdb_archived=0;
	# If a genomicsDB was modified
	if($genomics_db && ! $genomics_db_in) {
		# If genomicsDB was in a temporary location move to the final location
		if($genomics_db_final) {
			# If final location already exists make an archive copy first
			if(-e $genomics_db_final) {
				$gdb_archived=1;
				logentry("Archiving original genomicsDB",4);
				eval { archive_gdb($genomics_db_final) } or do {
					logentry("Could not archive original genomicsDB $genomics_db_final.",1);
					until(! -e $genomics_db_final) {
						$genomics_db_final="${genomics_db_final}.new";
					}
					logentry("Saving genomicsDB to alternate location $genomics_db_final.",1);
				};
			}
			logentry("Copying genomicsDB to output location",4);
			eval { dircopy($genomics_db,$genomics_db_final) } or do {
				$savetmp=1;
				logentry("Could not move working copy of genomicsDB ($genomics_db) to final location ($genomics_db_final). Retaining temporary files after execution.",1);
			};
		}	else {
			$genomics_db_final=$genomics_db;
		}
		gdb_history($genomics_db_final,$logOut,$gdb_archived,\@files_orig);
		#If final location for db is not in the output, make symbolic links
		if($genomics_db_final !~ /$outDir/) {
			$genomics_db_final =~ m{/?([^/]+)/?$};
			symlink($genomics_db_final,"$outDir/$1");
			copy("${genomics_db_final}.history","${outDir}/$1.history");
		}
	}

	unless($savetmp) {
		logentry("Removing tmp files",4);
		$sys=remove_tree($tmpdir);
	}
	else {
		logentry("Temporary files retained in ${tmpdir} on execution node (".hostname.") to prevent data loss.  See previous warning(s) for cause.",1);
	}

	# RETURN ENVIRONMENT VARIABLES TO INITIAL STATE
	if($tile_init_state) {
		$ENV{'TILEDB_DISABLE_FILE_LOCKING'}=$tile_init_state;
	}
	elsif($ENV{'TILEDB_DISABLE_FILE_LOCKING'}) {
		delete $ENV{'TILEDB_DISABLE_FILE_LOCKING'};
	}


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

 redrep-genotyper.pl --in FILE_OR_DIR_NAME --out DIRNAME --ref REFERENCE_FASTA [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts bam output from redrep-refmap.pl and fasta-formatted reference sequence(s) with indexes, then identifies variants.

=head1 OPTIONS

=head2 REQUIRED PARAMETERS

=over 3

=item B<-i, --in>=FILENAME

Input GVCF(s) or GenomicsDB.  (Required)

Can be one of the folllowing formats:

1. Path to one or more g.vcf files (comma separated; must have ".g.vcf" extension)
2. Path to one or more FOFN files (file with complete paths to one or more g.vcf files -- 1 per line; must have ".fofn" or ".fofn.list" extension)
3. Path to one or more directories of g.vcf files (comma separated)
4. Path to one or more FODN files (file with complete paths to directories containing one or more g.vcf files -- 1 directory per line; must have ".fodn" or ".fodn.list" extension)
5. Comma separated list of any combination of 1-4.
6. An existing GenomicsDB (must include the prefix gendb://)

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

If specified inputs will be imported into a GATK GenomicsDB before genotyping (typically speeds analysis significantly, but may fail with large numbers of inputs or intervals/contigs.).  If the directory does not exist, a new GenomicsDB will be created.  If the directory exists and is a valid genomicsDB, then inputs will be added to the existing database.  Parameter has no effect if the --in is a genomicsDB.

=item B<--db_load_only>

Use to load data to genomics db only (genotyping will be skipped).  If specified the parameter --genomicsdb is also required.

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

=item 2.3 - 8/7/2020: Fixed disk space reporting and gdb_batchsize bugs.  Added --db_load_only.  Added genomicsDB batching and staging capabilities.

=back

=head1 DEPENDENCIES

=head2 Requires the following external programs be in the system PATH:

=over 3

=item samtools (tested with version 1.4.1)

=item java (tested with version 1.8.0_151-b12 OpenJDK Runtime Environment)

=back

=head2 Requires the following external programs be installed in the locations defined by environment variables

=over 3

=item GATK (tested with version 4.1.8.1)

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
