#!/usr/bin/perl

my $ver="redrep-genotyper.pl Ver. 2.2beta [10/11/2019 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils qw(check_dependency check_jar cmd dir_get_files find_job_info find_tmpdir logentry logentry_then_die);
use Getopt::Long qw(:config no_ignore_case);
#use Parallel::ForkManager;
use Pod::Usage;
use POSIX qw(floor);
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);
use Filesys::Df qw(df);
use Sys::Hostname qw(hostname);

sub split_in_files;

### ARGUMENTS WITH NO DEFAULT
my($debug,$in,$outDir,$help,$manual,$force,$keep,$version,$refFasta,$intervals,$javaarg,$gatkarg,$stage);


### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu			=	1;
my $mem				= 50;						#in GB
our $verbose	= 0;
our $tmpdir		=	find_tmpdir();

GetOptions (
	"i|in=s"													=>	\$in,
	"o|out=s"													=>	\$outDir,
	"r|ref=s"													=>	\$refFasta,
	"l|log=s"													=>	\$logOut,

	"t|threads|ncpu=i"								=>	\$ncpu,
	"mem=i"														=>	\$mem,
	"T|tmpdir=s"											=>	\$tmpdir,
	"S|stage"													=>	\$stage,

	"java=s"													=>	\$javaarg,
	"gatk=s"													=>	\$gatkarg,
	"L|intervals=s"										=>	\$intervals,

	"f|force"													=>	\$force,
	"d|debug"													=>	\$debug,
	"k|keep|keep_temp"								=>	\$keep,
	"V|verbose+"											=>	\$verbose,

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
	require Data::Dumper; import Data::Dumper;
	require diagnostics; import diagnostics;
	$keep=1;
	$verbose=10;
}


### SET DEFAULT METHOD OF FILE PROPAGATION
my $mv = \&move;
if($keep) {
	$mv = \&copy;
}


### DECLARE OTHER GLOBALS
my $sys;												# system call variable
#my $manager = new Parallel::ForkManager( floor($ncpu/$hmm_threads) );
my $GATKargs=" ";
#my $gvcf_files;

### THROW ERROR IF OUTPUT DIRECTORY ALREADY EXISTS (unless $force is set)
if(-d $outDir) {
	if(! $force) {
		pod2usage( -msg  => "ERROR!  Output directory $outDir already exists.  Use --force flag to overwrite.", -exitval => 2);
	}
	else {
		$sys=remove_tree($outDir);
	}
}


### OUTPUT FILE LOCATIONS
my $intermed="$outDir/intermed";
my $fofn_out="$intermed/gvcf.fofn.list";
my $combined_gvcf_out="$intermed/combined.g.vcf";


### CREATE OUTPUT DIR
mkdir($outDir);
mkdir($outDir."/intermed");


### CREATE LOG FILES
$logOut="$outDir/log.geno.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)",0);
print $LOG "Command: $0 $script\n";
print $LOG "Executing on ".hostname."\n";
print $LOG find_job_info();
print $LOG "Temporary directory $tmpdir (" . (df("$tmpdir")->{'bfree'}/1,073,741,824) . " GB free)\n";
print $LOG "Running in Debug Mode\n" if($debug && $debug>0);
print $LOG "Keeping intermediate files.  WARNING: Can consume significant extra disk space\n" if($keep && $keep>0);
print $LOG "Log verbosity level: $verbose\n" if($verbose && $verbose>0);


### REDREP BIN/SCRIPT LOCATIONS
my $execDir=$ENV{'REDREPBIN'};
my $utilDir=$ENV{'REDREPUTIL'};
my $libDir=$ENV{'REDREPLIB'};
print $LOG "RedRep Bin: $execDir\n";
print $LOG "RedRep Utilities: $utilDir\n";
print $LOG "RedRep Libraries: $libDir\n";


### CHECK FOR AND SETUP EXTERNAL DEPENDENCIES
logentry("Checking External Dependencies",0);

	# samtools
	my $samtools=check_dependency("samtools","--version","s/\r?\n/ | /g");

	# java
	#our $java=check_dependency("java","-version","s/\r?\n/ | /g");    #not currently needed; retain for future
	my $java_opts = "-Xmx${mem}g $javaarg -d64";

	#GATK
	#	my $GATK=check_jar("GATK",$ENV{'GATKJAR'},"--version","s/\r?\n/ | /g","Environment Variable GATKJAR is not defined or is a not pointing to a valid jar file!  Please define valid GATK jar file location with command: export GATKJAR='PATH_TO_GATK_JAR");
	my $gatk  = check_dependency("gatk","--version 2>&1 | tail -n +4","s/\r?\n/ | /g");
	$gatk .= qq( --java-options "$java_opts");



############
### MAIN

	# SPLIT DIR, FILES, FOFN, FODN, etc
	logentry("PARSING INPUT",0);
	my @files=split_in_files($in);
	logentry("G.VCF FILES DETECTED: ".scalar(@files),1);

	if($stage) {
		logentry("STAGING INPUT FILES TO $tmpdir",0);
		my @tmpfiles=@files;
		@files=stage_files(@tmpfiles);
		logentry(scalar(@files)." files staged to $tmpdir",1);
	}

	# ref Fasta stub
	my $refFasta_stub=fileparse($refFasta, qr/\.[^.]*$/);

	if(! -e $refFasta.".fai") {
		logentry("REFERENCE FASTA INDEX NOT FOUND: BUILDING",0);
		$sys=cmd("$samtools faidx $refFasta","Building reference index");
	}
	if(! -e $refFasta_stub.".dict") {
		logentry("REFERENCE DICTIONARY NOT FOUND: BUILDING",0);
		$sys=cmd("$gatk CreateSequenceDictionary --REFERENCE $refFasta --OUTPUT $refFasta_stub.dict","Building reference dictionary");
	}


	if($intervals) {
		my @intervals=split(/,/,$intervals);
		foreach my $interval (@intervals) {
			$GATKargs.="-L $interval "
		}
	}

	$GATKargs.="$gatkarg " if($gatkarg);

	my $gvcf_files = "--variant ";
	$gvcf_files .= join(" --variant ",@files);
	logentry("BEGIN GATK CombineGCVFs",0);
	$sys=cmd("$gatk CombineGVCFs --reference $refFasta $gvcf_files --output $combined_gvcf_out --tmp-dir $tmpdir $GATKargs","Calling GATK CombineGCVFs");
	logentry("COMPLETE GATK CombineGVCFs",0);

	logentry("BEGIN GATK GenotypeGVCFs",0);
	$sys=cmd("$gatk GenotypeGVCFs --reference $refFasta --variant $combined_gvcf_out --output $outDir/combined.vcf --tmp-dir $tmpdir $GATKargs","Calling GATK GenotypeGCVFs");
	logentry("COMPLETE GATK GenotypeGVCFs",0);



	# FILE CLEAN UP
	logentry("FILE CLEAN UP",0);
	if(! $keep) {
		logentry("Remove intermediate file directory",0);
		$sys=remove_tree($intermed);
	}


logentry("SCRIPT COMPLETE",0);
close($LOG);

exit 0;


#######################################
############### SUBS ##################

#######################################
### split_in_files
# Split input files, directories, FOFNs, and/or FODNs
sub split_in_files {
	my $in=shift;
	my @in=split(/,/, $in);
	foreach my $innie(@in) {
		chomp;
		# DIRECTORY
		if(-d $innie) {
			logentry("Parsing g.vcf files from directory $innie",1);
			my @tmpfiles=dir_get_files($innie,qr/\.g\.vcf/);
			push(@files,@tmpfiles);
		}
		# FODN
		elsif (-f $innie && ($innie=~/.fodn$/i || $innie=~/.fodn.list$/i)) {
			logentry("Parsing directories from FODN $innie",1);
			open(my $FODN,"<",$innie);
			while(<$FODN>) {
				chomp;
				my $dirname=$_;
				if(-d $dirname) {
					logentry("Parsing g.vcf files from directory $innie",2);
					my @tmpfiles=dir_get_files($innie,qr/\.g\.vcf/);
					push(@files,@tmpfiles);
				}
				else {
					logentry_then_die("ERROR: Input FODN $innie contains an entry that is not a directory: $dirname");
				}
			}
			close($FODN);
		}
		# FOFN
		elsif (-f $innie && ($innie=~/.fofn$/i || $innie=~/.fofn.list$/i)) {
			logentry("Parsing g.vcf files from FOFN $innie",1);
			open(my $FOFN,"<",$innie);
			while(<$FOFN>) {
				chomp;
				if (-f $_) {
					push(@files,$_);
				}
				else {
					logentry_then_die("ERROR: File ($_) specified in input FOFN $innie does not exist!");
				}
			}
			close($FOFN);
		}
		# FILENAME
		elsif (-f $innie && $innie=~/.g.vcf$/i) {
			push(@files,$innie);
		}
		# ERROR
		else {
			pod2usage( -msg  => "ERROR!  Argument to -i ($innie), file/directory not found or if a file does not have a proper file extension: .g.vcf, .fofn, .fodn.\n", -exitval => 2) if (! $innie);
		}
		logentry("Writing combined g.vcf file paths to FOFN $fofn_out",0);
		open(my $FOFN_OUT, ">", "$fofn_out");
		foreach (@files) {
			print $FOFN_OUT "$_\n";
		}
		close($FOFN_OUT);
	}
	return(@files);
}




__END__

=pod

=head1 NAME

redrep-genotyper.pl -- Call genotypes from GVCF files

=head1 SYNOPSIS

 redrep-genotyper.pl --in FILE_OR_DIR_NAME --out DIRNAME [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts bam output from redrep-refmap.pl and fasta-formatted reference sequence(s) with indexes, then identifies variants.

=head1 OPTIONS

=over 3

=item B<-1, -i, --in, --in1>=FILENAME

Input GVCF's.  (Required)

Can be one of the folllowing formats:

1. Path to one or more g.vcf files (comma separated; must have ".g.vcf" extension)
2. Path to one or more FOFN files (file with complete paths to one or more g.vcf files -- 1 per line; must have ".fofn" or ".fofn.list" extension)
3. Path to one or more directories of g.vcf files (comma separated)
4. Path to one or more FODN files (file with complete paths to directories containing one or more g.vcf files -- 1 directory per line; must have ".fodn" or ".fodn.list" extension)
5. Comma separated list of any combination of 1-4.

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA with index and dict files. (Required)

=item B<-L, --intervals>=string

Passes argument of same name to GATK.  May consist of one or more ranges separated by a comma. (EXAMPLE: -L chr1:1-100,chr2:34-500,chr3).  GATK description:

Use this option to perform the analysis over only part of the genome. You can use samtools-style intervals either explicitly on the command line (e.g. -L chr1 or -L chr1:100-200) or by loading in a file containing a list of intervals (e.g. -L myFile.intervals). Additionally, you can also specify a ROD file (such as a VCF file) in order to perform the analysis at specific positions based on the records present in the file (e.g. -L file.vcf). Finally, you can also use this to perform the analysis on the reads that are completely unmapped in the BAM file (i.e. those without a reference contig) by specifying -L unmapped.

=item B<-l, --log>=FILENAME

Log file output path. [ Default output-dir/log.snp.txt ]

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep>

Retain temporary intermediate files.

=item B<--mem>=integer

Max RAM usage for GATK Java Virtual Machine in GB (default 50)


=item B<-T, --tmpdir>

Set directory (local directory recommended) for fast temporary and staged file operations.  Defaults to system environment variable REDREP_TMPDIR, TMPDIR, TEMP, or TMP in that order (if they exist).  Otherwise assumes /tmp.

=item B<-V, --verbose>

Produce detailed log.  Can be involed multiple times for additional detail levels.

=item B<-v, --ver, --version>

Displays the current version.

=item B<-h, --help>

Displays the usage message.

=item B<-m, --man, --manual>

Displays full manual.

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

=item 2.11 = 10/9/2019: Minor code cleanup.  Verbosity settings added.

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

=item GATK (tested with version 3.8-0); Location defined by `env GATKJARS`

=item picard-tools (tested with version 2.19.0 SNAPSHOT); Location defined by `env PICARDJAR`

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
