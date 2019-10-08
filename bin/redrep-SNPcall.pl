#!/usr/bin/perl

# MANUAL FOR redrep-SNPcall.pl

=pod

=head1 NAME

redrep-SNPcall.pl -- De novo clustering of reduced representation libraries

=head1 SYNOPSIS

 redrep-SNPcall.pl --in FILE_OR_DIR_NAME --out DIRNAME [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts bam output from redrep-refmap.pl and fasta-formatted reference sequence(s) with indexes, then identifies variants.

=head1 OPTIONS

=over 3

=item B<-1, -i, --in, --in1>=FILENAME

Input file in sorted bam format with index or directory of such files. (Required)

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-r, --ref>=FILENAME

Reference FASTA with index and dict files. (Required)

=item B<-l, --log>=FILENAME

Log file output path. [ Default output-dir/log.snp.txt ]

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep_temp>

Retain temporary intermediate files.

=item B<-t, --threads, --ncpu>=integer

Number of cpu's to use for threadable operations.

=item B<--mem>=integer

Max RAM usage for GATK Java Virtual Machine in GB (default 50)

=item B<-dcov, --downsample_to_coverage>=integer

Passes argument of same name to GATK SNP caller.  See also -dfrac, -dt.  GATK description:

Coverage [integer] to downsample to. For locus-based traversals (eg., LocusWalkers and ActiveRegionWalkers),this controls the maximum depth of coverage at each locus. For non-locus-based traversals (eg., ReadWalkers), this controls the maximum number of reads sharing the same alignment start position. Note that this downsampling option does NOT produce an unbiased random sampling from all available reads at each locus: instead, the primary goal of the to-coverage downsampler is to maintain an even representation of reads from all alignment start positions when removing excess coverage. For a true across-the-board unbiased random sampling of reads, use -dfrac instead. Also note that the coverage target is an approximate goal that is not guaranteed to be met exactly: the downsampling algorithm will under some circumstances retain slightly more coverage than requested.

=item B<-dfrac, --downsample_to_fraction>=integer

Passes argument of same name to GATK SNP caller.  See also -dcov, -dt.  GATK description:

Fraction [0.0-1.0] of reads to downsample to.

=item B<-dt, --downsampling_type>=NONE, ALL_READS, BY_SAMPLE

Passes argument of same name to GATK SNP caller.  See also -dcov, -dfrac.  GATK description:

Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from the pile based on the method described here.

=item B<-L, --intervals>=string

Passes argument of same name to GATK SNP caller.  May consist of one or more ranges separated by a comma. (EXAMPLE: -L chr1:1-100,chr2:34-500,chr3).  GATK description:

Use this option to perform the analysis over only part of the genome. You can use samtools-style intervals either explicitly on the command line (e.g. -L chr1 or -L chr1:100-200) or by loading in a file containing a list of intervals (e.g. -L myFile.intervals). Additionally, you can also specify a ROD file (such as a VCF file) in order to perform the analysis at specific positions based on the records present in the file (e.g. -L file.vcf). Finally, you can also use this to perform the analysis on the reads that are completely unmapped in the BAM file (i.e. those without a reference contig) by specifying -L unmapped.

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

=item 0.0 - 1/18/2013: Draft1

=item 1.0 - 1/29/2013: Release

=item 1.2 - 7/13/2013: Dependencies Updated

=item 1.4 - 10/14/2013: GATK -dt, -dcov, -dfrac parameters added.

=item 1.5 - 3/28/2014: add version output to log; use path verisons of samtools, java; check ENV for $GATKJAR and $REDREPBIN

=item 1.6 - 9/2/2014: added --mem option to change amount of memory allocated to GATK JVM

=item 2.0 - 11/18/2016: added support for Picard Tools version 2.  Major Release.

=item 2.01 - 12/13/2016: added GATK -L parameter pass through

=item 2.1 - 1/26/2017: added GATK haplotyper w/ ERC GVCF support, g.vcf workflow, RedRep::Utils library support, code cleanup

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

Copyright 2012-2017 Shawn Polson, Randall Wisser, Keith Hopper.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Please acknowledge authors and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

my $script=join(' ',@ARGV);

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);
use Sys::Hostname qw(hostname);

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils qw(check_dependency check_jar cmd find_job_info get_execDir logentry logentry_then_die);


### ARGUMENTS WITH NO DEFAULT
my($in,$outDir,$help,$manual,$force,$metaFile,$keep,$version,$refFasta,$dcov,$dfrac,$dt,$intervals,$javaarg,$gatkarg);
our($debug);

### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu		=	1;
my $mem			= 50;						#in GB

GetOptions (	"i|in=s"										=>	\$in,
				"o|out=s"													=>	\$outDir,
				"r|ref=s"													=>	\$refFasta,
				"l|log=s"													=>	\$logOut,

				"t|threads|ncpu=i"								=>	\$ncpu,
				"mem=i"														=>	\$mem,

				"java=s"													=>	\$javaarg,
				"gatk=s"													=>	\$gatkarg,
				"dcov|downsample_to_coverage=i"		=>	\$dcov,
				"dfrac|downsample_to_fraction=f"	=>	\$dfrac,
				"dt|downsampling_type=s"					=>	\$dt,
				"L|intervals=s"										=>	\$intervals,

				"f|force"													=>	\$force,
				"d|debug"													=>	\$debug,
				"k|keep|keep_temp"								=>	\$keep,

				"v|ver|version"										=>	\$version,
				"h|help"													=>	\$help,
				"m|man|manual"										=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
my $ver="redrep-SNPcall.pl Ver. 2.2 (1/26/2017 rev)";
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);
pod2usage( -msg  => "ERROR!  Arguments -dcov and -dfrac are incompatible, chose one.\n", -exitval => 2)  if ($dcov && $dfrac);


### DEBUG MODE
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
mkdir($outDir."/gvcf");
my(@files);
if(-d $in)
{	opendir(DIR,$in);
	@files=grep /\.bam$/, readdir(DIR);
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
$logOut="$outDir/log.SNP.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)");
print $LOG "Command: $0 $script\n";
print $LOG "Executing on ".hostname."\n";
print $LOG find_job_info()."\n";


### REDREP BIN/SCRIPT LOCATIONS
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

	# samtools
	my $samtools=check_dependency("samtools","--version","s/\r?\n/ | /g");

	#picard
	#picard
	# Picard tools (as of v 2.19.0) gives an error code of 1 when version number is checked.  To get around the $? check in &cmd(), added " 2>&1;v=0" to the version check flag to make error code 0, but still get version redirected to SDTOUT
	my $picard=check_jar("Picard Tools",$ENV{'PICARDJAR'},"AddCommentsToBam --version 2>&1;v=0","s/\r?\n/ | /g","Environment Variable PICARDJAR is not defined or is a not pointing to a valid file!  Please define valid picard tools location with command: export PICARDJAR='PATH_TO_PICARD_JAR");

	#GATK
	my $GATK=check_jar("GATK",$ENV{'GATKJAR'},"--version","s/\r?\n/ | /g","Environment Variable GATKJAR is not defined or is a not pointing to a valid jar file!  Please define valid GATK jar file location with command: export GATKJAR='PATH_TO_GATK_JAR");


### OUTPUT FILE LOCATIONS
#my $gatk_out="$outDir/combined.SNP.vcf";


############
### MAIN

	# ref Fasta stub
	my $stub2=fileparse($refFasta, qr/\.[^.]*$/);

	if(! -e $refFasta.".fai")
	{	logentry("REFERENCE FASTA INDEX NOT FOUND: BUILDING");
		$sys=cmd("$samtools faidx $refFasta","Building reference index");
	}
	if(! -e $stub2.".dict")
	{	logentry("REFERENCE DICTIONARY NOT FOUND: BUILDING");
		$sys=cmd("$java $picard CreateSequenceDictionary R=$refFasta O=$stub2.dict","Building reference dictionary");
	}

	my $bam_files=join(" -I ",@files);
	my $GATKargs.="-dcov $dcov " if($dcov);
	$GATKargs.="-dfrac $dfrac " if($dfrac);
	$GATKargs.="-dt $dt " if($dt);
	if($intervals)
	{	my @intervals=split(/,/,$intervals);
		foreach my $interval (@intervals)
		{	$GATKargs.="-L $interval "
		}
	}

	$GATKargs.="-intervals $dfrac " if($dfrac);
	$GATKargs.="$gatkarg " if($gatkarg);
  logentry("GATK VARIANT CALLING: ERC GVCF SINGLE-SAMPLE DISCOVERY MODE");

foreach my $file (@files)
{	my $stub=fileparse($file, qr/\.[^.]*(\.gz)?$/);
	$sys=cmd("$java $GATK -T HaplotypeCaller -R $refFasta -ERC GVCF -I $file -o $outDir/$stub.g.vcf -gt_mode DISCOVERY $GATKargs -nct $ncpu -rf BadCigar","Run HaplotypeCaller on ${file}");
}

	# FILE CLEAN UP
	logentry("FILE CLEAN UP");
	if(! $keep)
	{	logentry("Remove intermediate file directory");
		$sys=remove_tree($intermed);
	}


logentry("SCRIPT COMPLETE");
close($LOG);

exit 0;


#######################################
############### SUBS ##################




__END__