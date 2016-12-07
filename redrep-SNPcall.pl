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

=back

=head1 DEPENDENCIES

=item samtools (tested with version 1.3.1)

=item GATK (tested with version 3.5-0)

=item java (tested with version 1.8.0_111-b16 OpenJDK Runtime Environment)

=item picard-tools (tested with version 2.4.1)

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
my($in,$outDir,$help,$manual,$force,$metaFile,$keep,$debug,$version,$diff,$refFasta,$dcov,$dfrac,$dt);


### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
my $ncpu		=	1;
my $mem=50;


GetOptions (	"i|in=s"							=>	\$in,
				"o|out=s"							=>	\$outDir,
				"r|ref=s"							=>	\$refFasta,
				"l|log=s"							=>	\$logOut,

				"f|force"							=>	\$force,
				"d|debug"							=>	\$debug,
				"k|keep_temp"						=>	\$keep,

				"t|threads|ncpu=i"					=>	\$ncpu,
				"mem=i"								=>	\$mem,

				"dcov|downsample_to_coverage=i"	=> \$dcov,
				"dfrac|downsample_to_fraction=f"	=> \$dfrac,
				"dt|downsampling_type=s"			=> \$dt,

				"v|ver|version"						=>	\$version,
				"h|help"							=>	\$help,
				"m|man|manual"						=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
my $ver="redrep-SNPcall.pl Ver. 2.0 (11/18/2016 rev)";
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Argument -i (input file/directory) missing.\n", -exitval => 2) if (! $in);
pod2usage( -msg  => "ERROR!  Required argument -r (reference fasta) missing.\n", -exitval => 2) if (! $refFasta);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) missing.\n", -exitval => 2)  if (! $outDir);
pod2usage( -msg  => "ERROR!  Arguments -dcov and -dfrac are incompatible, chose one.\n", -exitval => 2)  if ($dcov && $dfrac);


if($debug)
{	require warnings; import warnings;
	require Data::Dumper; import Data::Dumper;
	$keep=1;
}


### DECLARE OTHER GLOBALS
my $sys;												# system call variable

#our $manager = new Parallel::ForkManager( $ncpu );


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


### BIN/SCRIPT LOCATIONS
my %ver;
my %path;
my $execDir=$ENV{'REDREPBIN'};
my $java = "java -Xmx${mem}g -d64 -jar ";
$ver{'java'}=`java -version 2>&1|grep version`;
$path{'java'}=`which java`;
my $picard = $ENV{'PICARDJARS'}."/picard.jar";
$ver{'picard'}=`$java $picard AddCommentsToBam --version 2>&1`;
my $GATK = $ENV{'GATKJAR'};
$ver{'gatk'}=`$GATK --version 2>&1`;
$path{'gatk'}=$GATK;
my $samtools = "samtools";
$ver{'samtools'}=`$samtools 2>&1|grep Version`;
$path{'samtools'}=`which $samtools`;
chomp %ver;
chomp %path;


### CREATE LOG FILES
$logOut="$outDir/log.SNP.txt" if (! $logOut);

open(our $LOG, "> $logOut");
print $LOG "$0 $script\n";
print $LOG "RedRep Scripts: $execDir\n";
print $LOG "Dependency: $path{'java'} ($ver{'java'})\n";
print $LOG "Dependency: $path{'samtools'} ($ver{'samtools'})\n";
print $LOG "Dependency: $path{'gatk'} ($ver{'gatk'})\n";
logentry("SCRIPT STARTED ($ver)");


### OUTPUT FILE LOCATIONS
my $gatk_out="$outDir/combined.SNP.vcf";

### CHECK FOR EXTERNAL DEPENDENCIES
# Not yet implemented



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
#	$bam_files="-I $bam_files";
	my $GATKargs.="-dcov $dcov " if($dcov);
	$GATKargs.="-dfrac $dfrac " if($dfrac);
	$GATKargs.="-dt $dt " if($dt);
    logentry("GATK VARIANT CALLING: MULTI-SAMPLE DISCOVERY MODE");
	system("$java $GATK -T UnifiedGenotyper -R $refFasta -I $bam_files -o $gatk_out -gt_mode DISCOVERY $GATKargs -nt $ncpu -rf BadCigar");


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
### logentry
# Enter time stamped log entry
sub logentry
{	my $message=shift;
	print $LOG POSIX::strftime("%m/%d/%Y %H:%M:%S > $message\n", localtime);
}



__END__
