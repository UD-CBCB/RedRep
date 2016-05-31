#!/usr/bin/perl

# MANUAL FOR redrep-cluster.pl

=pod

=head1 NAME

redrep-cluster.pl -- De novo clustering of reduced representation libraries

=head1 SYNOPSIS

 redrep-cluster.pl --in FILENAME [--in2 FILENAME]--out DIRNAME --meta FILENAME [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Accepts fastq output from redrep-qc.pl and performs clustering of reduced representaiton tags without a priori reference.
 
=head1 OPTIONS

=over 3

=item B<-1, -i, --in, --in1>=FILENAME

Input file (single file or first read) in fastq format. (Required) 

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required) 

=item B<-l, --log>=FILENAME

Log file output path. [ Default output-dir/log.txt ]

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k, --keep_temp>

Retain temporary intermediate files.

=item B<-p, --id>

Minimum identity cutoff for clustering (Default=0.98)

=item B<-n, --minlen>

Minimum length for a cluster seed sequence (Default=35)

=item B<-q, --minqual>

Minimum quality for a cluster seed sequence (Default=20)

=item B<--minsize>

Minimum cluster size

=item B<--hist_stats>

Run the histogram statistics.  May dramatically SLOW execution time.

=item B<--size_bin_start>

Starting size of histogram bins for cluster size statistics (default=1)

=item B<--size_bin_num>

Number of histogram bins for cluster size statistics (default=50)

=item B<--size_bin_width>

Width of histogram bins for cluster size statistics (default=250)


=item B<--len_bin_start>

Starting size of histogram bins for cluster seed length statistics (default=81)

=item B<--len_bin_num>

Number of histogram bins for cluster seed length statistics (default=4)

=item B<--len_bin_width>

Width of histogram bins for cluster seed length statistics (default=5)



=item B<-t, --threads, ==ncpu>=integer

Number of cpu's to use for threadable operations. (NOT CURRENTLY IMPLEMENTED)

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

=item 0.0 - 11/5/2012: Draft1

=item 1.0 - 12/5/2012: Release

=item 1.1 - 1/31/2013: Default is now to NOT run hist stats

=item 1.2 - 3/28/2014: Added dependency version to log; import $REDREPBIN from ENV

=item 1.3 - 11/3/2015: Add minsize (minimum cluster size) option


=back

=head1 DEPENDENCIES

=item usearch (tested with version usearch_i86linux32 v6.0.307)

=head2 Requires the following Perl libraries:

=over 3

=item strict

=item Getopt::Long

=item File::Basename

=item Pod::Usage

=item POSIX

=item

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

Copyright 2012-2014 Shawn Polson, Randall Wisser, Keith Hopper.  
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
#use Parallel::ForkManager;

sub binClass;
sub cmd;
sub logentry;


### ARGUMENTS WITH NO DEFAULT
my($inFile,$inFile2,$outDir,$help,$manual,$force,$metaFile,$keep,$debug,$version,$hist_stats);


### ARGUMENTS WITH DEFAULT
my $logOut;									# default post-processed
#my $statsOut;								# default post-processed
my $ncpu		=	1;
my $id			=	0.98;
my $minlen		=	35;
my $minqual		=	20;
my $minsize		=	1;
my $sizeBinStart=	1;
my $sizeBinSize = 	250;
my $sizeBinNum	=	50;
my $sizeBinEnd	=	999999999;  # hard-coded max
my $lenBinStart=	81;
my $lenBinSize = 	5;
my $lenBinNum	=	4;
my $lenBinEnd	=	999;		# hard-coded max



GetOptions (	
				"1|i|in|in1=s"				=>	\$inFile,
				"2|in2=s"					=>	\$inFile2,
				"o|out=s"					=>	\$outDir,
#				"c|meta=s"					=>	\$metaFile,
				"l|log=s"					=>	\$logOut,
#				"s|stats=s"					=>	\$statsOut,

				"f|force"					=>	\$force,
				"d|debug"					=>	\$debug,
				"k|keep_temp"				=>	\$keep,

				"p|id=s"					=>	\$id,
				"n|minlen=i"				=>	\$minlen,
				"q|minqual=f"				=>	\$minqual,
				"minsize=i"				=>	\$minsize,
				
				"hist_stats"				=>	\$hist_stats,
				"size_bin_start=i"			=>	\$sizeBinStart,
				"size_bin_num=i"			=>	\$sizeBinNum,
				"size_bin_width=i"			=>	\$sizeBinSize,
				"len_bin_start=i"			=>	\$sizeBinStart,
				"len_bin_num=i"				=>	\$lenBinNum,
				"len_bin_width=i"			=>	\$lenBinSize,
								
				"t|threads|ncpu=i"			=>	\$ncpu,
				
				"v|ver|version"				=>	\$version,
				"h|help"					=>	\$help,
				"m|man|manual"				=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
my $ver="redrep-cluster.pl Ver. 1.3 (11/3/2015 rev)";
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Required argument -i (input file 1) not found.\n", -exitval => 2) if (! $inFile);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) not found.\n", -exitval => 2)  if (! $outDir);
#pod2usage( -msg  => "ERROR!  Required argument -m (metadata file) not found.\n", -exitval => 2)  if (! $metaFile);

if($debug)
{	require warnings; import warnings;
	require Data::Dumper; import Data::Dumper;
	$keep=1;
}

# Make sure hist bin size parameters in spec
$sizeBinEnd=($sizeBinSize*$sizeBinNum)+$sizeBinStart+99 if(($sizeBinSize*$sizeBinNum)+$sizeBinStart > $sizeBinEnd);
$lenBinEnd=($lenBinSize*$lenBinNum)+$lenBinStart+99 if(($lenBinSize*$lenBinNum)+$lenBinStart > $lenBinEnd);


### DECLARE OTHER GLOBALS
my $sys;												# system call variable
my $stub=fileparse($inFile, qr/\.[^.]*(\.gz)?$/);
#my $stub2=fileparse($inFile2, qr/\.[^.]*(\.gz)?$/);

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


### CREATE OUTPUT DIR
mkdir($outDir);
mkdir($outDir."/intermed");


### BIN/SCRIPT LOCATIONS
my %ver;
my %path;
my $execDir=$ENV{'REDREPBIN'};
my $sort_fastq_sh="$execDir/utilities/sort-fastq.sh";
my $fastq2fasta_pl="$execDir/utilities/fastq2fasta.pl";
my $fasta_abbrev_pl="$execDir/utilities/fasta-abbrev.pl";
my $usearch5="usearch5";
my $usearch8="usearch8";
$ver{'usearch'}=`$usearch8 --version 2>&1`;
$path{'usearch'}=$usearch8;
my $clstr_sort_by_pl="clstr_sort_by.pl";
my $clstr_sort_prot_by_pl="clstr_sort_prot_by.pl";
my $plot_len1_pl="plot_len1.pl";
my $CDHitClustComp_pl="$execDir/utilities/CDHitClustComp.pl";
chomp %ver;
chomp %path;




### CREATE LOG FILES
$logOut="$outDir/log.cluster.txt" if (! $logOut);

open(our $LOG, "> $logOut");
print $LOG "$0 $script\n";
print $LOG "RedRep Scripts: $execDir\n";
print $LOG "Dependency: $path{'usearch'} ($ver{'usearch'})\n";
logentry("SCRIPT STARTED ($ver)");



### CHECK FOR EXTERNAL DEPENDENCIES
# Not yet implemented


### OUTPUT FILE LOCATIONS
my $intermed="$outDir/intermed";
my $fq_sorted="$intermed/$stub.sort.fastq";
my $fasta_sorted="$intermed/$stub.sort.fasta";
my $fasta_header="$intermed/$stub.sort.header.txt";
my $fasta_sorted_abbrev="$intermed/$stub.sort.abbrev.fasta";
my $splitDir="$intermed/split";
my $uc_centroid="$intermed/$stub.uc.centroid.fasta";
my $uc_centroid_filter="$intermed/$stub.uc.centroid.filter.fasta";
my $uc_consensus="$intermed/$stub.uc.consensus.fasta";
my $uc_seeds="$intermed/$stub.uc.seeds.fasta";
my $uc_uc="$outDir/$stub.uc";
my $uc_uc_seeds="$intermed/$stub.seeds.uc";
my $uc_clstr="$intermed/$stub.uc.clstr";
my $uc_hist_stats="$outDir/$stub.hist_stats.tsv";
my $uc_cluster_freq="$outDir/$stub.cluster_freq.tsv";


############
### MAIN


	logentry("SORT FASTQ FILE");
	$sys=cmd("$sort_fastq_sh $inFile $fq_sorted $minlen $minqual","Sort fastq file");
	logentry("CONVERT FASTQ TO FASTA FORMAT");
	$sys=cmd("$fastq2fasta_pl $fq_sorted > $fasta_sorted","Convert fastq to fasta");
	logentry("ABBREVIATE FASTA HEADERS");
	$sys=cmd("$fasta_abbrev_pl $fasta_sorted $fasta_sorted_abbrev $fasta_header","Abbreviate fasta headers");

	# Split Fasta File if larger than 4GB
	my $filesize=-s $fasta_sorted_abbrev;
	my $slice=$fasta_sorted_abbrev;
	my $upperlim=4000000000;
	my $lowerlim=3950000000;
	my $sizeinout="-sizeout";

	if($filesize > $upperlim)
	{	open(DAT,"$fasta_sorted_abbrev");
		my $it=1;
		my $slicesize;
		$sys=cmd("mkdir $splitDir$it","Make split $it directory");
		$slice="$splitDir$it/$it.fasta";
		logentry("PRODUCING FASTA SPLIT $it");
		open(SUB, ">$slice");
		while(<DAT>)
		{	print SUB $_;
			$slicesize+=length($_);
			my $tmp=<DAT>;
			print SUB $tmp;
			$slicesize+=length($tmp);
			if($slicesize>$lowerlim)
			{	close(SUB);
				logentry("RUNNING UCLUST SPLIT $it");
				$sys=cmd("$usearch8 --cluster_smallmem $slice --leftjust $sizeinout --id $id --centroids $splitDir$it/$it.centroid --uc $splitDir$it/$it.uc --log $splitDir$it/uclust_$it.log","UCLUST FASTA SPLIT $slice");
				$slice="$splitDir$it/$it.centroid";
				$it++;
				$sys=cmd("mkdir $splitDir$it","Make split $it directory");
				$sys=cmd("cp $slice $splitDir$it/$it.fasta","Roll over slice input");
				$slice="$splitDir$it/$it.fasta";
				logentry("PRODUCING FASTA SPLIT $it");
				$sizeinout="-sizein -sizeout";
				open(SUB, ">>$slice");
				$slicesize=-s $slice;
				pod2usage( -msg  => "Size of input file exceeds capability of 32-bit usearch.\n", -exitval => 2) if ($slicesize>$lowerlim);
			}
		}
		close(SUB);
	}

	# UCLUST IF SMALLER THAN 4GB OR FINAL SPLIT
	logentry("RUNNING UCLUST FINAL");
	$sys=cmd("$usearch8 --cluster_smallmem $slice --leftjust $sizeinout --id $id --centroids $uc_centroid --uc $uc_uc --log $outDir/uclust_final.log","UCLUST FINAL RUN: $slice");
	$sys=cmd("$usearch8 --sortbysize $uc_centroid --fastaout $uc_centroid_filter --minsize $minsize","UCLUST CENTROID FILTER FINAL RUN: $slice") if($minsize>1);

	$sys=cmd("grep '^S' $uc_uc > $uc_uc_seeds","Make seeds uc file");
	$sys=cmd("$usearch5 --uc2fasta $uc_uc_seeds --input $fasta_sorted_abbrev --output $uc_seeds","Make seeds fasta file");
	

	# CONVERT UC TO CLSTR
	logentry("CONVERT UC FILE TO CLSTR FORMAT");
	$sys=cmd("$usearch5 --uc2clstr $uc_uc --output $uc_clstr","Convert uc to clstr format");

	# SORT CLUSTERS
	logentry("SORT CLUSTERS BY SIZE");
	$sys=cmd("$clstr_sort_by_pl $uc_clstr no > $uc_clstr.tmp","Sort clusters by size");
	logentry("SORT CLUSTERS BY SEED LENGTH");
	$sys=cmd("$clstr_sort_prot_by_pl len $uc_clstr.tmp > $uc_clstr.sort","Sort clusters by sequence length");
	
	
	
	if($hist_stats)
	{
		
		# DEFINE CLUSTER SIZE BINS FOR STATS
		my $sizeBinStr;
		for (my $j=0;$j<$sizeBinNum; $j++)
		{	$sizeBinStr .= "," if($sizeBinStr);
			$sizeBinStr .= binClass($j,$sizeBinSize,$sizeBinStart);
		}
		{	$sizeBinStr .= "," if($sizeBinStr);
			my $bottom=($sizeBinNum*$sizeBinSize)+$sizeBinStart;
			my $top=$sizeBinEnd;
			$sizeBinStr .= $bottom."-".$top;
		}
				
		# DEFINE CLUSTER LENGTH BINS FOR STATS
		my $lenBinStr;
		for (my $j=0;$j<$lenBinNum; $j++)
		{	$lenBinStr .= "," if($lenBinStr);
			$lenBinStr .= binClass($j,$lenBinSize,$lenBinStart);
		}
		{	$lenBinStr .= "," if($lenBinStr);
			my $bottom=($lenBinNum*$lenBinSize)+$lenBinStart;
			my $top=$lenBinEnd;
			$lenBinStr .= $bottom."-".$top;
		}
		
		# RUN HISTOGRAM STATS
		logentry("RUN HISTOGRAM (CLUSTER SIZE/LENGTH) STATS");
		$sys=cmd("$plot_len1_pl $uc_clstr.sort \\$sizeBinStr \\$lenBinStr > $uc_hist_stats","Run Cluster Stats (Histogram)");
	}
	else
	{	logentry("SKIPPING HISTOGRAM (CLUSTER SIZE/LENGTH) STATS: --no_hist_stats flag used");
	}
		
	# PARSE CLUSTER COUNTS
	logentry("PARSE CLUSTER COUNTS");
	$sys=cmd("$CDHitClustComp_pl $uc_clstr $uc_cluster_freq","Parse Cluster Counts");
	
	
	# FILE CLEAN UP
	logentry("FILE CLEAN UP");
	$sys=cmd("mv $uc_clstr.sort $outDir","Move clstr sort file");
	$sys=cmd("mv $uc_centroid $outDir","Move centroid fasta file");
	$sys=cmd("mv $uc_centroid_filter $outDir","Move filtered centroid fasta file") if(-e $uc_centroid_filter);
	if(! $keep | !$debug)
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
	{	return $sys;
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

