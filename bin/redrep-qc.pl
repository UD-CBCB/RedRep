#!/usr/bin/perl

my $ver="redrep-qc.pl Ver. 2.3 [08/13/2020 rev]";
my $script=join(' ',@ARGV);

use strict;

die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPLIB.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPLIB'} && -e $ENV{'REDREPLIB'} && -d $ENV{'REDREPLIB'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPUTIL.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPUTIL'} && -e $ENV{'REDREPUTIL'} && -d $ENV{'REDREPUTIL'});
die "ERROR: RedRep Installation Environment Variables not properly defined: REDREPBIN.  Please check your redrep.profile REDREP_INSTALL setting and make sure that file is sourced.\n" unless($ENV{'REDREPBIN'} && -e $ENV{'REDREPBIN'} && -d $ENV{'REDREPBIN'});

use lib $ENV{'REDREPLIB'};
use RedRep::Utils;
use RedRep::Utils qw(avgQual concat countFastq IUPAC2regexp IUPAC_RC);
use Getopt::Long qw(:config no_ignore_case);
use Carp qw(cluck confess longmess);
use Parallel::ForkManager;
use Pod::Usage;
use File::Basename qw(fileparse);
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dircopy);
use File::Path qw(make_path remove_tree);
use Filesys::Df qw(df);
use Sys::Hostname qw(hostname);

sub BC_File;
sub BCStats;
sub fastq_pair_repair;
sub filter;

### ARGUMENTS WITH NO DEFAULT
our($outDir,$intermed,$tmpdir,$debug,$no_stage_intermed,$tmp_in_outdir); 	# globals to properly handle file unstaging in case of fatal errors
my($inFile,$inFile2,$help,$manual,$force,$sepQC,$preQC,$postQC,$metaFile,$version,$stage,$int_rs_keep,$miss_hang_keep);
our($keep);


### ARGUMENTS WITH DEFAULT
my $fpAdapt		=	"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC";
my $tpAdapt		=	"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
my $logOut;										# default post-processed
my $statsOut;									# default post-processed
our $maxLen		=	9999999;
our $minLen		=	35;
my $qual			=	30;
our $maxN			=	2;
my $part			=	1;						#  Partial alignment max for BC deconv
my $mismatch	=	1;
my $ncpu			=	1;
our $verbose	= 4;

GetOptions (
	"1|i|in|in1=s"			=>	\$inFile,
	"2|in2=s"						=>	\$inFile2,
	"o|out=s"						=>	\$outDir,
	"c|meta=s"					=>	\$metaFile,
	"l|log=s"						=>	\$logOut,
	"s|stats=s"					=>	\$statsOut,

	"f|force"						=>	\$force,
	"d|debug:+"					=>	\$debug,
	"k|keep|keep_temp"	=>	\$keep,
	"V|verbose:+"				=>	\$verbose,
	"pre_qc!"						=>	\$preQC,
	"post_qc!"					=>	\$postQC,

	"per_barcode_qc"		=>	\$sepQC,

	"5|5p_adapt=s"			=>	\$fpAdapt,    			# 5' sequencing adapter
	"3|3p_adapt=s"			=>	\$tpAdapt,					# 3' sequencing adapter

	"x|maxLen=i"				=>	\$maxLen,
	"w|minLen=i"				=>	\$minLen,
	"n|max_N_run=i"			=>	\$maxN,
	"p|part=i"					=>	\$part,
	"e|mismatch=i"			=>	\$mismatch,
	"q|qual=i"					=>	\$qual,
	"keep_int_rs"				=>	\$int_rs_keep,
	"keep_miss_hang"		=>	\$miss_hang_keep,

	"t|threads|ncpu=i"	=>	\$ncpu,
	"T|tmpdir=s"				=>	\$tmpdir,
	"S|stage"						=>	\$stage,
	"tmp_in_outdir"			=>	\$tmp_in_outdir,

	"v|ver|version"			=>	\$version,
	"h|help"						=>	\$help,
	"m|man|manual"			=>	\$manual);


### VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage(-verbose => 1)  if ($help);
die "\n$ver\n\n" if ($version);
pod2usage( -msg  => "ERROR!  Required argument -i (input file 1) not found.\n", -exitval => 2) if (! $inFile);
pod2usage( -msg  => "ERROR!  Required argument -o (output directory) not found.\n", -exitval => 2)  if (! $outDir);
pod2usage( -msg  => "ERROR!  Required argument -m (metadata file) not found.\n", -exitval => 2)  if (! $metaFile);


### DEBUG MODE
if($debug) {
	require warnings; import warnings;
	$keep=1;
	$verbose=10;
	if($debug>1) {
		$verbose=100;
		require Data::Dumper; import Data::Dumper;
		require diagnostics; import diagnostics;
		$tmp_in_outdir=1 if($debug>2);
	}
}


### SET DEFAULT METHOD OF FILE PROPAGATION
my $mv = \&move;
if($keep)
{	$mv = \&copy;
}


### DECLARE OTHER GLOBALS
my $sys;												# system call variable
my $stub=fileparse($inFile, qr/\.[^.]*(\.gz)?$/);
my $stub2=fileparse($inFile2, qr/\.[^.]*(\.gz)?$/) if $inFile2;
our $manager = new Parallel::ForkManager( $ncpu );


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
mkdir($outDir) || confess("ERROR: Can't create output directory $outDir");
$logOut="$outDir/log.qc.txt" if (! $logOut);
open(our $LOG, "> $logOut");
logentry("SCRIPT STARTED ($ver)",2);


### OUTPUT FILE LOCATIONS
$tmpdir								=	get_tmpdir($tmpdir);
$intermed							=	"$tmpdir/intermed";
$intermed							=	"$outDir/intermed" if($tmp_in_outdir);
my $stage_dir					=	$tmpdir."input/";

mkdir($intermed);
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

	# cutadapt
	my $cutadapt=check_dependency("cutadapt","--version","s/\r?\n/ | /g","Available from http://code.google.com/p/cutadapt/",1);

	# fastqc
	my $fastqc=check_dependency("fastqc","--version","s/\r?\n/ | /g","Available from http://www.bioinformatics.babraham.ac.uk/projects/fastqc/",1);

	# fastx_barcode_splitter
	my $fastx_bc_split=check_dependency("fastx_barcode_splitter.pl"," |head -n 1","s/\r?\n/ | /g","Part of fastx toolkit, available from http://hannonlab.cshl.edu/fastx_toolkit/commandline.html",1);


### STAT FILE
$statsOut="$outDir/stats.txt" if (! $statsOut);
open(our $STAT, "> $statsOut");


### FILE LOCATIONS
my $BCFile_p1_base=$outDir."/barcodes_p1.txt";
my $BCFile_p2_base=$outDir."/barcodes_p2.txt";
my @BCFile_p1;
my @BCFile_p2;
my $file_concat_paired_p1="final.paired_p1.fastq";
my $file_concat_paired_p2="final.paired_p2.fastq";
my $file_concat_single_p1="final.single_p1.fastq";
my $file_concat_single_p2="final.single_p2.fastq";
my $file_concat_single="final.single.fastq";

my $dir_preQC=$outDir."/pre-fastqc";
my $dir_trim1=$outDir."/5-trim";
my $dir_deconv_p1=$outDir."/6-deconv_p1";
my $dir_deconv_p2=$outDir."/6-deconv_p2" if ($inFile2);
my $dir_recomb=$outDir."/7-recomb-mates";
my $dir_filter=$outDir."/8-filter";
my $dir_final_deconv=$outDir."/final-trimmed-deconv";
my $dir_final_concat=$outDir."/final-trimmed-concat";
my $dir_postQC=$outDir."/post-fastqc";


### STAGE Files
if($stage) {
	logentry("STAGING FILES TO $tmpdir",3);
	my @f;
	push(@f, $inFile);
	push(@f, $inFile2) if($inFile2);
	@f=stage_files($stage_dir,\@f);
	$inFile=shift(@f);
	$inFile2=shift(@f) if ($inFile2);
	logentry("Input fastq file(s) staged to $stage_dir",4);
}



### STEP 1 -- READ METADATA FILE
# Reads in metadata file and parses out informative fields into a hash called %meta

logentry("BEGIN STEP 1: PROCESS METADATA FILE",3);

my %meta;		# metadata hash of hashes with unique_id as primary key, @ targetfields as secondary keys
my %index;		# hash with p1_index_seq,p2_index_seq as key, unique_id as value
my %p1_p2_ind;	# hash of arrays mapping p1_index_seq (key) to p2_index_seq (value array);


{
	open(META, $metaFile) or pod2usage( -msg  => "ERROR!  Metadata File $metaFile not found.\n", -exitval => 2);
	my @fieldnames;
	my @targetfields=("unique_id","p1_recog_site","p1_hang_seq","p1_index_seq","p2_recog_site","p2_hang_seq","p2_index_seq");
	my @excludes=("none","None","NONE","N/A","n/a","");
	while(<META>) {
		if(@fieldnames) {
			chomp;
			my @fields=split(/\t/,$_);
			my %temp;
			for(my $i=0; $i<scalar(@fieldnames); $i++) {
				if(grep(/^$fieldnames[$i]$/,@targetfields)) {
					if(! grep(/$fields[$i]/,@excludes)) {
						$temp{$fieldnames[$i]}=$fields[$i];
					}
				}
			}

			# Assign %meta
			%{$meta{$temp{'unique_id'}}}=%temp;

			# Assign %index and # p1_p2_ind.  Needed for 2 barcode deconvolution.
			my $temp_index;
			if($temp{p2_index_seq}) {
				$temp_index=$temp{p1_index_seq}.",".$temp{p2_index_seq};
				push(@{$p1_p2_ind{$temp{p1_index_seq}}},$temp{p2_index_seq});
			} else {
				$temp_index=$temp{p1_index_seq}.",no_p2_index";
			}
			$index{$temp_index}=$temp{unique_id};
		} else {
			chomp;
			@fieldnames=split(/\t/,$_);

			#Check for required fields
			my @missing_fields;
			for(my $i=0; $i<scalar(@targetfields); $i++) {
				unless(grep(/^$targetfields[$i]$/,@fieldnames)) {
					push(@missing_fields,$targetfields[$i]);
				}
			}
			logentry_then_die("ERROR: Required field(s) ".join(', ',@missing_fields)." not contained in the metadata file header row.") if (@missing_fields);
		}
	}
	close(META);
}


### STEP 2 -- MAKE BARCODE FILES
# Uses metadata file to produce barcode file(s) required by fastx_barcode_splitter.pl

logentry("BEGIN STEP 2: CREATE BARCODE FILES",3);
my $BC_p2;
my @BCLen_p1;
my @BCLen_p2;
(@BCLen_p1)=BC_File($BCFile_p1_base,"p1",\%meta);
foreach my $l (@BCLen_p1) {
	push(@BCFile_p1,$l);
}
(@BCLen_p2)=BC_File($BCFile_p2_base,"p2",\%meta);

if($BCLen_p2[0] && @BCLen_p2!=0 && $BCLen_p2[0]!=0) {
	$BC_p2=1;
	foreach my $l (@BCLen_p2) {
		push(@BCFile_p2,$l);
	}
}

### STEP 3 -- INITIAL STATS
logentry("BEGIN STEP 3: PRODUCE INITIAL STATISTICS",3);
print $STAT "INITIAL STATISTICS\n";
$sys=countFastq("$inFile", "Count input file 1 fastq file");
print $STAT "Sequence count in input file 1: $sys";
if($inFile2) {
	$sys=countFastq("$inFile2", "Count input file 2 fastq file");
	print $STAT "Sequence count in input file 2: $sys";
}
print $STAT "\n";


### STEP 4 -- PRE-FASTQC
if($preQC) {
	logentry("BEGIN STEP 4: PRE-FASTQC",3);
	mkdir($dir_preQC);
	$sys=cmd("$fastqc --outdir $dir_preQC --format fastq --threads $ncpu --extract --quiet $inFile", "Run Pre-fastqc File 1 (P1)");
	$sys=cmd("$fastqc --outdir $dir_preQC --format fastq --threads $ncpu --extract --quiet $inFile2", "Run Pre-fastqc File 2 (P2)") if ($inFile2);
} else {
	logentry("OMITTING STEP 4: PRE-FASTQC",2);
}


### STEP 5 -- TRIMMING
	logentry("BEGIN STEP 5: TRIMMING",3);
	mkdir($dir_trim1);
	$sys=cmd("$cutadapt --quality-base 33 -q ${qual} -a $tpAdapt -m 1 -o '$dir_trim1/$stub.trim1.fastq' $inFile","Trim1-p2 (qual/3' seq adapter)");
	$sys=countFastq("$dir_trim1/$stub.trim1.fastq", "Trim1 count fastq");
	print $STAT "=================================================\nSTATISTICS AFTER TRIMMING (step 5)\n";
	print $STAT "File1 sequence count after trimming/filter step (quality/3'adapter trim): $sys";
	if($inFile2) {  #6/18/14 for some reason 3' adapter was being trimed in next step . . . changed to 5'
		$sys=cmd("$cutadapt --quality-base 33 -q ${qual} -a $fpAdapt -m 1 -o '$dir_trim1/$stub2.trim1.fastq' $inFile2","Trim1-p2 (qual/5' seq adapter)");
		$sys=countFastq("$dir_trim1/$stub2.trim1.fastq", "Trim1 count fastq");
		print $STAT "File2 sequence count after trimming/filter step (quality/3'adapter trim): $sys";
	}
	print $STAT "\n";


### STEP 6 -- BARCODE DECONVOLUTION
	logentry("BEGIN STEP 6: BARCODE DECONVOLUTION",3);

	mkdir($dir_deconv_p1);

	# PARSE P1 - SE or PE sequencing same
	foreach my $f (@BCFile_p1) {
		my $tmpFilename="$dir_trim1/$stub.trim1.fastq";
		if(-e "$dir_deconv_p1/unmatched") {
			$sys=&$mv("$dir_deconv_p1/unmatched","$dir_deconv_p1/unmatched_tmp");
			$tmpFilename="$dir_deconv_p1/unmatched_tmp";
		}
		my $temp="cat $tmpFilename | $fastx_bc_split --bol --bcfile ${BCFile_p1_base}_${f} --mismatches $mismatch --prefix '$dir_deconv_p1/' ";
		$temp.="--partial $part" if ($part);
		$sys=cmd("$temp","BC deconvolution read1, $f bp barcodes");
		if(-e "$dir_deconv_p1/unmatched_tmp") {
			#$sys=cmd("rm $dir_deconv_p1/unmatched_tmp","Remove temp file");
			$sys=unlink("$dir_deconv_p1/unmatched_tmp");
		}
	}

	# PARSE P2 - If paired end reads
	if($inFile2) {
		mkdir($dir_deconv_p2);
		foreach my $f (@BCFile_p2) {
			my $tmpFilename="$dir_trim1/$stub2.trim1.fastq";
			if(-e "$dir_deconv_p2/unmatched") {
				$sys=&$mv("$dir_deconv_p2/unmatched","$dir_deconv_p2/unmatched_tmp");
				$tmpFilename="$dir_deconv_p1/unmatched_tmp";
			}
			# If PE with 2-sided barcodes
			if($BC_p2) {
				my $temp="cat $tmpFilename | $fastx_bc_split --bol --bcfile ${BCFile_p2_base}_${f} --mismatches $mismatch --prefix '$dir_deconv_p2/' ";
				$temp.="--partial $part" if ($part);
				$sys=cmd("$temp","BC deconvolution read2");
				if(-e "$dir_deconv_p2/unmatched_tmp") {
					#$sys=cmd("rm $dir_deconv_p2/unmatched_tmp","Remove temp file");
					$sys=unlink("$dir_deconv_p2/unmatched_tmp");
				}
			} else {
			# If PE with 1-sided barcodes . . . always assume P1 is the barcode and P2 is unbarcoded . . .
			# Originally a cp, made symbolic link to save space.  NEEDS TESTING
				#$sys=cmd("ln -s $dir_trim1/$stub2.trim1.fastq $dir_deconv_p2/no_p2_index","P2 not barcoded.  Using p2 file without deconvolution.");
				$sys=&$mv("$dir_trim1/$stub2.trim1.fastq","$dir_deconv_p2/no_p2_index");
			}
		}
	}

	print $STAT "=================================================\nSTATISTICS AFTER BARCODE DECONVOLUTION (step 6)\n";
	BCStats($dir_deconv_p1);
	BCStats($dir_deconv_p2) if($inFile2);
	print $STAT "\n";


### STEP 7 -- MERGE BARCODES
logentry("BEGIN STEP 7: MERGE BARCODES",3);


{	mkdir($dir_recomb);
	opendir(DIR,"$dir_deconv_p1");
	my @files=grep { (!/^\./) } readdir DIR;
	close(DIR);

	foreach my $file_p1 (@files) {
		unless($file_p1 eq "unmatched") {
			if(! $inFile2) {  # not paired end
				#Copy file $file_p1 to $dir_recomb
				$sys=&$mv("$dir_deconv_p1/$file_p1","$dir_recomb/".$index{"$file_p1,no_p2_index"}.".singles");
			}	elsif($p1_p2_ind{$file_p1}) {  # paired end with two sided barcodes
				foreach my $file_p2 (@{$p1_p2_ind{$file_p1}}) {
					if(-e "$dir_deconv_p2/$file_p2.single_p2") {
						fastq_pair_repair("$dir_deconv_p1/$file_p1","$dir_deconv_p2/$file_p2.single_p2",$dir_recomb,$index{"$file_p1,$file_p2"},0);
						#Move temp singles file $file_p2
						$sys=&$mv("$dir_recomb/$file_p2.single_p2","$dir_deconv_p2/");
					}	else {
						fastq_pair_repair("$dir_deconv_p1/$file_p1","$dir_deconv_p2/$file_p2",$dir_recomb,$index{"$file_p1,$file_p2"},0);
						#Move temp singles file $file_p2
						$sys=&$mv("$dir_recomb/$file_p2.single_p2","$dir_deconv_p2/");
					}
				}
			} else {     #paired end with one sided barcodes
				if(-e "$dir_deconv_p2/no_p2_index.single_p2") {
					fastq_pair_repair("$dir_deconv_p1/$file_p1","$dir_deconv_p2/no_p2_index.single_p2",$dir_recomb,$index{"$file_p1,no_p2_index"},0);
					#Rename p1 singles
					$sys=move("$dir_recomb/$file_p1.single_p1","$dir_recomb/$index{$file_p1.',no_p2_index'}.single_p1");
					#Move temp singles file
					$sys=&$mv("$dir_recomb/no_p2_index.single_p2","$dir_deconv_p2/");
				}	else {
					fastq_pair_repair("$dir_deconv_p1/$file_p1","$dir_deconv_p2/no_p2_index",$dir_recomb,$index{"$file_p1,no_p2_index"},0);
					#Rename p1 singles
					$sys=move("$dir_recomb/$file_p1.single_p1","$dir_recomb/$index{$file_p1.',no_p2_index'}.single_p1");
					#Move temp p2 singles file
					$sys=&$mv("$dir_recomb/no_p2_index.single_p2","$dir_deconv_p2/");
				}
			}
		}
	}
	mkdir("$dir_recomb/singles");
	#Move broken pairs
	$sys=&$mv($_,"$dir_recomb/singles/") for glob "$dir_recomb/*.single*";
	if($inFile2) {
		mkdir("$dir_recomb/paired_p1");
		#Move singles files to recomb directory
		mkdir("$dir_recomb/paired_p2");
		#Move p1 pairs
		$sys=&$mv($_,"$dir_recomb/paired_p1/") for glob "$dir_recomb/*.paired_p1";
		#Move p2 pairs
		$sys=&$mv($_,"$dir_recomb/paired_p2/") for glob "$dir_recomb/*.paired_p2";
	}
}


### STEP 8 -- FILTER
# TRIM BARCODE AND FILTER SEQUENCES WITHOUT 5' HANG OR WITH INTERNAL RESTRICTION SITE

logentry("BEGIN STEP 8: FILTER",3);

print $STAT "=================================================\nSTATISTICS AFTER FILTERS\n";

{	mkdir("$dir_filter");

	if(! $inFile2) {		# not paired end
		mkdir("$dir_filter/singles");
		filter("p1","$dir_recomb/singles","$dir_filter/singles","single",\%meta);
	} elsif($inFile2) {		# paired end
		mkdir("$dir_filter/paired_p1");
		filter("p1","$dir_recomb/paired_p1","$dir_filter/paired_p1","paired_p1",\%meta);
		mkdir("$dir_filter/paired_p2");
		filter("p2","$dir_recomb/paired_p2","$dir_filter/paired_p2","paired_p2",\%meta);
	}
	if($inFile2 && ! $BC_p2) {	# collect paired end singles if resolvable (i.e. one barcode system)
		mkdir("$dir_filter/singles");
		filter("p1","$dir_recomb/singles","$dir_filter/singles","single_p1",\%meta);
	}
}


### STEP 9 -- FINAL RECONCILE PAIRS

logentry("BEGIN STEP 9: RECONCILE PAIRS",3);

{	mkdir($dir_final_deconv);
	mkdir($dir_final_concat);
	if($inFile2) {
		opendir(DIR,"$dir_filter/paired_p1");
		my @files=grep { (!/^\./) } readdir DIR;
		close(DIR);
		foreach my $file_p1 (@files) {
			#$manager->start and next;
			unless($file_p1 eq "discarded") {
				if($inFile2) {
					my $stub=fileparse($file_p1, qr/\.[^.]*(\.gz)?$/);
					fastq_pair_repair("$dir_filter/paired_p1/$file_p1","$dir_filter/paired_p2/$stub.paired_p2",$dir_final_deconv,$stub,0);
					cmd("cat $dir_filter/singles/$stub.single_p1 $dir_final_deconv/$stub.single_p1 > $dir_final_deconv/$stub.single_p1_new","Merge filter/final singles") if(-e "$dir_filter/singles/$stub.single_p1");
					#Delete temporary singles
					$sys=unlink("$dir_final_deconv/$stub.single_p1") if(-e "$dir_filter/singles/$stub.single_p1");
					#Rename single file $dir_final_deconv/$stub.single_p1_new
					$sys=move("$dir_final_deconv/$stub.single_p1_new","$dir_final_deconv/$stub.single_p1") if(-e "$dir_filter/singles/$stub.single_p1");
				}
			}
			#$manager->finish;
		}
	}	else {
		#Move filter/final singles
		$sys=&$mv($_,"$dir_final_deconv/") for glob "$dir_filter/singles/*single";
	}
	print $STAT "=================================================\nFINAL STATISTICS\n";
	BCStats($dir_final_deconv);
}

{	opendir(DIR,"$dir_final_deconv");
	my @files=grep { (!/^\./) } readdir DIR;
	close(DIR);
	foreach my $file (@files) {
		#Add fastq extensions
		$sys=move("$dir_final_deconv/$file","$dir_final_deconv/$file.fastq");
	}
}


### STEP 10 -- FINAL CONCATENATE

logentry("BEGIN STEP 10: PRODUCE CONCATENATED FASTQs",3);

{	opendir(DIR,"$dir_final_deconv");
	my @files=grep { (!/^\./) } readdir DIR;
	close(DIR);

	if($inFile2) {
		my @files_paired_p1 = grep(/paired_p1/,@files);
		concat("$dir_final_deconv","$dir_final_concat","$file_concat_paired_p1",\@files_paired_p1);
		my @files_paired_p2 = grep(/paired_p2/,@files);
		concat("$dir_final_deconv","$dir_final_concat","$file_concat_paired_p2",\@files_paired_p2);
		my @files_single_p1 = grep(/single_p1/,@files);
		concat("$dir_final_deconv","$dir_final_concat","$file_concat_single_p1",\@files_single_p1);
		my @files_single_p2 = grep(/single_p2/,@files);
		concat("$dir_final_deconv","$dir_final_concat","$file_concat_single_p2",\@files_single_p2);
	}
	if(! $BC_p2 || ! $inFile2) {
		my @files_single = grep(/single/,@files);
		concat("$dir_final_deconv","$dir_final_concat","$file_concat_single",\@files_single);
	}
}

### STEP 11 -- POST-fastqc

if($postQC) {
	logentry("BEGIN STEP 11: POST-FASTQC",3);

	my @files;

	mkdir($dir_postQC);
	if($sepQC) {
		opendir(DIR,"$dir_final_deconv");
		@files=grep { (!/^\./) } readdir DIR;
		close(DIR);
		foreach my $file (@files) {
			$sys=cmd("$fastqc --outdir $dir_postQC --format fastq --threads $ncpu --noextract --quiet $dir_final_deconv/$file", "Run Post-fastqc $dir_final_deconv/$file");
		}
	} else {
		opendir(DIR,"$dir_final_concat");
		@files=grep { (!/^\./) } readdir DIR;
		close(DIR);
		foreach my $file (@files) {
			$sys=cmd("$fastqc --outdir $dir_postQC --format fastq --threads $ncpu --noextract --quiet $dir_final_concat/$file", "Run Post-fastqc $dir_final_concat/$file");
		}
	}
} else {
	logentry("OMITTING STEP 11: POST-FASTQC",2);
}


### STEP 12 -- CLEAN UP

logentry("BEGIN STEP 12: CLEAN UP",3);

if(! $keep) {
	foreach my $f (@BCFile_p1) {
		#Remove Fastx Barcode Files p1
		$sys=unlink("${BCFile_p1_base}_${f}");
	}
	#Remove trim1, deconv_p1, recomb, filter directories
	$sys=remove_tree("$dir_trim1","$dir_deconv_p1","$dir_recomb","$dir_filter");
	if ($inFile2) {
		#Remove deconv_p2 directory
		$sys=remove_tree("$dir_deconv_p2");
		foreach my $f (@BCFile_p2) {
			#Remove Fastx Barcode Files p2
			$sys=unlink("${BCFile_p2_base}_${f}") if(-e $f);
		}
	}
}


### STANDARD FILE CLEAN UP
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
close($STAT);


exit 0;


#######################################
############### SUBS ##################


#######################################
### BC_File
# make fastx_barcode_splitter.pl compatible barcode file
sub BC_File {
	my $file=shift;
	my $direction=shift;
	my %meta=%{(shift)};
	my $BC=0;
	my @BCLen;
	my %seen;
	my %uniq_bc;
	logentry("Processing $direction barcode file",4);
	foreach my $sample (sort {$a cmp $b} keys %meta) {
		if ($meta{$sample}{$direction."_index_seq"}) {
			$BC=1;
			# Add index to HOA unless if has been seen previously
			unless($seen{$meta{$sample}{$direction."_index_seq"}}++) {
				push(@{$uniq_bc{length($meta{$sample}{$direction."_index_seq"})}}, $meta{$sample}{$direction."_index_seq"});
			}
		}	elsif($BC==1) {
			# All samples must have an index if any has an index
			pod2usage( -msg  => "ERROR!  If ".$direction."_index_seq is set for one sample, it must be set for all samples.  No ".$direction."_index_seq for sample $sample.\n", -exitval => 2);
		}
	}
	foreach my $len (sort { $b <=> $a } keys %uniq_bc) {
		push (@BCLen, $len);
		open(BC,"> ".$file."_".$len);
		foreach my $bc (@{$uniq_bc{$len}}) {
			print BC $bc."\t".$bc."\n";
		}
		close(BC);
	}
	return (@BCLen);
}


#######################################
### BCStats
# Count and record per barcode sequence counts from a directory of deconvoluted fastq files
sub BCStats {
	my $inDir=shift;
	my $sys;
	my $total;
	logentry("Processing Sequence Count Stats for $inDir",4);

	opendir(DIR,"$inDir");
	my @files=grep { (!/^\./) } readdir DIR;
	close(DIR);

	print $STAT "Sequence count ($inDir):\n"; #cleanup ba
	foreach my $file (sort @files) {
		if(-s $inDir."/".$file) {
			$sys=countFastq("$inDir/$file","File count $inDir/$file");
		}	else {
			$sys="0\n";
		}
		print $STAT "$file\t$sys";
		$total+=$sys;
	}
	print $STAT "TOTAL\t$total\n";
	print $STAT "\n";
}


#######################################
### fastq_pair_repair
# matches pairs and segregated broken pair sequences from a pair of input fastq files
sub fastq_pair_repair {
	my $PE1=shift;
	my $PE2=shift;
	my $outDir=shift;
	my $outName=shift;
	my $no_singles=shift;

	logentry("Merging paired end mates: $PE1 $PE2",5);

	# PARSE OUTPUT FILEBASE
	my $out1=fileparse($PE1, qr/\.[^.]*(\.gz)?$/);
	my $out2=fileparse($PE2, qr/\.[^.]*(\.gz)?$/);

	# FILE HANDLES
	my($DATA,$PAIR1,$PAIR2,$SNGL1,$SNGL2);

	# OPEN PE1
	open ($DATA,$PE1) or logentry_then_die("Cannot open $PE1. $!");

	my %seqs;

	# REG EXP FOR FASTQ HEADERS
	my $hdr_ptrn;
	# Illumina 1-1.7:	@HWUSI-EAS100R:6:73:941:1973#0/1
	# Illumina 1.8+:	@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	$hdr_ptrn='^\@(\S+)[\/ ][12]';

	# PROCESS PE1
	while(<$DATA>) {
		if(/$hdr_ptrn/) {
			$seqs{$1}=$_;
			$seqs{$1}.=<$DATA>.<$DATA>.<$DATA>;
		} else {
			logentry_then_die("File format error in $PE1 near line ".$..". value=$_");
		}
	}
	close $DATA;

	open($DATA,$PE2) or logentry_then_die("Cannot open $PE2.$!");
	open($PAIR1,"> $outDir/$outName.paired_p1") or logentry_then_die("Cannot open output file $outDir/$outName.paired_p1. $!");
	open($PAIR2,"> $outDir/$outName.paired_p2") or logentry_then_die("Cannot open output file $outDir/$outName.paired_p2. $!");
	if(! $no_singles) {
		open($SNGL1,"> $outDir/$out1.single_p1") or logentry_then_die("Cannot open output file $outDir/$out1.single_p1. $!");
		open($SNGL2,"> $outDir/$out2.single_p2") or logentry_then_die("Cannot open output file $outDir/$out1.single_p2. $!");
	}


	# PROCESS PE2 AND OUTPUT PAIRS/SINGLES FROM PE2
	while(<$DATA>) {
		if (/$hdr_ptrn/) {
			if ($seqs{$1}) {
				print $PAIR1 $seqs{$1};
				undef $seqs{$1};
				$_.=<$DATA>.<$DATA>.<$DATA>;
				print $PAIR2 $_;
			}	else {
				$_.=<$DATA>.<$DATA>.<$DATA>;
				print $SNGL2 $_ if (! $no_singles);
			}
		}	else {
			logentry_then_die("File format error in $PE2 near line ".$..".value=$_");
		}
	}

	# PRINT SINGLES FROM PE1
	if(! $no_singles) {
		foreach my $key(keys %seqs)	{
			if (($seqs{$key})&&($seqs{$key} ne ""))	{
				print $SNGL1 $seqs{$key};
			}
		}
		close $SNGL1;
		close $SNGL2;
		my $sys
	}

	close $DATA;
	close $PAIR1;
	close $PAIR2;
}


#######################################
### filter
# perform step 4 filtering on sequences
sub filter {
	my $direction=shift;  #p1 or p2
	my $in_dir=shift;
	my $out_dir=shift;
	my $suffix=shift;
	my %meta=%{(shift)};

	logentry("Filtering $suffix",4);

	opendir(DIR,"$in_dir");
	my @files=grep { (!/^\./) } readdir DIR;
	close(DIR);
	mkdir("$out_dir/discarded");
	mkdir("$out_dir/tmp");

	print $STAT "Per sample sequence counts ($suffix):\n";
	print $STAT "sample\tafter_deconv\tinternal_RS\tmissing_hang_seq\tN_run\tmax_length\tmin_length\ttotal_removed\tfinal_total\n";

	foreach my $file (sort @files) {
		$manager->start and next;
		open(TMP,"> $out_dir/tmp/$file.stats.tmp");
		unless($file eq "unmatched" || $file eq "discarded") {
			my $stub=fileparse($file, qr/\.[^.]*(\.gz)?$/);
			my $BCLen=length($meta{$stub}{$direction."_index_seq"});
			my $BCLen_2=$BCLen-1 if($BCLen);

			open(DAT, "$in_dir/$file");
			open(OUT, ">$out_dir/$stub.$suffix");
			open(DISC, ">$out_dir/discarded/$stub.discard_$suffix");

			my $count=0;
			my $keep_seq=0;
			my $trimLen=0;
			my $seq="";
			my $hang=$meta{$stub}{$direction.'_hang_seq'};
			# Convert hang seq to regexp pattern (allows degenerate base matching)
			$hang = IUPAC2regexp($hang);
			my $stat_cnt=0;
			my $stat_intRS=0;
			my $stat_fphang_miss=0;
			my $stat_Nrun=0;
			my $stat_maxLen=0;
			my $stat_minLen=0;
			my $stat_tot=0;
			my $stat_disc=0;
			my $max_N_run=$main::maxN+1;

			while(<DAT>) {
				my $datIn=$_;
				$count++;
				if($count==2) {
					chomp($seq);
					# missing hang seq
					if($BCLen && $datIn=~ /^.{$BCLen_2,$BCLen}$hang/) {
						$datIn=~ s/^(.{$BCLen_2,$BCLen})($hang)/$2/;
						$trimLen=length($1);
						$seq .= " EXPECT_BC=".$meta{$stub}{$direction."_index_seq"}." FOUND_BC=$1 SAMPLE=$stub";
						$keep_seq=1;
					} elsif(! $BCLen && ($datIn =~ /^$hang/)) {
						$seq .= " EXPECT_BC=NONE FOUND_BC=NONE SAMPLE=$stub";
						$keep_seq=1;
					}	else {
						$stat_fphang_miss++;
						if($miss_hang_keep) {
							$keep_seq=1;
							my($exp_bc, $found_bc);
							if($BCLen) {
								$datIn=~ s/^(.{$BCLen})// ;
								$trimLen=length($1);
								$found_bc=$1;
								$exp_bc=$meta{$stub}{$direction."_index_seq"};
							}	else {
								$exp_bc="NONE";
								$found_bc="NONE";
							}
							$seq .= " EXPECT_BC=$exp_bc FOUND_BC=$found_bc SAMPLE=$stub MISS_HANG=YES";
						}
					}

					# internal restriction site
					#FIXED in 3.0 added P2 as separate conditional push since putting value of "NA" was causing that pattern to get filtered
					my @RS = ($meta{$stub}{p1_recog_site});
					push(@RS,$meta{$stub}{p2_recog_site}) if($meta{$stub}{p2_recog_site});
					foreach my $r (@RS) {
						my $r1=IUPAC2regexp($r);
						my $r2=IUPAC2regexp(IUPAC_RC($r));
						if($datIn =~ /$r1/ || $datIn =~ /$r2/) {
							$stat_intRS++;
							if($int_rs_keep) {
								$seq .= " INT_RS=YES";
							}	else {
								$keep_seq=0 unless($int_rs_keep);
							}
						}
					}

					# internal runs of N's
					if($datIn =~ /N{$max_N_run}/i) {
						$stat_Nrun++;
						$keep_seq=0;
					}

					# too long
					if(length($datIn) > $main::maxLen) {
						$stat_maxLen++;
						$keep_seq=0;
					}

					# too short
					if(length($datIn)<$main::minLen) {
						$stat_minLen++;
						$keep_seq=0;
					}
					$seq .= "\n";
				}
				if($count==4 && $trimLen) {
					$datIn =~ s/^.{$trimLen}//;
				}
				$seq.=$datIn;
				if($count==4) {
					$stat_cnt++;
					my $qualIn=$datIn;
					chomp($qualIn);
					my ($qavg,$len);

					if($qualIn) {
						($qavg,$len)=avgQual($qualIn);
					}	else {
						($qavg,$len)=('-',0);
					}
					$seq =~ s/(SAMPLE=\S+)/$1 LENGTH=$len MEAN_QUAL=$qavg/;

					if($keep_seq==1) {
						print OUT $seq;
						$stat_tot++;
					} else {
						print DISC $seq;
						$stat_disc++;
					}
					$keep_seq=0;
					$count=0;
					$trimLen=0;
					$seq="";
				}
			}
			print TMP "$stub\t$stat_cnt\t$stat_intRS\t$stat_fphang_miss\t$stat_Nrun\t$stat_maxLen\t$stat_minLen\t$stat_disc\t$stat_tot\n";

			close(OUT);
			close(DAT);
			close(DISC);
		}
		close(TMP);
		$manager->finish;
	}
	$manager->wait_all_children;
	sleep 1;
	$sys=cmd("cat $out_dir/tmp/*.stats.tmp > $out_dir/tmp/all.stats.tmp","Concatenate temporary stats files");
	my $counts=cmd("awk 'BEGIN {FS=OFS=".'"\t"'."} NR == 1 { n2 =\$2; n3 = \$3; n4 = \$4; n5 = \$5; n6 = \$6; n7 = \$7; n8 = \$8; n9 = \$9; next } { n2 += \$2; n3 += \$3; n4 += \$4; n5 += \$5; n6 +=\$6; n7 += \$7; n8 += \$8; n9 += \$9 } END { print n2, n3, n4, n5, n6, n7, n8, n9 }' $out_dir/tmp/all.stats.tmp","Format Stats File");
	$sys=cmd("sort -k1,1 $out_dir/tmp/all.stats.tmp","Sort temp stats file");
	print $STAT $sys;
	print $STAT "TOTAL\t$counts\n";
	print $STAT "\n";
	$sys=remove_tree("$out_dir/tmp") unless ($main::keep);
}


__END__


#######################################
########### DOCUMENTATION #############
=pod

=head1 NAME

redrep-qc.pl -- Illumina QC, trimming, filtering for Reduced Representation Analysis

=head1 SYNOPSIS

 redrep-qc.pl --in FILENAME [--in2 FILENAME]--out DIRNAME --meta FILENAME [PARAMETERS]
                     [--help] [--manual]
=head1 DESCRIPTION

Performs quality evaluation, trimming, and filtering of Illumina sequenced reduced representation libraries.

=head1 OPTIONS

=head2 REQUIRED PARAMETERS

=over 3

=item B<-1, -i, --in, --in1>=FILENAME

Input file (single file or first read) in fastq format. (Required)

=item B<-2, --in2>=FILENAME

Input file (second read) in fastq format. (Required, if paired end).  NOTE: Paired end mode is experimental, use at own risk.

=item B<-o, --out>=DIRECTORY_NAME

Output directory. (Required)

=item B<-c, --meta>=FILENAME

Metadata file in tab delimited format.  Must contain header row with at least the following column headings:  unique_id,p1_recog_site,p1_hang_seq,p1_index_seq,p2_recog_site,p2_hang_seq,[p2_index_seq].  Empty fields in a required column may be left blank or filled with: N/A. (Required)

=back

=head2 OPTIONAL PARAMETERS

=head3 Program Specific Parameters

=head4 Adapter Parameters

=over 3

=item B<-5, --5p_adapt>

Specify 5' sequencing adapter.  Default 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC'.

=item B<-3, --3p_adapt>

Specify 3' sequencing adapter.  Default 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'.

=back

=head4 Barcode Parameters

=over 3

=item B<-p, --part>=integer

Allow barcode to be -p bp shorter than the specified barcode.  -e parameter must be greater than or equal to -p.  (default=1)

=item B<-e, --mismatch>=integer

Allow barcode to have -e mismatches from the specified barcode.  Gaps and missing ends (-p) count as mismatches.  (default=1)

=back

=head4 QC Reporting Parameters

=over 3

=item B<--pre_qc>

Perform pre-trimming quality report.

=item B<--post_qc>

Perform post-trimming quality report.

=item B<--per_barcode_qc>

Run post-trimming quality reports on each barcode instead of on the full dataset.

=back

=head4 Sequence Trimming and Filtering Parameters

=over 3

=item B<-x, --maxLen>=integer

Maximum sequence length cutoff.  default=9999999

=item B<-w, --minLen>=integer

Minimum sequence length cutoff.  default=35

=item B<-n, --max_N_run>=integer

Maximum number of consecutive N's to allow in middle of trimmed sequence.  default=2

=item B<-q, --qual>=integer(0-93)

Quality cut-off for end-trimming.  Performed using the BWA algorithm.  default=30

=item B<--keep_miss_hang>

Keep sequences that lack the expected overhang sequence.

=item B<--keep_int_rs>

Keep sequences with internal restriction site.

=back

=head3 Program Behavior and Resource Control

=over 3

=item B<-f, --force>

If output directory exists, force overwrite of previous directory contents.

=item B<-k,--keep>

Retain temporary intermediate files.

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

Log file output path. [ Default output-dir/log.qc.txt ]

=item B<-s, --stats>=FILENAME

Stats file output path. [ Default output-dir/stats.txt ]

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

=item 0.1 - 10/10/2012: Stable base functionality for single end reads

=item 1.0 - 10/29/2012: Major overhaul adding paired end read handling, improved statistics, file handling, logging

=item 1.1 - 10/31/2012: Added capability for gzipped input files, fixed fastq header bug, added TOTAL row to barcode statistics

=item 1.2 - 10/31/2012: Added basic multi-threading for filtering

=item 1.3 - 11/1/2012: Added sort to stats to compensate for multithreading

=item 1.4 - 11/2/2012: Bugfix, -1 and -2 full file path handled

=item 1.5 - 11/5/2012: Bugfix, reworked wait on children for threads

=item 1.6 - 1/21/2013: Bugfix, -m flag had two purposes

=item 1.7 - 3/28/2014: add version output to log; check ENV for $REDREPBIN

=item 1.8 - 6/18/2014: add support for variable length barcodes, bug with trimming of 5p adapter on second PE read fixed

=item 1.81 - 10/7/2014: fixed bug in 2.8 that was masking hangseq with "..."

=item 1.9 - 3/23/2015: variable length barcodes, included reworking of barcode removal and hang seq setection, moved from system calls to perl libraries for file handling (with a few exceptions)

=item 2.0 - 11/21/2016: additional variable length barcode functionality.  Major version release.

=item 2.1 - 2/1/2017: Code cleanup.

=item 2.11 - 10/9/2019: Minor code cleanup.  Verbosity settings added.  Last version with GATK 3.x compatibility.

=item 2.2 - 10/28/2019: Reversed QC report parameters, now default is no pre or post QC report.  Many enhancements including parallelization, documentation, and logging.

=item 2.3 - 8/1/2020: Fixed input file staging and disk space reporting bugs.

=back

=head1 DEPENDENCIES

=head2 Requires the following external programs be in the system PATH:

=over 3

=item fastqc v0.11.5

=item cutadapt v1.14

=item fastx_barcode_splitter.pl v0.0.13

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
