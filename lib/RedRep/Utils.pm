package RedRep::Utils;

my $ver="RedRep::Utils Ver. 2.3 [08/13/2020 rev]";

use strict;
use lib $ENV{'REDREPLIB'};
use Exporter qw(import);
use Pod::Usage;
use POSIX qw(strftime);
use Carp qw(cluck confess longmess);
use Cwd qw(abs_path cwd);
use File::Basename qw(fileparse dirname);
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dircopy);
use File::Path qw(make_path remove_tree);
use File::stat;
use File::Which qw(which);
use String::Random;
use Time::Seconds;
#use IO::Zlib;

our @EXPORT    = qw(check_dependency check_ref_fasta cmd find_job_info get_dir_filelist get_path_bwa get_path_gatk get_path_samtools get_tmpdir logentry logentry_then_die script_time stage_files stage_gdb);
our @EXPORT_OK = qw(archive_gdb avgQual build_argument_list build_contig_list build_fasta_dictionary build_fasta_index build_fofn build_listfile build_vcf_index check_jar cmd_STDOUT concat countFastq check_genomicsdb gdb_history get_timestamp IUPAC2regexp IUPAC_RC read_fodn read_fofn read_listfile retrieve_contigs round split_in_files tarball_dir);

sub archive_gdb;
sub avgQual;
sub build_argument_list;
sub build_contig_list;
sub build_fasta_dictionary;
sub build_fasta_index;
sub build_fofn;
sub build_listfile;
sub build_vcf_index;
sub check_dependency;
sub check_jar;
sub check_genomicsdb;
sub check_ref_fasta;
sub cmd;
sub cmd_STDOUT;
sub concat;
sub countFastq;
sub find_job_info;
sub gdb_history;
sub get_dir_filelist;
#sub get_execDir;
sub get_path_bwa;
sub get_path_gatk;
sub get_path_samtools;
sub get_timestamp;
sub get_tmpdir;
sub IUPAC2regexp;
sub IUPAC_RC;
sub logentry;
sub logentry_then_die;
sub read_fodn;
sub read_fofn;
sub read_listfile;
sub retrieve_contigs;
sub round;
sub script_time;
sub split_in_files;
sub stage_files;
sub stage_gdb;
sub tarball_dir;

#internal functions
sub _sub_info;


#######################################
### avgQual
# determine mean quality score for a fastq quality string
sub avgQual {
	my $qstr=shift;
	chomp($qstr);
	my $qstr_sum;
	$qstr_sum += (ord $_)-33 foreach split //, $qstr;
	my $len=length($qstr);
	my $qavg=round($qstr_sum/$len,1);
	return ($qavg,$len);
}


#######################################
### archive_gdb
# Create an archived copy of a GenomicsDB and remove original.
sub archive_gdb {
	my $source=shift;
	my $error;
	my $tarball;

	$source=~s!/+$!!;  # remove any trailing slashes
	my ($source_dir,$source_stub)=($source =~ m{^(.*?)/?([^/]+)$}); # split to path and final dir
	$source_dir="./" if (! $source_dir);

	#set last modified as timestamp
	my $fstat=stat($source);
	my $timestamp=get_timestamp($fstat->mtime);

	my $archive_dir="${source_dir}/${source_stub}-${timestamp}";

	eval { make_path($archive_dir) } or do { logentry_then_die("Could not create temporary directory $archive_dir for archiving.",1); };
	if(-e "${source}.history") {
		eval { copy("${source}.history",$archive_dir) } or do { logentry_then_die("Could not copy original genomicsDB history ${source}.history to $archive_dir for archiving.",1); };
	} else {
		logentry("No history file located for original GenomicsDB.  No history is included in the archive file.",1);
	}
	eval { move($source,"${archive_dir}/${source_stub}") } or do { logentry_then_die("Could not move original genomicsDB $source to $archive_dir for archiving.",1); };
	eval { $tarball=tarball_dir($archive_dir,$source_dir) } or do { logentry_then_die("Could not create tarball of $archive_dir ($@).",1); };
	eval { remove_tree($archive_dir) } or do { logentry_then_die("Could not delete temporary directory ($archive_dir).",1); };

	if (-e $tarball && ! -e $source) {
		logentry("The original genomicsDB $source has been successfully archived to $tarball."),4;
		return 1;
	} else {
		logentry_then_die("Something unexpected happenned while archiving $source, not clear if operation succeeded.",1);
	}
}


#######################################
### build_argument_list
# Build list of arguments from array
sub build_argument_list {
	my @args=@{(shift)};
	my $delim=shift;
	my $arg_out=" ";
	foreach my $arg (@args) {
		$arg_out.="$delim $arg ";
	}
	return $arg_out;
}


#######################################
### build_contig_list
# Builds a list of contigs from a fasta file
sub build_contig_list {
	my $fasta=shift;
	my $out_path=shift;
	my @contigs=retrieve_contigs($fasta);
	open(my $OUT,">",$out_path);
	foreach my $contig (@contigs) {
		print $OUT $contig."\n";
	}
	close($OUT);
}


#######################################
### build_fasta_dictionary
# Confirm existence of build fasta dictionary.
sub build_fasta_dictionary {
	my $fasta=shift;
	my $fasta_stub=shift;
	my $gatk=get_path_gatk();
	logentry("BUILDING REFERENCE FASTA INDEX $fasta.fai",4);
	my $sys=cmd("$gatk CreateSequenceDictionary --REFERENCE $fasta --OUTPUT $fasta_stub.dict","Building reference dictionary");
}


#######################################
### build_fasta_index
# Build fasta index.
sub build_fasta_index {
	my $fasta=shift;
	my $samtools=get_path_samtools();
	logentry("BUILDING REFERENCE FASTA INDEX $fasta.fai",4);
	my $sys=cmd("$samtools faidx $fasta","Building reference index");
}


#######################################
### build_fofn
# Builds fofn file from array of files
sub build_fofn {
	my @files=@{(shift)};
	my $fofn_out=shift;
	logentry("Writing file paths to FOFN $fofn_out",4);
	foreach my $file (@files) {
		unless(-e $file && -f $file) {
			logentry("File path $file does not exist.  Not including in FOFN $fofn_out",1);
		}
	}
	build_listfile(\@files,$fofn_out);
}


#######################################
### build_listfile
# Builds list file from array
sub build_listfile {
	my @items=@{(shift)};
	my $file_out=shift;
	logentry("Writing list to $file_out",5);
	open(my $LIST_OUT, ">", "$file_out") or logentry_then_die("Cannot create list file $file_out");
	foreach my $item (@items) {
			print $LIST_OUT "$item\n";
	}
	close($LIST_OUT);
}


#######################################
### build_vcf_index
# Check existence and build vcf indexes
sub build_vcf_index {
	_sub_info((caller(0))[3],\@_);
	my $vcf_file=shift;
	my $gatk=get_path_gatk();
	my $idx_file=$vcf_file.".idx";
	logentry("Checking vcf index ($idx_file) for $vcf_file",11);
	if(-f $vcf_file && $vcf_file=~/vcf$/) {
		unless(-f $idx_file) {
			logentry("VCF index not found for $vcf_file: Building.",4);
			my $sys=cmd("$gatk IndexFeatureFile --feature-file $vcf_file","Calling GATK to build VCF index for $vcf_file") or logentry_then_die("Building VCF index for $vcf_file.  Cannot create $idx_file");
		}
	}
	else {
		logentry_then_die("File $vcf_file does not exist or does not appear to be a vcf file.");
	}
	return $idx_file;
}


#######################################
### check_dependency
# Confirm existence of dependency in system PATH, identify version, log, and return path.
sub check_dependency {
	my $prog=shift;       # binary name
	my $ver_flag=shift;   # flag, etc. to use to get version
	my $ver_regexp=shift; # regxp to extract version
	my $err_msg=shift;    # optional error message information (e.g. where to get software)
	my $report_out=shift;	# If defined . . . print version info to log
	$err_msg="" if (! $err_msg);
	my $path=which $prog or logentry_then_die("External dependency '${prog}' not installed in system PATH. ${err_msg}\n");
	chomp($path);
	my $ver="unknown";
	if($report_out) {
		$ver=cmd("${path} ${ver_flag}","Checking version of ${prog}") if($ver_flag);
		chomp($ver);
		eval "\$ver=~$ver_regexp" if($ver_regexp);
		logentry("Dependency ${prog} found: ${path} (${ver})",2);
	}
	return $path;
}


#######################################
### check_genomicsdb
# Confirm existence of genomicsDB and validate that it seems to have a valid structure
# Also removes gendb prefix if present and returns new path
sub check_genomicsdb {
	my $in=shift;
	my $db_input=shift;
	my $err_msg;
	my $err_msg_pre;
	if($db_input) {
		$err_msg_pre="Directory ($in) specified as --in, ";
	}
	else {
		$err_msg_pre="GenomicsDB location specified in --genomicsDB ($in), ";
	}
	$in =~ s!^gendb://!!;
	if(-e $in && -d $in) {
	 	unless(-f "$in/callset.json") {
			my $err_msg="${err_msg_pre}, location exists but does not seem to be a valid genomicsDB.";
			$err_msg.="  If this is meant to be a new genomicsDB, the directory specified must not exist." if (! $db_input);
			logentry_then_die($err_msg);
		}
	}
	else {
		logentry_then_die("${err_msg_pre}, but does not exist or is not a directory.");
	}
	return $in;
}


#######################################
### check_ref_fasta
# Confirm existence of reference fasta file and its indexes.  Build if necessary.
sub check_ref_fasta {
	my $refFasta=shift;
	my $gatk=get_path_gatk();
	my $samtools=get_path_samtools();
	my $sys;
	if(-e $refFasta && -f $refFasta) {
		my ($file,$path)=fileparse($refFasta, qr/\.[^.]*$/);
		my $refFasta_stub=$path.$file;
		if(! -e $refFasta.".fai") {
			logentry("REFERENCE FASTA INDEX NOT FOUND.",4);
			build_fasta_index($refFasta);
		}
		if(! -e $refFasta_stub.".dict") {
			logentry("REFERENCE FASTA DICTIONARY NOT FOUND.",4);
			build_fasta_dictionary($refFasta,$refFasta_stub);
		}
	}
	else {
		logentry_then_die("Reference fasta file $refFasta not found.");
	}
}



#######################################
### check_jar
# Confirm existence of dependency in system PATH, identify version, log, and return path.
sub check_jar {
	my $prog=shift;
	my $path=shift;
	my $ver_flag=shift;
	my $ver_regexp=shift;
	my $err_msg=shift;
	my $ver;

	chomp($path);
	if(-e $path && -f $path) {
		$ver=cmd("${main::java} ${path} ${ver_flag}","Checking version of ${prog}");
		chomp($ver);
		eval "\$ver=~$ver_regexp";
		logentry("Dependency ${prog} found: ${path} (${ver})",2);
	}
	else {
		logentry_then_die("External dependency '${prog}' not installed at ${path}. ${err_msg}\n");
	}
	return $path;
}


#######################################
### cmd
# run system command and collect output and error states
sub cmd {
	my $cmd=shift;
	my $message=shift;
	my $no_die=shift;

	logentry("System call ($message): $cmd",7);

	my $sys=`$cmd 2>&1`;
	my $err=$?;

	if ($err) {
		logentry_then_die("While trying to run a shell command.  Details:\ncommand: $cmd\nError Code: $err\nError message:\n$sys");
	}
	return $sys;
}


#######################################
### cmd_STDOUT
# run system command and collect output and error states (when std out must be used separately)
sub cmd_STDOUT {
	my $cmd=shift;
	my $message=shift;

	logentry("System call ($message): $cmd",7);

	my $sys=`$cmd 2> err`;
	my $err=$?;
	if ($err) {
		$sys=`cat err`;
		unlink("err");
		logentry_then_die("While trying to run a shell command ($message).  Details:\ncommand: $cmd\nError Code: $err\nError message:\n$sys");
		return 1;
	}
	else {
		unlink("err");
		logentry("STDOUT ($message): $sys",9);
	}
	return $sys;
}


#######################################
### concat
# concatenate multiple files specified in array
sub concat {
	my $inDir=shift;		# input directory
	my $outDir=shift;		# output directory
	my $outFile=shift;		# output filename
	my @files=@{(shift)};	# array of filenames

	my $file_str=join(' ',map { "$inDir/$_" } @files);
	my $sys=cmd("cat $file_str > $outDir/$outFile","Make concatenated fastq file $outFile");
}


#######################################
### countFastq
# count sequences in fastq file
sub countFastq {
	my $path=shift;
	my $message=shift;
	my $sys;
	if($path =~ /gz$/) {
		$sys=cmd('expr $(zcat '.$path.'| wc -l) / 4',$message);
		#my $fh = IO::Zlib->new($path) or logentry_then_die("Can't open file $path: $!");
		#$count++ while <$fh>;
		#undef($fh);
	}
	else {
		$sys=cmd('expr $(cat '.$path.'| wc -l) / 4',$message);
		#open(FILE, "< $path") or logentry_then_die("Can't open file $path: $!");
		#$count++ while <FILE>;
		#close(FILE);
	}
	return $sys;
}


#######################################
### find_job_info
# Determine if script is running as a job and return details
sub find_job_info {
	my $info;
	#PBS PBS_JOBNAME PBS_JOBID PBS_TASKNUM
	if($ENV{'PBS_JOBNAME'}) {
		$info="PBS scheduler found.  Job ID: $ENV{'PBS_JOBID'}; Job Name: $ENV{'PBS_JOBNAME'}; Number of Procs: $ENV{'PBS_NP'};\n";
	}
	#SLURM SLURM_JOB_ID SLURM_JOB_NAME SLURM_NTASKS
	elsif($ENV{'SLURM_JOB_NAME'}) {
		$info="SLURM scheduler found.  Cluster: $ENV{'SLURM_CLUSTER_NAME'}; Job ID: $ENV{'SLURM_JOB_ID'}; Job Name: $ENV{'SLURM_JOB_NAME'}; Number of Procs: $ENV{'SLURM_CPUS_PER_TASK'};\n";
		if($main::verbose>4) {
			my $sys=`scontrol show job $ENV{'SLURM_JOB_ID'}`;
			chomp($sys);
			$info .= "Detailed SLURM Job Info:\n   ".$sys if($sys);
		}
	}
	#SGE JOB_NAME JOB_ID NSLOTS
	elsif($ENV{'JOB_NAME'}) {
		$info="SGE scheduler found.  Job ID: $ENV{'JOB_ID'}; Job Name: $ENV{'JOB_NAME'}; Number of Slots: $ENV{'NSLOTS'};\n";
	}
	#LSF LSB_JOBID LSB_JOBNAME
	elsif($ENV{'LSB_JOBNAME'}) {
		$info="LSF scheduler found.  Job ID: $ENV{'LSB_JOBID'}; Job Name: $ENV{'LSB_JOBNAME'};\n";
	}
	else {
		$info="No task scheduler information detected.\n"
	}
return $info;
}


#######################################
### gdb_history
# Make and/or append to a genomicsDB history file.
sub gdb_history {
	my $gdb=shift;
	my $log_path=shift;
	my $gdb_archived=shift;
	my @in=@{(shift)};

	$gdb =~ s!/+$!!;  # remove any trailing slashes

	my $sep;
	my $message="Adding Files ".join(", ", @in).".  Job logged at $log_path.";
	my $warn;

	if (-e "${gdb}.history") {
		$sep="Appending to genomicsDB $gdb: ";
	}
	elsif ($gdb_archived) {
		$warn="WARNING: Original GenomicsDB does not appear to have an associated history file.  This history file is incomplete!\n";
		logentry("Original GenomicsDB does not appear to have an associated history file.  Creating new history file: history chain will be incomplete.",2);
		$sep="Appending to genomicsDB $gdb: Adding files ";
	}
	else {
		$sep="Creating genomicsDB $gdb: ";
	}

	open(my $HIST,">> ${gdb}.history");
	print $HIST $warn if($warn);
	print $HIST POSIX::strftime("%m/%d/%Y %H:%M:%S $sep $message\n", localtime);
	close($HIST);
}


#######################################
### get_dir_filelist
# open directory and return list of file paths matching a regexp pattern (case insensitive)
sub get_dir_filelist {
	my $dirpath=shift;
	my $pattern=shift;
	opendir(my $DIR,$dirpath) or logentry_then_die("Could not open directory path: $dirpath");
	my @files=grep /$pattern/i, readdir($DIR);
	close($DIR);
	s/^/$dirpath\// for @files;  #prepend $in (dirpath) to each element
	return @files;
}


#######################################
### get_path_bwa
# get bwa PATH
sub get_path_bwa {
	my $report_out=shift;
	my $bwa_path=check_dependency("bwa",' 2>&1 |grep "Version"',"s/\r?\n/ | /g",$report_out);
	return $bwa_path;
}


#######################################
### get_path_gatk
# get samtools PATH
sub get_path_gatk {
	my $report_out=shift;
	my $gatk_path  = check_dependency("gatk","--version 2>&1 | tail -n +4","s/\r?\n/ | /g",$report_out);
	return $gatk_path;
}


#######################################
### get_path_samtools
# get samtools PATH
sub get_path_samtools {
	my $report_out=shift;
	my $samtools_path=check_dependency("samtools","--version","s/\r?\n/ | /g",$report_out);
	return $samtools_path;
}


#######################################
### get_timestamp
# get timestamp suitable for filename
# example get_timestamp(time); would get current time
sub get_timestamp {
		my $time=shift;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime($time);
    my $timestamp = sprintf ( "%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $timestamp;
}


#######################################
### get_tmpdir
# Determine location of temporary directory and setup folder
sub get_tmpdir {
	my $tmpdir=shift;

	unless($tmpdir) {
		$tmpdir="/tmp";
		if ($ENV{'REDREP_TMPDIR'}) {
			$tmpdir=$ENV{'REDREP_TMPDIR'};
		}
		elsif ($ENV{'TMPDIR'}) {
			$tmpdir=$ENV{'TMPDIR'};
		}
		elsif ($ENV{'TEMP'}) {
			$tmpdir=$ENV{'TEMP'};
		}
		elsif ($ENV{'TMP'}) {
			$tmpdir=$ENV{'TMP'};
		}
	}
	unless(-e $tmpdir && -d $tmpdir) {
		logentry_then_die("Temporary directory $tmpdir does not exist.")
	}

	my $randStr=String::Random->new;
	$tmpdir.="/".$randStr->randregex("\\w{20}")."/";
	mkdir($tmpdir) or logentry_then_die("Could not create temporary directory $tmpdir.");
	return $tmpdir;
}


#######################################
### get_execDir
# Determine the RedRep bin location
#sub get_execDir {
#	my $execDir;
#	if($ENV{'REDREPBIN'}) {
#		if(-e $ENV{'REDREPBIN'} and -d $ENV{'REDREPBIN'}) {
#			$execDir=$ENV{'REDREPBIN'};
#			logentry("REDREP Installation location ($execDir) determined from REDREPBIN Environment Variable",2);
#		}
#		else {
#			my $death_throw="RedRep Installation location indicated in REDREPBIN Environment Variable is not a valid system directory!  Please define valid installation location with command: export REDREPBIN='PATH_TO_REDREP_INSTALLATION'\n";
#			logentry_then_die($death_throw);
#		}
#	}
#	else {
#		my($filename, $dirs, $suffix) = fileparse(abs_path($0));
#		$execDir=$dirs;
#		logentry("REDREPBIN Environment Variable not found.  RedRep Installation location assumed from location of executing script ($execDir),2")
#	}
#	return $execDir;
#}


#######################################
### IUPAC2regexp
# Convert sequence with IUPAC approved degenerate bases to regexp compatible pattern and return.
sub IUPAC2regexp {
	my $seq_in=shift;
	chomp($seq_in);
	my $seq_out;
	foreach my $char (split //, $seq_in) {
		if($char =~ /a/i) {
			$char = "[anvhdmwrANVHDMWR]";
		}
		elsif($char =~ /c/i) {
			$char = "[cnvhbmsyCNVHBMSY]";
		}
		elsif($char =~ /t/i) {
			$char = "[tunhdbkwyTUNHDBKWY]";
		}
		elsif($char =~ /g/i) {
			$char = "[gnvdbksrGNVDBKSR]";
		}
		elsif($char =~ /u/i) {
			$char = "[tunhdbkwyTUNHDBKWY]";
		}
		elsif($char =~ /n/i) {
			$char = "[actugnvhdbmkwsyrACTUGNVHDBMKWSYR]";
		}
		elsif($char =~ /v/i) {
			$char = "[acgnvmsrACGNVMSR]";
		}
		elsif($char =~ /h/i) {
			$char = "[actunhmywACTUNHMYW]";
		}
		elsif($char =~ /d/i) {
			$char = "[atugndwkrATUGNDWKR]";
		}
		elsif($char =~ /b/i) {
			$char = "[ctugnbyksCTUGNBYKS]";
		}
		elsif($char =~ /m/i) {
			$char = "[acmvhnACMVHN]";
		}
		elsif($char =~ /k/i) {
			$char = "[tugdbknTUGKDBN]";
		}
		elsif($char =~ /w/i) {
			$char = "[atuwhdnATUWHDN]";
		}
		elsif($char =~ /s/i) {
			$char = "[cgsvbnCGSVBN]";
		}
		elsif($char =~ /y/i) {
			$char = "[ctuyhbnCTUYHBN]";
		}
		elsif($char =~ /r/i) {
			$char = "[agrvdnAGRVDN]";
		}
		else {
			logentry_then_die("sequence $seq_in from metadata file contains illegal non-IUPAC base: $char.");
		}
		$seq_out .= $char;
	}
	return $seq_out;
}


#######################################
### IUPAC_RC
# Reverse Complement sequence all IUPAC Bases allowed
sub IUPAC_RC {
		my $seq = shift;
		chomp($seq);
  	my $seq = reverse($seq);
    $seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $seq;
}


#######################################
### logentry
# Enter time stamped log entry
sub logentry {
	my $message=shift;	# message to print
	my $level=shift;		# verbosity level.  0 always printed, 1 printed if verbosity is set. (optional, default=0)
	$level=4 unless($level);
	my $sep;
	if($level==0) {
		$sep="ERROR:";
	}
	elsif($level==1){
		$sep="WARNING:";
		cluck($message);
	}
	elsif($level==2){
		$sep="INFO:";
	}
	else {
		$sep=">" x ($level -2);
	}
	chomp($message);
	print $main::LOG POSIX::strftime("%m/%d/%Y %H:%M:%S $sep $message\n", localtime) if($main::verbose >= $level);
}


#######################################
### logentry_then_die
# Enter time stamped log entry
sub logentry_then_die {
	my $message=shift;
	my $no_unstage=shift;  # if set don't remove or unstage temp/intermed files -- used when die is expected to be caught
	logentry(longmess($message), 0); #log error with Carp backtrace
	logentry("Error encountered moving intermediate files to $main::outDir", 0);
	logentry("Cleaning up temporary files", 0);
	unless($no_unstage) { #do nothing
		unless($main::no_stage_intermed) {  #no need to move as its already in outDir
			logentry("Moving intermediate result files to output directory ${main::outDir}/intermed for troubleshooting",4);
			dircopy(${main::intermed},"${main::outDir}/intermed") if(${main::intermed} ne "${main::outDir}/intermed");
		}
		if($main::tmp_in_outdir || $main::debug) { # remove tmp
			logentry("Cleaning up temporary files",4);
			remove_tree($main::tmpdir);
		}	else {
			logentry("Temporary files have been left behind on the execution node in $main::tmpdir for troubleshooting because the --debug flag was used.",1) if($main::debug);
			remove_tree($main::intermed); # no need to leave two copies of this folder
		}
	}
	confess($message); #die with backtrace to STDERR
}


#######################################
### read_fodn
# read entries from a FODN and return array of validated filenames
sub read_fodn {
	my $infile=shift;
	my $file_pattern=shift;
	my @dirs=read_listfile($infile);
	my @files;
	#logentry("Parsing files from FODN $infile",4);
	foreach my $dir (@dirs) {
		chomp;
		if(-e $dir && -d $dir) {
			logentry("FODN $infile: Parsing files from directory $dir",5);
			push(@files,get_dir_filelist($dir,$file_pattern));
		}
		else {
			logentry_then_die("Input FODN $infile contains an entry that is not a directory: $dir");
		}
	}
	return(@files);
}


#######################################
### read_fofn
# read entries from a FOFN and return array of validated filenames
sub read_fofn {
	my $infile=shift;
	#logentry("Parsing files from FOFN $infile",4);
	my @files=read_listfile($infile);
	foreach my $file (@files) {
		unless (-e $file && -f $file) {
			logentry_then_die("File ($file) specified in FOFN $infile does not exist!");
		}
	}
	return @files;
}


#######################################
### read_listfile
# read entries from a list file and return array of items
sub read_listfile {
	my $file=shift;
	my @items;
	open(my $LIST,"<",$file);
	while(<$LIST>) {
		chomp;
		push(@items,$_);
	}
	close($LIST);
	return @items;
}


#######################################
### retrieve_contigs
# Retrieves contigs from a fasta index and returns as an array.
sub retrieve_contigs {
	my $fasta=shift;
	my $fai_path="$fasta.fai";
	build_fasta_index($fasta) if(! -e $fai_path);
	my @contigs;
	open(my $FAI,"<","$fai_path");
	while(<$FAI>) {
		chomp;
		my @cols=split(/\t/);
		push(@contigs,$cols[0]);
	}
	close($FAI);
	return @contigs;
}


#######################################
### round
# Round float ($number) to $dec digits
sub round {
	my $number = shift || 0;
	my $dec = 10 ** (shift || 0);
	return int( $dec * $number + .5 * ($number <=> 0)) / $dec;
}


#######################################
### script_time
# Returns string with total time for script execution
sub script_time {
	my $runtime = time - $^T;
	my $t = Time::Seconds->new($runtime);
	my $total=$t->pretty;
	return $total;
}


#######################################
### split_in_files
# Split input files, directories, FOFNs, and/or FODNs
sub split_in_files {
	my $in=shift;							# comma separated input string
	my $file_pattern=shift;		# file match pattern (probably based on file extension)
	my $fofn_out=shift;				# output file path for combined.fofn (optional)
	my @in;
	if($in=~/\,/) {
		@in=split(/\,/, $in);
	}
	else {
		push(@in,$in);
	}
	my @files;
	foreach my $innie (@in) {
		chomp $innie;
		# DIRECTORY
		if(-d $innie) {
			logentry("Parsing files from directory $innie",4);
			push(@files,get_dir_filelist($innie,$file_pattern));
		}
		# FODN
		elsif (-f $innie && ($innie=~/.fodn$/i || $innie=~/.fodn.list$/i)) {
			logentry("Parsing directories from FODN $innie",4);
			push(@files,read_fodn($innie,$file_pattern));
		}
		# FOFN
		elsif (-f $innie && ($innie=~/.fofn$/i || $innie=~/.fofn.list$/i)) {
			logentry("Parsing files from FOFN $innie",4);
			push(@files,read_fofn($innie));
		}
		# FILENAME
		elsif (-f $innie && $innie=~/$file_pattern/i) {
			push(@files,$innie);
		}
		# ERROR
		else {
			pod2usage( -msg  => "Argument to -i ($innie), file/directory not found or a file does not have a proper file extension.\n", -exitval => 2) if (! $innie);
		}
	}
	if($fofn_out) {
		build_fofn(\@files,$fofn_out);
	}
	return(@files);
}


#######################################
### stage_files
# Stage files to a target location
sub stage_files {
	my $target = shift;
	my @stage_files=@{(shift)};
	my @new_stage_files;
	foreach my $file (@stage_files) {
		copy($file,$target) or logentry_then_die("File $file could not be copied to tmpdir $target.  Ensure file exists and that tmpdir has sufficient space.");
		my $filename=fileparse($file);
		push(@new_stage_files, $target."/".$filename);
	}
	return(@new_stage_files);
}


#######################################
### stage_gdb
# Stage genomicsDB (or any directory) to a target location
sub stage_gdb {
	my $source = shift;
	my $target_base = shift;
	$source =~ m{/?([^/]+)/?$};
	my $target=$target_base."/".$1;

	#If exists stage the directory . . . if not only return the staged path
	if(-e $source) {
		dircopy($source,$target) or logentry_then_die("Directory $source could not be copied to tmpdir $target.  Ensure directory exists and that tmpdir has sufficient space.");
	}
	else {
		logentry("GenomicsDB $source does not exist.  Setting staged path for new GenomicsDB: $target.",4);
	}
	return($target);
}


#######################################
### _sub_info
# Can be added at beginning of function to log information about function calls.
# Can be useful for debugging.  Must use --verbose=12 or higher
# expected usage: _sub_info((caller(0))[3],\@_);
sub _sub_info {
	my $sub=shift;
	my @arr=@{(shift)};
	no warnings;
	my $args=join(",",@arr);
	logentry("$sub($args)",12);
}


#######################################
### tarball_dir
# Create an archived copy of a directory
# Uses pigz to compress if present, otehrwise gzip
# system must have tar command
sub tarball_dir {
	my $source=shift;
	my $target_dir=shift;
	my $ncpu=${main::ncpu};
	my ($source_dir,$source_stub)=($source =~ m{^(.*?)/?([^/]+)/?$});  # *? makes the .* non-greedy
	my $error;
	my $gzip="gzip";
	my $pigz=check_dependency("pigz","--version");
	$gzip="${pigz} --best --recursive -p $ncpu" if($pigz && $ncpu>1);
	eval { my $sys=cmd("tar --use-compress-program='$gzip' -C $source_dir -cf ${target_dir}/${source_stub}.tar.gz $source","Compressing copy of $source",1) } or logentry_then_die($@);
	return "${target_dir}/${source_stub}.tar.gz";
}



__END__


#######################################
########### DOCUMENTATION #############
=pod

=head1 NAME

RedRep/utils.pm

=head1 DESCRIPTION

Core utility modules for RedRep package

=head1 VERSION

version 2.11

=head1 SYNOPSIS

 use RedRep::utils;    #Exports check_dependency cmd logentry logentry_then_die
 use RedRep::utils qw(avgQual check_dependency check_jar cmd cmd_STDOUT concat countFastq find_job_info IUPAC2regexp IUPAC_RC logentry logentry_then_die round);

=head1 MODULES

=head2 archive_gdb

=head2 avgQual

=head2 build_argument_list

=head2 build_vcf_index

=head2 check_dependency

=head2 check_genomicsdb

=head2 check_jar

=head2 cmd

=head2 cmd_STDOUT

=head2 concat

=head2 countFastq

=head2 get_dir_filelist

=head2 find_job_info

=head2 get_timestamp

=head2 get_tmpdir

=head2 IUPAC2regexp

=head2 IUPAC_RC

=head2 logentry

=head2 logentry_then_die

=head2 retrieve_contigs

=head2 script_time

=head2 stage_files

=head2 stage_gdb

=head2 round

=head1 VERSION HISTORY

=over 3

=item 2.1 - 2/1/2017: Added utils.pm

=item 2.11 - 10/9/2019: Added expanded logging options with verbosity settings added.  Last version with GATK 3.x compatibility.

=item 2.2 - 10/28/2019: Added GATK v4 compatibility.  Added genomicsDB compatibility and many other enhancements including parallelization, documentation, and logging.

=item 2.3 - 8/1/2020: Fixed disk space reporting bugs.  Added stage_gdb function.

=back

=head1 DEPENDENCIES

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
