package RedRep::Utils;
use strict;
use lib $ENV{'REDREPLIB'};
use Exporter qw(import);
use Pod::Usage;
use POSIX qw(strftime);
use Cwd qw(abs_path cwd);
use File::Basename qw(fileparse dirname);
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);
use File::Which qw(which);
#use IO::Zlib;

our @EXPORT    = qw(check_dependency cmd dir_get_files find_job_info find_tmpdir logentry logentry_then_die stage_files);
our @EXPORT_OK = qw(avgQual check_jar cmd_STDOUT concat countFastq find_job_info IUPAC2regexp IUPAC_RC round);


sub avgQual;
sub check_dependency;
sub check_jar;
sub cmd;
sub cmd_STDOUT;
sub concat;
sub countFastq;
sub dir_get_files;
sub find_job_info;
#sub get_execDir;
sub IUPAC2regexp;
sub IUPAC_RC;
sub logentry;
sub logentry_then_die;
sub round;
sub stage_files;


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
### check_dependency
# Confirm existence of dependency in system PATH, identify version, log, and return path.
sub check_dependency {
	my $prog=shift;       # binary name
	my $ver_flag=shift;   # flag, etc. to use to get version
	my $ver_regexp=shift; # regxp to extract version
	my $err_msg=shift;    # optional error message information (e.g. where to get software)
	$err_msg="" if (! $err_msg);
	my $path=which $prog || logentry_then_die("ERROR: External dependency '${prog}' not installed in system PATH. ${err_msg}\n");
	chomp($path);
	my $ver="unknown";
	$ver=cmd("${path} ${ver_flag}","Checking version of ${prog}") if($ver_flag);
	chomp($ver);
	eval "\$ver=~$ver_regexp" if($ver_regexp);
	logentry("Dependency ${prog} found: ${path} (${ver})",0);
	return $path;
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
		logentry("Dependency ${prog} found: ${path} (${ver})",0);
	}
	else {
		logentry_then_die("ERROR: External dependency '${prog}' not installed at ${path}. ${err_msg}\n");
	}
	return $path;
}


#######################################
### cmd
# run system command and collect output and error states
sub cmd {
	my $cmd=shift;
	my $message=shift;

	logentry("System call ($message): $cmd",1);

	my $sys=`$cmd 2>&1`;
	my $err=$?;
	if ($err) {
		logentry_then_die("ERROR while trying to run a shell command.  Details:\ncommand: $cmd\nError Code: $err\nError message:\n$sys");
	}
	else {
		logentry("STDOUT ($message): $sys",4);
	}
	return $sys;
}


#######################################
### cmd_STDOUT
# run system command and collect output and error states (when std out must be used separately)
sub cmd_STDOUT {
	my $cmd=shift;
	my $message=shift;

	logentry("System call ($message): $cmd",1);

	my $sys=`$cmd 2> err`;
	my $err=$?;
	if ($err) {
		$sys=`cat err`;
		unlink("err");
		logentry_then_die("ERROR while trying to run a shell command ($message).  Details:\ncommand: $cmd\nError Code: $err\nError message:\n$sys");
		return 1;
	}
	else {
		unlink("err");
		logentry("STDOUT ($message): $sys",4);
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
		#my $fh = IO::Zlib->new($path) or logentry_then_die("ERROR: Can't open file $path: $!");
		#$count++ while <$fh>;
		#undef($fh);
	}
	else {
		$sys=cmd('expr $(cat '.$path.'| wc -l) / 4',$message);
		#open(FILE, "< $path") or logentry_then_die("ERROR: Can't open file $path: $!");
		#$count++ while <FILE>;
		#close(FILE);
	}
	return $sys;
}


#######################################
### dir_get_files
# open directory and return list of file paths matching a regexp pattern (case insensitive)
sub dir_get_files {
	my $dirpath=shift;
	my $pattern=shift;
	opendir(my $DIR,$dirpath) || logentry_then_die("ERROR: Could not open directory path: $dirpath");
	my @tmpfiles=grep /\$pattern$/i, readdir($DIR);
	close($DIR);
	s/^/$dirpath\// for @tmpfiles;  #prepend $in (dirpath) to each element
	return @tmpfiles;
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
		if($main::verbose>1) {
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
### find_tmpdir
# Determine dafault location of temporary directory
sub find_tmpdir {
	my $tmpdir='/tmp';
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
#			logentry("REDREP Installation location ($execDir) determined from REDREPBIN Environment Variable",0);
#		}
#		else {
#			my $death_throw="ERROR: RedRep Installation location indicated in REDREPBIN Environment Variable is not a valid system directory!  Please define valid installation location with command: export REDREPBIN='PATH_TO_REDREP_INSTALLATION'\n";
#			logentry_then_die($death_throw);
#		}
#	}
#	else {
#		my($filename, $dirs, $suffix) = fileparse(abs_path($0));
#		$execDir=$dirs;
#		logentry("WARNING: REDREPBIN Environment Variable not found.  RedRep Installation location assumed from location of executing script ($execDir),0")
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
			logentry_then_die("ERROR: sequence $seq_in from metadata file contains illegal non-IUPAC base: $char.");
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
	$level=0 unless($level);
	my $sep=">" x ($level + 1);
	chomp($message);
	print $main::LOG POSIX::strftime("%m/%d/%Y %H:%M:%S $sep $message\n", localtime) if($main::verbose >= $level);
}


#######################################
### logentry_then_die
# Enter time stamped log entry
sub logentry_then_die {
	my $message=shift;
	logentry($message,0);
	die($message);
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
### stage_files
# Stage files to a target location
sub stage_files {
	my $target = shift;
	my @files=shift;
	my @newFiles;
	foreach my $file (@files) {
		copy($file,$target) || logentry_then_die("ERROR: File $file could not be copied to tmpdir $target.  Ensure file exists and that tmpdir has sufficient space.");
		my $filename=fileparse($file);
		push(@newFiles, $filename);
	}
	return(@newFiles);
}


__END__

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

=head1 VERSION HISTORY

=over 3

=item 2.1 = 2/1/2017: Added utils.pm

=item 2.11 = 10/9/2019: Added expanded logging options with verbosity settings added.

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

=item Exporter

=item File::Basename

=item File::Copy

=item File::Path

=item Getopt::Long

=item Pod::Usage

=item POSIX

=item Sys::Hostname

=back

=head1 MODULES

=head2 avgQual

=head2 check_dependency

=head2 check_jar

=head2 cmd

=head2 cmd_STDOUT

=head2 concat

=head2 countFastq

=head2 find_job_info

=head2 IUPAC2regexp

=head2 IUPAC_RC

=head2 logentry

=head2 logentry_then_die

=head2 round

=head1 AUTHOR

Written by Shawn Polson, University of Delaware

=head1 REPORTING BUGS

Report bugs to polson@udel.edu

=head1 COPYRIGHT

Copyright 2012-2019 Shawn Polson, Randall Wisser, Keith Hopper.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Please acknowledge author and affiliation in published work arising from this script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut
