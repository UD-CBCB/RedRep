## Configuration for RedRep
## Change to match your system

###########################
###  REQUIRED SETTINGS  ###
###########################

# Location of RedRep installation.  Main scripts (redrep-xxxx.pl) should be in a subdirectory named "bin",
# accessory scripts in a subdirectory named "bin/util", and libraries in a subdirectory named "lib".
REDREP_INSTALL=~/00-scripts/redrep/active-git/GATK4-compatibility/RedRep

# shared/network file system (NFS, lustre, etc)?  0=no, 1=yes.  If running on a cluster leave at 1.  If running on a stand-alone machine that never accesses external storage, may be set to 0.  If in doubt, leave at 1.
SHARED_FS=1


##########################
###  PERL -- REQUIRED  ###
##########################
# PERL must be installed on the system and located or symlinked at /usr/bin/perl.


########################################################################################
### SETTINGS BELOW CAN BE USED TO SETUP DEPENDENDIES IF NOT ALREADY IN SYSTEM PATH.  ###
########################################################################################

# Location of FASTQC
FASTQC_PATH=/usr/local/FastQC

# Location of cutadapt
CUTADAPT_PATH=/usr/local/cutadapt/bin

# Location of bwa
BWA_PATH=/usr/local/bwa

# Location of samtools
SAMTOOLS_PATH=/usr/local/samtools

# location of cd-hit
CDHIT_PATH=/usr/local/cd-hit-4.7/

# location of usearch
USEARCH_PATH=/usr/local/usearch/

# location of GATK (v4 or higher)
GATK_PATH=/usr/local/GATK/


####################################################################################
###  DEPRECATED SETTINGS -- RETAINED FOR COMPATIBILITY WITH OLDER RELEASES ONLY  ###
####################################################################################

# DEPRECATED Location of picard-tools jarfile (jarfile name should be included)
#PICARD_JAR_PATH=/usr/local/picard-tools/picard.jar

# DEPRECATED Location of GATK jarfile (jarfile name should be included)
#GATK_JAR_PATH=/usr/local/GATK3/GenomeAnalysisTK.jar

# Location of picard-tools jarfile directory (DEPRECATED in ver 2.1, only set with version 2.0 and previous)
#PICARD_DIR_PATH=


######################################################
######################################################
###  NO CHANGES SHOULD BE NEEDED BELOW THIS POINT  ###
######################################################
######################################################

export REDREPBIN=${REDREP_INSTALL}/bin
export REDREPLIB=${REDREP_INSTALL}/lib
export REDREPUTIL=${REDREP_INSTALL}/util

# ${REDREPUTIL}/check_perl_dependencies.pl

export PATH=${REDREPBIN}:${FASTQC_PATH}:${CUTADAPT_PATH}:${BWA_PATH}:${SAMTOOLS_PATH}:${CDHIT_PATH}:${USEARCH_PATH}:${GATK_PATH}:${PATH}

#export PICARDJAR=${PICARD_JAR_PATH}
#export PICARDJARS=${PICARD_DIR_PATH}
#export GATKJAR=${GATK_JAR_PATH}

# Set at 1 if files may be stored on a shared file system (NFS, lustre, and others). Prevents issues with GATK database commands.
export TILEDB_DISABLE_FILE_LOCKING=${SHARED_FS}
