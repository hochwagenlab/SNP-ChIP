#! /bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=4:00:00
#SBATCH --job-name=Spike-in_subset_mapper
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lv38@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out
#SBATCH --error=/scratch/%u/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#

# Generate subsets of spike-in experiment raw reads, map to hybrid genome and
# quantify number of mapped reads to each of the SK1 and S288C genomes.
# The ultimate goal is to determine the robustness of the calculated spike-in
# normalization factor (SPINF) across sequencing depths.

### Argument options:
# EXPID             Custom ID for output files
# RUNDIR            Path to directory to run script and save output in
# FQ                Path to input fastq file
# NREADS            Number of reads in the subsample (in millions)
# GZIP [GZIP=T/?]   If GZIP=T FASTQ file is gzipped at the end
#                   (only done if input file was originally gzipped)

### EXAMPLE:
# sbatch --export \
# EXPID="ah119spikea-062817_2M",\
# RUNDIR="/scratch/lv38",\
# FQ="HLYHHAFXX_n01_ah119spikea-062817_trimmed.fq.gz",\
# NREADS=2,GZIP=T \
# Spike-in_subset_mapping.slurm

#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

function elapsed_time() {
    ENDTIME=$(date +%s)

    TIME=$(($ENDTIME - $1))
    if [ $TIME -lt 60 ]
    then
        echo "$TIME sec"
    elif [ $TIME -ge 60 ]  && [ $TIME -lt 3600 ]
    then
        echo "$(($TIME / 60)) min"
    else
        echo "$(($TIME / 60 / 60)) hr"
    fi
}

function check_arg() {
    if [ -z "$1" ]
    then
        echo " ERROR: Please provide values for all required arguments"
        exit 2
    fi
}

function sam_to_sorted_and_indexed_bam() {
    # Convert to SAM file to BAM format, sort, and index

     local BASE=${1%.sam}
     local BAM=${BASE}.bam
     local S_BAM=${BASE}_sorted.bam

     samtools view -bS $1 > $BAM
     samtools sort -o $S_BAM $BAM
     samtools index $S_BAM

     # Clean up
     rm $1
     rm $BAM

     echo "$S_BAM"
}

function check_gzip() {
    if file --mime-type "$1" | grep -q gzip$
    then
        gzip -d $1
        ZIP_FQ=
    fi
}
#------------------------------------------------------------------------------#
#                                  IO checks                                   #
#------------------------------------------------------------------------------#
# Check arguments
check_arg $EXPID
check_arg $RUNDIR 
check_arg $FQ
check_arg $NREADS
check_arg $GZIP

cd $RUNDIR

# Check input file and run dir
[ -f $FQ ] || { echo " ERROR: Could not find file $FQ"; exit 2; }
[ -d $RUNDIR ] || { echo " ERROR: Could not find directory $RUNDIR"; exit 2; }

#------------------------------------------------------------------------------#
#                                                                              #
#                                Run pipeline                                  #
#                                                                              #
#------------------------------------------------------------------------------#

STARTTIME=$(date +%s)
echo \
"------------------------------------------------------------------------------"
echo ">>>>> Started pipeline: $EXPID"
echo \
"------------------------------------------------------------------------------"
date

# Bowtie 1.2.0 does not accept gzipped files
# Unzip them if gzipped (set empty variable to later know to rezip them)
check_gzip $FQ

[ ! -z ${ZIP_FQ+x} ] && echo ">>>>> Unzipped FASTQ file:" && echo $FQ

#------------------------------------------------------------------------------#
#                               Subsample FASTQ                                #
#------------------------------------------------------------------------------#

module purge
module load seqtk/intel/1.2-r94

echo ">>>>> Get sample of $NREADS million reads"

#f gzipped, remove the gzip extension from the stored file name
[ ! -z ${ZIP_FQ+x} ] && FEXT=${FQ##*.} && FQ=${FQ%.$FEXT}

FEXT=${FQ##*.}
FQBASE=${FQ%.$FEXT}
SUBFQ=${FQBASE}_${NREADS}Mreads.${FEXT}
#seqtk sample -s42 $FQ $(($NREADS * 1000000)) > $SUBFQ
seqtk sample -s42 $FQ ${NREADS}000000 > ${SUBFQ}

#------------------------------------------------------------------------------#
#                  Align reads to reference genome with Bowtie                 #
#------------------------------------------------------------------------------#

echo ">>>>> Map reads with Bowtie"
module purge
module load bowtie/gnu/1.2.1.1

REFGEN='/home/lv38/Library/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue'
OUTSAM=${EXPID}_${NREADS}Mreads.sam

bowtie -q -m 1 -v 0 -p 8 -S \
    --un ${EXPID}_${NREADS}Mreads_PM-Unaligned.fastq \
    --max ${EXPID}_${NREADS}Mreads_PM-Max.fastq \
    $REFGEN $SUBFQ $OUTSAM


if [ "$GZIP" = T ]
then
    # Was the file originally gzipped?
    if [ ! -z ${ZIP_FQ+x} ]
    then
        echo ">>>>> gzip FASTQ file:"
        echo $FQ
        gzip $FQ
    fi
fi

# Convert to BAM, sort and index
# In the process, reset file variables to map files (to use in next step)
echo ">>>>> Convert SAM to sorted and indexed BAM"
module purge
module load samtools/intel/1.3.1

SBAM=$(sam_to_sorted_and_indexed_bam $OUTSAM)

#------------------------------------------------------------------------------#
#                           Quantify read alignment                            #
#------------------------------------------------------------------------------#
echo ">>>>> Quantify reads aligned to SK1 and S88C genomes"
module purge
module load samtools/intel/1.3.1

samtools idxstats $SBAM | cut -f 1,3 > stats_${SBAM%_sorted.bam}.txt


#------------------------------------------------------------------------------#
#                                  Clean up                                    #
#------------------------------------------------------------------------------#
echo ">>>>> Delete unnecessary files..."
rm $SUBFQ
rm $SBAM
rm ${SBAM}.bai
rm ${EXPID}_${NREADS}Mreads_PM-Unaligned.fastq
rm ${EXPID}_${NREADS}Mreads_PM-Max.fastq

#------------------------------------------------------------------------------#
ELAPSEDTIME=$(elapsed_time $STARTTIME)
echo
echo "-----"
echo "-----"
echo "Completed pipeline in $ELAPSEDTIME"
echo \
"------------------------------------------------------------------------------"

