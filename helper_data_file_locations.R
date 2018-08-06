# This file contains paths to all data files. It is used by the other files as #
# the centralized source of required file paths.                               #

# File locations need to point to the actual files.                            #

library(here)
library(stringr)

# Lrge fileData locations
bedGraph_dir <- '~/Documents/Datasets/Spike-in_bedGraph'
average_reps <- file.path(bedGraph_dir, 'Ave_Reps')

pileup_dir <- file.path('~/Documents/Datasets/Spike-in_Bedtools_Pileup')

#------------------------------------------------------------------------------#
#                            MACS2 fragment pileup                             #
#                               bedGraph files                                 #
#------------------------------------------------------------------------------#

# Non-spiked
# (data from:
#      Sun et al. Elife 2015,
#      Markowitz et al. PLoS Genet 2017)
nonspiked_rep_pileup <- list(
  'WT'=file.path(
    average_reps,
    'Red1-wildtype-29-34-199-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz'),
  'red1_ycs4S'=file.path(
    average_reps,
    'Red1-red1_ycs4S-AH7011-83-96-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
)

# Spiked
WT_1st_pileup <- file.path(
  bedGraph_dir, 'AH119spikein-060717_trimmed_S288c_SK1_Yue_PM_FE.bdg')
WT_2nd_pileup <- file.path(
  bedGraph_dir, 'AH119spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

WT_rep_pileup <- file.path(
  average_reps, 'Red1-wildtype-29-34-199-Reps-SK1_S288c_Yue-PM_SPMR_FE.bdg')

red1_pG162A_1st_pileup <- file.path(
  bedGraph_dir, 'AH9048spikein-060717_trimmed_S288c_SK1_Yue_PM_FE.bdg')
red1_pG162A_2nd_pileup <- file.path(
  bedGraph_dir, 'AH9048spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

hop1_1st_pileup <- file.path(
  bedGraph_dir, 'AH9120spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')
hop1_2nd_pileup <- file.path(
  bedGraph_dir, 'ah9120spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

hop1_rep_pileup <- file.path(
  average_reps, 'Red1-hop1-520-590-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg')

rec8_1st_pileup <- file.path(
  bedGraph_dir, 'ah8115spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

rec8_2nd_pileup <- file.path(
  bedGraph_dir, 'AH8115spike-Red1-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg')

dot1_1st_pileup <- file.path(
  bedGraph_dir, 'ah8104spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

dot1_2nd_pileup <- file.path(
  bedGraph_dir, 'AH8104spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

set1_1st_pileup <- file.path(
  bedGraph_dir, 'ah8584spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

set1_2nd_pileup <- file.path(
  bedGraph_dir, 'AH8584spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

set1dot1_1st_pileup <- file.path(
  bedGraph_dir, 'ah8583spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

set1dot1_2nd_pileup <- file.path(
  bedGraph_dir, 'AH8583spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')


gH2AX_pileup <- list(
  'WT'=file.path(
    bedGraph_dir,
    'AH119spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/RED1'=file.path(
    bedGraph_dir,
    'AH8218spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'red1/RED1'=file.path(
    bedGraph_dir,
    'AH8220spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/ycs4S'=file.path(
    bedGraph_dir,
    'AH7011spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/red1'=file.path(
    bedGraph_dir,
    'AH8219spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'spo11-YF/spo11-YF'=file.path(
    bedGraph_dir,
    'AH4206spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg')
)

#------------------------------------------------------------------------------#
#                                MACS2 called                                  #
#                                   peaks                                      #
#------------------------------------------------------------------------------#
peak_dir <- here('data/MACS2_peaks')

WT_spike_rep_peaks <- list(
  narrow=file.path(
    peak_dir,
    'Red1-WT-467-516-634-reps_S288C_SK1_Yue_PM_SPMR_peaks.narrowPeak'),
  broad=file.path(
    peak_dir,
    'Red1-WT-467-516-634-reps_S288C_SK1_Yue_PM_SPMR_peaks.broadPeak')
)

nonspiked_rep_peaks <- list(
  narrow=file.path(
    peak_dir,
    'Red1-wildtype-29-34-199-Reps-S288CYue-PM_SPMR_peaks.narrowPeak'),
  broad=file.path( # Data from Markowitz et al. PLoS Genet 2017
    peak_dir,
    'Red1-wildtype-29-34-199-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
)

#------------------------------------------------------------------------------#
#                                     GFF                                      #
#                                    files                                     #
#------------------------------------------------------------------------------#
gff_dir <- here('data/GFF')


SK1_gff_file <- file.path(gff_dir, 'SK1.all_feature.gff')
S288C_gff_file <- file.path(gff_dir, 'S288c.all_feature.gff')


#------------------------------------------------------------------------------#
#                               Read counts per                                #
#                            chromosome and genome                             #
#------------------------------------------------------------------------------#
read_count_dir <- here('data/Read_counts_per_chr/')

read_counts_per_chr <- as.list(sapply(list.files(read_count_dir),
                                      function(x) file.path(read_count_dir, x)))
# Make element names more readable
names(read_counts_per_chr) <- sapply(names(read_counts_per_chr),
                                     function(x) unlist(str_split(x, '_'))[4])

# Were the reads trimmed?
for (i in seq_along(read_counts_per_chr)) {
  names(read_counts_per_chr)[i] <- 
    ifelse(str_detect(read_counts_per_chr[[i]], 'trimmed'),
           yes=paste(names(read_counts_per_chr)[i], 'trimmed', sep='_'),
           no=names(read_counts_per_chr)[i])
}


#------------------------------------------------------------------------------#
#                            Aligned reads versus                              #
#                              sequencing depth                                #
#------------------------------------------------------------------------------#
Subsample_aligned_read_counts <- here(
  'data/Spike-in_subsample_aligned_read_counts/')


#------------------------------------------------------------------------------#
#                             Bedtools read pileup                             #
#                               bedGraph files                                 #
#------------------------------------------------------------------------------#
bedtools_pileups <- as.list(sapply(list.files(pileup_dir),
                                   function(x) file.path(pileup_dir, x)))
# Make element names more readable
names(bedtools_pileups) <- sapply(names(bedtools_pileups),
                                  function(x) unlist(str_split(x, '_'))[4])

#------------------------------------------------------------------------------#
#                                     SNP                                      #
#                                     list                                     #
#------------------------------------------------------------------------------#
SNPs <- here('data/S288c_v_SK1.snp')
