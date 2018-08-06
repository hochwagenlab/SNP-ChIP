# This file contains paths to all data files. It is used by the other files as #
# the centralized source of required file paths.                               #

library(here)
library(stringr)

#------------------------------------------------------------------------------#
#                            MACS2 fragment pileup                             #
#                               bedGraph files                                 #
#------------------------------------------------------------------------------#
base_dir <- '/Volumes/LabShare/HTGenomics/HiSeqOutputs'
base_dir <- '~/Desktop/SNP-ChIP_paper_data'

# Non-spiked
average_reps_SK1 <- file.path(base_dir, 'AveReps_SK1Yue_MACS2_FE')

nonspiked_rep_pileup <- list(
  'WT'=file.path(
    average_reps_SK1, 'Red1-wildtype-29-34-199-Reps-SK1Yue-B3W3-MACS2',
    'Red1-wildtype-29-34-199-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz'),
  'red1_ycs4S'=file.path( # Data from Markowitz et al. PLoS Genet 2017
    average_reps_SK1, 'Red1-red1_ycs4S-AH7011-83-96-Reps-SK1Yue-B3W3-MACS2',
    'Red1-red1_ycs4S-AH7011-83-96-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
)

# Spiked
first_exp <- file.path(base_dir, '2017-06-28_HLYHHAFXX_repeated_in_2017-08/')
second_exp <- file.path(base_dir, '2017-10-12_HNKY2AFXX')
third_exp <- file.path(base_dir, '2017-11-15_HVJH5AFXX')
fourth_exp <- file.path(base_dir, '2017-12-20_HW373AFXX')
fifth_exp <- file.path(base_dir, '2018-01-16_HWKGGAFXX_Red1_dosage_series')
average_reps <- file.path(base_dir, 'AveReps_SK1_S288c_Yue_hybrid_MACS2_FE')

WT_1st_pileup <- file.path(
  first_exp,
  'AH119spikein-060717_trimmed_S288c_SK1_Yue_PM',
  'AH119spikein-060717_trimmed_S288c_SK1_Yue_PM_FE.bdg')
WT_2nd_pileup <- file.path(
  second_exp,
  'AH119spikein-100917_trimmed_S288c_SK1_Yue_PM',
  'AH119spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

WT_rep_pileup <- file.path(
  average_reps,
  'Red1-wildtype-29-34-199-Reps-SK1_S288c_Yue-PM_SPMR',
  'Red1-wildtype-29-34-199-Reps-SK1_S288c_Yue-PM_SPMR_FE.bdg')

red1_pG162A_1st_pileup <- file.path(
  first_exp,
  'AH9048spikein-060717_trimmed_S288c_SK1_Yue_PM',
  'AH9048spikein-060717_trimmed_S288c_SK1_Yue_PM_FE.bdg')
red1_pG162A_2nd_pileup <- file.path(
  second_exp,
  'AH9048spikein-100917_trimmed_S288c_SK1_Yue_PM',
  'AH9048spikein-100917_trimmed_S288c_SK1_Yue_PM_FE.bdg')

hop1_1st_pileup <- file.path(
  second_exp,
  'AH9120spikein-100917_S288c_SK1_Yue_PM',
  'AH9120spikein-100917_S288c_SK1_Yue_PM_FE.bdg')
hop1_2nd_pileup <- file.path(
  third_exp,
  'ah9120spikein_11082017_S288c_SK1_Yue_PM',
  'ah9120spikein_11082017_S288c_SK1_Yue_PM_FE.bdg')

hop1_rep_pileup <- file.path(
  average_reps,
  'Red1-hop1-520-590-reps_S288C_SK1_Yue_PM_SPMR',
  'Red1-hop1-520-590-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg')

rec8_1st_pileup <- file.path(
  third_exp,
  'ah8115spikein_11082017_trimmed_S288c_SK1_Yue_PM',
  'ah8115spikein_11082017_trimmed_S288c_SK1_Yue_PM_FE.bdg')

rec8_2nd_pileup <- file.path(
  fifth_exp,
  'AH8115spike-Red1-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
  'AH8115spike-Red1-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg')

gH2AX_base_path <- file.path(base_dir,
                             '2018-01-16_HWKGGAFXX_Red1_dosage_series')
gH2AX_pileup <- list(
  'WT'=file.path(
  gH2AX_base_path, 'AH119spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
  'AH119spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/RED1'=file.path(
    gH2AX_base_path, 'AH8218spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
    'AH8218spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'red1/RED1'=file.path(
    gH2AX_base_path, 'AH8220spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
    'AH8220spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/ycs4S'=file.path(
    gH2AX_base_path, 'AH7011spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
    'AH7011spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'ycs4S/red1'=file.path(
    gH2AX_base_path, 'AH8219spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
    'AH8219spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg'),
  'spo11-YF/spo11-YF'=file.path(
    gH2AX_base_path, 'AH4206spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM',
    'AH4206spike-gh2ax-ChIP_011618_trimmed_S288c_SK1_Yue_PM_FE.bdg')
)

#------------------------------------------------------------------------------#
#                                MACS2 called                                  #
#                                   peaks                                      #
#------------------------------------------------------------------------------#
WT_spike_rep_peaks <- list(
  narrow=file.path(
    'data',
    'Red1-WT-467-516-634-reps_S288C_SK1_Yue_PM_SPMR_peaks.narrowPeak'),
  broad=file.path(
    'data',
    'Red1-WT-467-516-634-reps_S288C_SK1_Yue_PM_SPMR_peaks.broadPeak')
)

nonspiked_rep_peaks <- list(
  narrow=file.path(
    'data',
    'Red1-wildtype-29-34-199-Reps-S288CYue-PM_SPMR_peaks.narrowPeak'),
  broad=file.path( # Data from Markowitz et al. PLoS Genet 2017
    'data',
    'Red1-wildtype-29-34-199-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
)

#------------------------------------------------------------------------------#
#                                     GFF                                      #
#                                    files                                     #
#------------------------------------------------------------------------------#
base_dir <- file.path('~/Documents/Work/2015-2018_PostDoc_NYU/LabWork',
                      'GenomeSequences/S288C_SK1_Yue_hybrid_genome')

SK1_gff_file <- file.path(base_dir, 'SK1.all_feature.gff')
S288C_gff_file <- file.path(base_dir, 'S288c.all_feature.gff')


#------------------------------------------------------------------------------#
#                               Read counts per                                #
#                            chromosome and genome                             #
#------------------------------------------------------------------------------#
base_dir <- here('data/Read_counts_per_chr/')

read_counts_per_chr <- as.list(sapply(list.files(base_dir),
                                      function(x) file.path(base_dir, x)))
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
base_dir <- file.path('~/Documents/Datasets/Spike-in_Bedtools_Pileup')

bedtools_pileups <- as.list(sapply(list.files(base_dir),
                                   function(x) file.path(base_dir, x)))
# Make element names more readable
names(bedtools_pileups) <- sapply(names(bedtools_pileups),
                                  function(x) unlist(str_split(x, '_'))[4])

#------------------------------------------------------------------------------#
#                                     SNP                                      #
#                                     list                                     #
#------------------------------------------------------------------------------#
SNPs <- here('data/S288c_v_SK1.snp')


#------------------------------------------------------------------------------#
#                              In silico spike-in                              #
#                                 simulation                                   #
#------------------------------------------------------------------------------#
# base_dir <- file.path('~/Google_Drive_NYU/LabShare_Luis/Large_files',
#                       '2018.03_Spike-in_simulation')
# 
# insilico_pileup <- file.path(
#   base_dir, 'Simulated_spike-in_YueSK1_S288C_PM_SPMR',
#   'Simulated_spike-in_YueSK1_S288C_PM_SPMR_FE.bdg')
# 
# insilico_narrow_peaks <- file.path(
#   base_dir, 'Simulated_spike-in_YueSK1_S288C_PM_SPMR',
#   'Simulated_spike-in_YueSK1_S288C_PM_SPMR_peaks.narrowPeak')
# 
# insilico_broad_peaks <- file.path(
#   base_dir, 'Simulated_spike-in_YueSK1_S288C_PM_SPMR',
#   'Simulated_spike-in_YueSK1_S288C_PM_SPMR_peaks.broadPeak')
# 
# insilico_ctrl_pileup <- file.path(
#   base_dir, 'AH119_Red1_chip-SK1Yue-B3W3-MACS2',
#   'AH119_Red1_chip-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
# 
# insilico_ctrl_narrow_peaks <- file.path(
#   base_dir, 'AH119_Red1_chip-SK1Yue-B3W3-MACS2',
#   'AH119_Red1_chip-SK1Yue-PM_B3W3_MACS2_peaks.narrowPeak')
# 
# insilico_ctrl_broad_peaks <- file.path(
#   base_dir, 'AH119_Red1_chip-SK1Yue-B3W3-MACS2',
#   'AH119_Red1_chip-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak')
