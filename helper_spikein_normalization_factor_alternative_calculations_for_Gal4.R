# This file produces a table containing spike-in normalization factors
# calculated using the different methods for the Gal4 experiment samples.

library(here)
library(readr)
library(GenomicRanges)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load IO code
source(here('helper_io.R'))
# Load code to compute spike-in normalization factor
source(here('helper_spikein_normalization_factor.R'))

#------------------------------------------------------------------------------#
#                                   Using                                      #
#                                read counts                                   #
#------------------------------------------------------------------------------#
# Locate samples
# Locate samples
names(read_counts_per_chr)[
  str_detect(read_counts_per_chr, fixed('ah2243', ignore_case=TRUE))]

# List of samples ordered correctly
files_for_calculation <- list(
  '0'=c(read_counts_per_chr[['ah2243-0nm-estradiol-gal4-chip_trimmed']],
        read_counts_per_chr[['ah2243-0nm-estradiol-input_trimmed']]),
  '10'=c(read_counts_per_chr[['ah2243-10nm-estradiol-gal4-chip_trimmed']],
         read_counts_per_chr[['ah2243-10nm-estradiol-input_trimmed']]),
  '100'=c(read_counts_per_chr[['ah2243-100nm-estradiol-gal4-chip_trimmed']],
          read_counts_per_chr[['ah2243-100nm-estradiol-input_trimmed']]),
  '10000'=c(
    read_counts_per_chr[['ah2243-1000nm-estradiol-gal4-chip_trimmed']],
    read_counts_per_chr[['ah2243-1000nm-estradiol-input_trimmed']])
)

# Define reference (WT strain)
WT_chip_counts <- files_for_calculation[['0']][1]
WT_input_counts <- files_for_calculation[['0']][2]

sinf_using_counts <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_counts <- c(
    sinf_using_counts,
    spikein_normalization_factor_from_counts(
      ref_chip=WT_chip_counts,
      ref_input=WT_input_counts,
      test_chip=files_for_calculation[[i]][1],
      test_input=files_for_calculation[[i]][2]
    )
  )
}


#------------------------------------------------------------------------------#
#                                   Using                                      #
#                              bedtools pileup                                 #
#------------------------------------------------------------------------------#
# Locate samples
names(bedtools_pileups)[
  str_detect(bedtools_pileups, fixed('ah2243', ignore_case=TRUE))]

# List of samples ordered correctly
files_for_calculation <- list(
  '0'=c(bedtools_pileups[['ah2243-0nm-estradiol-gal4-chip']],
        bedtools_pileups[['ah2243-0nm-estradiol-input']]),
  '10'=c(bedtools_pileups[['ah2243-10nm-estradiol-gal4-chip']],
         bedtools_pileups[['ah2243-10nm-estradiol-input']]),
  '100'=c(bedtools_pileups[['ah2243-100nm-estradiol-gal4-chip']],
          bedtools_pileups[['ah2243-100nm-estradiol-input']]),
  '10000'=c(bedtools_pileups[['ah2243-1000nm-estradiol-gal4-chip']],
            bedtools_pileups[['ah2243-1000nm-estradiol-input']])
)

# Define reference (WT strain)
WT_chip_pileups <- files_for_calculation[['0']][1]
WT_input_pileups <- files_for_calculation[['0']][2]

t0 <- proc.time()[3]

sinf_using_pileups <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_pileups <- c(
    sinf_using_pileups,
    spike_in_normalization_factor_from_pileup(
      ref_chip=WT_chip_pileups,
      ref_input=WT_input_pileups,
      test_chip=files_for_calculation[[i]][1],
      test_input=files_for_calculation[[i]][2]
    )
  )
}
message()
message('---')
message('---')
message('Completed in ', elapsed_time(t0, proc.time()[3]))

#------------------------------------------------------------------------------#
#                                   Using                                      #
#                      bedtools pileup overlapping SNPs                        #
#------------------------------------------------------------------------------#
# Load list of SNPs
SNPs <- read_tsv(SNPs)
tail(SNPs)

# Convert to GRanges
SNPs <- with(SNPs, GRanges(chr, IRanges(position, position)))
SNPs

# Rerun analysis keeping only ranges overlapping SNPs
t0 <- proc.time()[3]
sinf_using_pileups_on_SNPs <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_pileups_on_SNPs <- c(
    sinf_using_pileups_on_SNPs,
    spike_in_normalization_factor_from_pileup(
      ref_chip=WT_chip_pileups,
      ref_input=WT_input_pileups,
      test_chip=files_for_calculation[[i]][1],
      test_input=files_for_calculation[[i]][2],
      positions_to_keep=SNPs 
    )
  )
}
message()
message('---')
message('---')
message('Completed in ', elapsed_time(t0, proc.time()[3]))


#------------------------------------------------------------------------------#
#                                   Using                                      #
#               bedtools pileup overlapping SNPs in Red1 peaks                 #
#------------------------------------------------------------------------------#
# Load Red1 peaks
peaks <- import_bedGraph(here(
  'data/AH2243_1000nm_estradiol_spike-in_06032018_trimmed_S288c_SK1_Yue_PM_peaks.narrowPeak'))

peaks <- import_bedGraph(
  '/Volumes/LabShare/HTGenomics/HiSeqOutputs/2018-03-06_HYJGYAFXX/AH2243_0nm_estradiol_spike-in_06032018_trimmed_S288c_SK1_Yue_PM/AH2243_0nm_estradiol_spike-in_06032018_trimmed_S288c_SK1_Yue_PM_peaks.narrowPeak')

sort(peaks)
sort(SNPs)

# Subset SNPs by overlap of Red1 peaks (use the subset to calculate SINFs)
SNPs_in_peaks <- subsetByOverlaps(SNPs, peaks)
message('Kept ', length(SNPs_in_peaks), ' SNPs (overlapping Red1 peaks)')

# Rerun analysis keeping only ranges overlapping SNPs within Red1 peaks
t0 = proc.time()[3]
sinf_using_pileups_on_SNPs_in_peaks <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_pileups_on_SNPs_in_peaks <- c(
    sinf_using_pileups_on_SNPs_in_peaks,
    spike_in_normalization_factor_from_pileup(
      ref_chip=WT_chip_pileups,
      ref_input=WT_input_pileups,
      test_chip=files_for_calculation[[i]][1],
      test_input=files_for_calculation[[i]][2],
      positions_to_keep=SNPs_in_peaks
    )
  )
}
message()
message('---')
message('---')
message('Completed in ', elapsed_time(t0, proc.time()[3]))


#------------------------------------------------------------------------------#
#                              Save final table                                #
#                                  to file                                     #
#------------------------------------------------------------------------------#

sinf_using_counts
sinf_using_pileups
sinf_using_pileups_on_SNPs
sinf_using_pileups_on_SNPs_in_peaks

sinfs <- rbind(
  data.frame(sample=names(files_for_calculation),
             sinf=sinf_using_counts,
             method='read_counts'),
  data.frame(sample=names(files_for_calculation),
             sinf=sinf_using_pileups,
             method='read_pileups'),
  data.frame(sample=names(files_for_calculation),
             sinf=sinf_using_pileups_on_SNPs,
             method='read_pileups_on_SNPs'),
  data.frame(sample=names(files_for_calculation),
             sinf=sinf_using_pileups_on_SNPs_in_peaks,
             method='read_pileups_on_SNPs_in_peaks')
)

write_csv(sinfs, here(
  'data/spikein_normalization_factors_using_different_methods_for_Gal4.csv'))
