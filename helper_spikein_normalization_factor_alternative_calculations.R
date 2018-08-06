# This file produces a table containing spike-in normalization factors
# calculated using the different methods for all available samples.

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
names(read_counts_per_chr)[
  str_detect(read_counts_per_chr, fixed('ah8584', ignore_case=TRUE))]

# List of samples ordered correctly
files_for_calculation <- list(
  AH9048_1=c(read_counts_per_chr[['ah9048spikeb-062817_trimmed']],
             read_counts_per_chr[['ah9048spikea-062817_trimmed']]),
  AH9048_2=c(read_counts_per_chr[['ah9048spikeb-100917_trimmed']],
             read_counts_per_chr[['ah9048spikea-100917_trimmed']]),
  AH8104_1=c(read_counts_per_chr[['ah8104spikeb-100917_trimmed']],
             read_counts_per_chr[['ah8104spikea-100917_trimmed']]),
  AH8104_2=c(read_counts_per_chr[['ah8104spike-chip_trimmed']],
             read_counts_per_chr[['ah8104spike-inp_trimmed']]),
  AH8584_1=c(read_counts_per_chr[['ah8584spikeb-100917_trimmed']],
             read_counts_per_chr[['ah8584spikea-100917_trimmed']]),
  AH8584_2=c(read_counts_per_chr[['ah8584spike-chip_trimmed']],
             read_counts_per_chr[['ah8584spike-inp_trimmed']]),
  AH8583_1=c(read_counts_per_chr[['ah8583spikeb-100917_trimmed']],
             read_counts_per_chr[['ah8583spikea-100917_trimmed']]),
  AH8583_2=c(read_counts_per_chr[['ah8583spike-chip_trimmed']],
             read_counts_per_chr[['ah8583spike-inp_trimmed']]),
  AH9120_1=c(read_counts_per_chr[['ah9120spikeb-100917_trimmed']],
             read_counts_per_chr[['ah9120spikea-100917_trimmed']]),
  AH9120_2=c(read_counts_per_chr[['ah9120spike-chip_trimmed']],
             read_counts_per_chr[['ah9120spike-inp_trimmed']]),
  AH8115_1=c(read_counts_per_chr[['ah8115spike-chip_trimmed']],
             read_counts_per_chr[['ah8115spike-inp_trimmed']]),
  AH8115_2=c(read_counts_per_chr[['ah8115spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah8115spike-red1-inp_trimmed']]),
  AH8151_1=c(read_counts_per_chr[['ah8151spike-chip_trimmed']],
             read_counts_per_chr[['ah8151spike-inp_trimmed']]),
  AH8151_2=c(read_counts_per_chr[['ah8151spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah8151spike-red1-inp_trimmed']]),
  AH8218_1=c(read_counts_per_chr[['ah8218spike-chip_trimmed']],
             read_counts_per_chr[['ah8218spike-inp_trimmed']]),
  AH8218_2=c(read_counts_per_chr[['ah8218spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah8218spike-red1-inp_trimmed']]),
  AH8220_1=c(read_counts_per_chr[['ah8220spike-chip_trimmed']],
             read_counts_per_chr[['ah8220spike-inp_trimmed']]),
  AH8220_2=c(read_counts_per_chr[['ah8220spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah8220spike-red1-inp_trimmed']]),
  AH7011_1=c(read_counts_per_chr[['ah7011spike-chip_trimmed']],
             read_counts_per_chr[['ah7011spike-inp_trimmed']]),
  AH7011_2=c(read_counts_per_chr[['ah7011spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah7011spike-red1-inp_trimmed']]),
  AH8219_1=c(read_counts_per_chr[['ah8219spike-chip_trimmed']],
             read_counts_per_chr[['ah8219spike-inp_trimmed']]),
  AH8219_2=c(read_counts_per_chr[['ah8219spike-red1-chip_trimmed']],
             read_counts_per_chr[['ah8219spike-red1-inp_trimmed']])
)

# Use all WT replicates
WT_chip_rep_counts <- list(
  read_counts_per_chr[['ah119spikeb-062817_trimmed']],
  read_counts_per_chr[['ah119spikeb-100917_trimmed']],
  read_counts_per_chr[['ah119spike-chip_trimmed']])

WT_input_rep_counts <- list(
  read_counts_per_chr[['ah119spikea-062817_trimmed']],
  read_counts_per_chr[['ah119spikea-100917_trimmed']],
  read_counts_per_chr[['ah119spike-inp_trimmed']])

sinf_using_counts <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_counts <- c(
    sinf_using_counts,
    spikein_normalization_factor_from_counts(
      ref_chip=WT_chip_rep_counts,
      ref_input=WT_input_rep_counts,
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
  str_detect(bedtools_pileups, fixed('ah9048', ignore_case=TRUE))]

# List of samples ordered correctly
files_for_calculation <- list(
  AH9048_1=c(bedtools_pileups[['ah9048spikeb-062817']],
             bedtools_pileups[['ah9048spikea-062817']]),
  AH9048_2=c(bedtools_pileups[['ah9048spikeb-100917']],
             bedtools_pileups[['ah9048spikea-100917']]),
  AH8104_1=c(bedtools_pileups[['ah8104spikeb-100917']],
             bedtools_pileups[['ah8104spikea-100917']]),
  AH8104_2=c(bedtools_pileups[['ah8104spike-chip']],
             bedtools_pileups[['ah8104spike-inp']]),
  AH8584_1=c(bedtools_pileups[['ah8584spikeb-100917']],
             bedtools_pileups[['ah8584spikea-100917']]),
  AH8584_2=c(bedtools_pileups[['ah8584spike-chip']],
             bedtools_pileups[['ah8584spike-inp']]),
  AH8583_1=c(bedtools_pileups[['ah8583spikeb-100917']],
             bedtools_pileups[['ah8583spikea-100917']]),
  AH8583_2=c(bedtools_pileups[['ah8583spike-chip']],
             bedtools_pileups[['ah8583spike-inp']]),
  AH9120_1=c(bedtools_pileups[['ah9120spikeb-100917']],
             bedtools_pileups[['ah9120spikea-100917']]),
  AH9120_2=c(bedtools_pileups[['ah9120spike-chip']],
             bedtools_pileups[['ah9120spike-inp']]),
  AH8115_1=c(bedtools_pileups[['ah8115spike-chip']],
             bedtools_pileups[['ah8115spike-inp']]),
  AH8115_2=c(bedtools_pileups[['ah8115spike-red1-chip']],
             bedtools_pileups[['ah8115spike-red1-inp']]),
  AH8151_1=c(bedtools_pileups[['ah8151spike-chip']],
             bedtools_pileups[['ah8151spike-inp']]),
  AH8151_2=c(bedtools_pileups[['ah8151spike-red1-chip']],
             bedtools_pileups[['ah8151spike-red1-inp']]),
  AH8218_1=c(bedtools_pileups[['ah8218spike-chip']],
             bedtools_pileups[['ah8218spike-inp']]),
  AH8218_2=c(bedtools_pileups[['ah8218spike-red1-chip']],
             bedtools_pileups[['ah8218spike-red1-inp']]),
  AH8220_1=c(bedtools_pileups[['ah8220spike-chip']],
             bedtools_pileups[['ah8220spike-inp']]),
  AH8220_2=c(bedtools_pileups[['ah8220spike-red1-chip']],
             bedtools_pileups[['ah8220spike-red1-inp']]),
  AH7011_1=c(bedtools_pileups[['ah7011spike-chip']],
             bedtools_pileups[['ah7011spike-inp']]),
  AH7011_2=c(bedtools_pileups[['ah7011spike-red1-chip']],
             bedtools_pileups[['ah7011spike-red1-inp']]),
  AH8219_1=c(bedtools_pileups[['ah8219spike-chip']],
             bedtools_pileups[['ah8219spike-inp']]),
  AH8219_2=c(bedtools_pileups[['ah8219spike-red1-chip']],
             bedtools_pileups[['ah8219spike-red1-inp']])
)

# Use all WT replicates
WT_chip_rep_pileups <- list(
  bedtools_pileups[['ah119spikeb-062817']],
  bedtools_pileups[['ah119spikeb-100917']],
  bedtools_pileups[['ah119spike-chip']])

WT_input_rep_pileups <- list(
  bedtools_pileups[['ah119spikea-062817']],
  bedtools_pileups[['ah119spikea-100917']],
  bedtools_pileups[['ah119spike-inp']])

t0 <- proc.time()[3]

sinf_using_pileups <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_pileups <- c(
    sinf_using_pileups,
    spike_in_normalization_factor_from_pileup(
      ref_chip=WT_chip_rep_pileups,
      ref_input=WT_input_rep_pileups,
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
t0 = proc.time()[3]
sinf_using_pileups_on_SNPs <- vector()
for (i in seq_along(files_for_calculation)) {
  sinf_using_pileups_on_SNPs <- c(
    sinf_using_pileups_on_SNPs,
    spike_in_normalization_factor_from_pileup(
      ref_chip=WT_chip_rep_pileups,
      ref_input=WT_input_rep_pileups,
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
peaks_SK1 <- import_bedGraph(
  here('data/Red1-wildtype-29-34-199-Reps-SK1Yue-PM_B3W3_MACS2_peaks.narrowPeak'))
peaks_S288C <- import_bedGraph(
  here('data/Red1-wildtype-29-34-199-Reps-S288CYue-PM_SPMR_peaks.narrowPeak'))

# Make unified peaks GRanges
peaks_S288C <- renameSeqlevels(
  peaks_S288C, paste0(GenomeInfoDb::seqlevels(peaks_S288C), '_S288C'))
peaks_SK1 <- renameSeqlevels(
  peaks_SK1, paste0(GenomeInfoDb::seqlevels(peaks_SK1), '_SK1'))
peaks <- c(peaks_SK1, peaks_S288C)

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
      ref_chip=WT_chip_rep_pileups,
      ref_input=WT_input_rep_pileups,
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

# write_csv(
#     sinfs,
#     here('data/spikein_normalization_factors_using-different_methods.csv'))
