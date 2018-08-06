# This file contains the functions to compute the spike-in normalization factor

library(stringr)
library(readr)

# Load IO code
source(here('helper_io.R'))

#' Compute spike-in normalization factor from total read counts
#'
#' Computes spike-in normalization factor between two spiked-in samples using
#' total counts of aligned reads. Inputs paths to text files containing counts
#' of aligned reads per chromosome of a hybrid SK1:S288C genome.
#' @param ref_chip_counts Either a single or a list of paths to reference ChIP
#' samples' read counts file. No default.
#' @param ref_input_counts Either a single or a list of paths to reference input
#' samples' read counts file. No default.
#' @param test_chip_counts Either a single or a list of paths to test ChIP
#' samples' read counts file. No default.
#' @param test_input_counts Either a single or a list of paths to test input
#' samples' read counts file. No default.
#' @param return_counts Logical indicating whether to return the computed read
#' counts instead of the normalization factor. Defaults to \code{FALSE}.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spikein_normalization_factor_from_counts(
#'      ref_chip_counts='Counts_AH119_chip.txt',
#'      ref_input_counts='Counts_AH119_input.txt',
#'      test_chip_counts='Counts_AH8104_chip.txt',
#'      test_input_counts='Counts_AH8104_input.txt')
#'
#' spikein_normalization_factor_from_counts(
#'     ref_chip_counts=list('Counts_AH119_chip_1.txt',
#'                          'Counts_AH119_chip_2.txt',
#'                          'Counts_AH119_chip_3.txt'),
#'     ref_input_counts=list('Counts_AH119_inp_1.txt',
#'                           'Counts_AH119_inp_2.txt',
#'                           'Counts_AH119_inp_3.txt'),
#'     test_chip_counts='Counts_AH8104_chip.txt',
#'     test_input_counts='Counts_AH8104_input.txt')
#' }
#' @export
spikein_normalization_factor_from_counts <- function(
  ref_chip_counts, ref_input_counts, test_chip_counts, test_input_counts,
  return_counts=FALSE) {
  
  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)
  
  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }
  
  # Print files to read to console
  message('>>> Read alignment count files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    
  
  message()
  # Read files into tibble in list
  tables <- list()
  for (i in seq_along(files)) {
    tables[[i]] <- sapply(files[[i]], FUN=read_tsv, col_names=F,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(tables) <- names(files)
  
  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- list()
  for (i in seq_along(tables)) {
    counts[[i]] <- sapply(tables[[i]], FUN=sum_per_genome,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(counts) <- names(tables)
  
  # Add-up counts for replicates (results in nested lists)
  for (i in seq_along(counts)) {
    if (length(counts[[i]]) > 1) {
      total <- counts[[i]][[1]]
      for (j in 2:length(counts[[i]])) {
        total <- total + counts[[i]][[j]]
      }
      counts[[i]] <- total
    } else counts[[i]] <- unlist(counts[[i]])
  }
  
  if (return_counts) {
    message('---')
    message('Done!')
    return(counts)
  }
  
  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)
  
  message('---')
  message('Done!')
  
  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
  # Compute sum of reads aligned to each genome
  S288C <- sum(
    df[apply(df, 1, function(x) str_detect(x[1],'_S288C')), 2])
  SK1 <- sum(
    df[apply(df, 1, function(x) str_detect(x[1], '_SK1')), 2])
  
  # Print result to console
  message('  S288C: ', formatC(S288C, big.mark=",",
                               drop0trailing=TRUE, format="f"))
  message('  SK1: ', formatC(SK1, big.mark=",",
                             drop0trailing=TRUE, format="f"))
  message('      ', round(S288C * 100 / (SK1 + S288C), 1), '% spike-in reads')
  
  # Return result as named vector
  c('S288C'=S288C, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
  # Compute Q values
  Q_ctrl_input <- ctrl_input['S288C'] / ctrl_input['SK1']
  Q_ctrl_chip <- ctrl_chip['S288C'] / ctrl_chip['SK1']

  Q_test_input <- test_input['S288C'] / test_input['SK1']
  Q_test_chip <- test_chip['S288C'] / test_chip['SK1']

  # Compute normalization factors
  a_ctrl <- Q_ctrl_input / Q_ctrl_chip
  a_test <- Q_test_input / Q_test_chip

  # Return reference strain-centric normalization factor
  a_test/ a_ctrl
}


#' Compute spike-in normalization factor from read pileups
#'
#' Computes spike-in normalization factor between two spiked-in samples using
#' read coverage per genomic position (read pileups). Inputs paths to bedgraph
#' files containing pileups of reads aligned to hybrid SK1:S288C genome.
#' @param fun Function used to summarise coverage values per genome bacgkground.
#' Defaults to \code{mean}.
#' @param ref_chip Either a single or a list of paths to reference ChIP samples'
#' read pileup bedgraph file. No default.
#' @param ref_input Either a single or a list of paths to reference input
#' samples' read pileup bedgraph file. No default.
#' @param test_chip Either a single or a list of paths to test ChIP samples'
#' read pileup bedgraph file. No default.
#' @param test_input Either a single or a list of paths to test input samples'
#' read pileup bedgraph file. No default.
#' @param positions_to_keep \code{GRanges} object containing genomic positions 
#' to keep for the analysis. If provided, pileup \code{GRanges} will be subset
#' by overlaps with these positions. No default.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spike_in_normalization_factor_from_pileup(
#'      ref_chip='Counts_AH119_chip.bdg',
#'      ref_input='Counts_AH119_input.bdg',
#'      test_chip='Counts_AH8104_chip.bdg',
#'      test_input='Counts_AH8104_input.bdg')
#'      
#' spike_in_normalization_factor_from_pileup(
#'      ref_chip='Counts_AH119_chip.bdg',
#'      ref_input='Counts_AH119_input.bdg',
#'      test_chip='Counts_AH8104_chip.bdg',
#'      test_input='Counts_AH8104_input.bdg',
#'      positions_to_keep=SNPs)
#'
#' spike_in_normalization_factor_from_pileup(
#'     ref_chip_counts=list('Counts_AH119_chip_1.bdg',
#'                          'Counts_AH119_chip_2.bdg',
#'                          'Counts_AH119_chip_3.bdg'),
#'     ref_input_counts=list('Counts_AH119_inp_1.bdg',
#'                           'Counts_AH119_inp_2.bdg',
#'                           'Counts_AH119_inp_3.bdg'),
#'     test_chip_counts='Counts_AH8104_chip.bdg',
#'     test_input_counts='Counts_AH8104_input.bdg')
#' }
#' @export

spike_in_normalization_factor_from_pileup <- function(
  fun=mean, ref_chip, ref_input, test_chip, test_input, positions_to_keep) {
  t0 <- proc.time()[3]
  
  # Put paths in list
  files <- list(ref_chip=ref_chip, ref_input=ref_input,
                test_chip=test_chip, test_input=test_input)
  
  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }
  
  # Print files to read to console
  message('>>> Load bedgraph files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    
  
  message()
  # Read files into GRanges
  grs <- list()
  for (i in seq_along(files)) {
    grs[[i]] <- suppressMessages(sapply(files[[i]], FUN=import_bedGraph,
                                        simplify=FALSE, USE.NAMES=TRUE))
  }
  names(grs) <- names(files)
  
  # Concatenate replicate GRanges objects
  for (i in seq_along(grs)) {
    grs[[i]] <- do.call('c', grs[[i]])
  }
  
  if (!missing(positions_to_keep)) {
    message('>>> Subset by overlaps with ',
            deparse(substitute(positions_to_keep)))
    grs <- lapply(grs, function(gr) subsetByOverlaps(gr, positions_to_keep))
  }
  
  get_genome_scores <- function(data, genome_label) {
    subset(data, str_detect(data@seqnames, genome_label))$score
  }
  
  message('>>> Compute ', deparse(substitute(fun)), '...')
  ref_chip_SK1 <- fun(get_genome_scores(grs[['ref_chip']],
                                        genome_label='_SK1'))
  ref_input_SK1 <- fun(get_genome_scores(grs[['ref_input']],
                                         genome_label='_SK1'))
  test_chip_SK1 <- fun(get_genome_scores(grs[['test_chip']],
                                         genome_label='_SK1'))
  test_input_SK1 <- fun(get_genome_scores(grs[['test_input']],
                                          genome_label='_SK1'))
  
  ref_chip_S288C <- fun(get_genome_scores(grs[['ref_chip']],
                                          genome_label='_S288C'))
  ref_input_S288C <- fun(get_genome_scores(grs[['ref_input']],
                                           genome_label='_S288C'))
  test_chip_S288C <- fun(get_genome_scores(grs[['test_chip']],
                                           genome_label='_S288C'))
  test_input_S288C <- fun(get_genome_scores(grs[['test_input']],
                                            genome_label='_S288C'))
  
  message('>>> Compute spike-in normalization factor...')
  # Compute Q values
  Q_ctrl_input <- ref_input_S288C / ref_input_SK1
  Q_ctrl_chip <- ref_chip_S288C / ref_chip_SK1
  
  Q_test_input <- test_input_S288C / test_input_SK1
  Q_test_chip <- test_chip_S288C / test_chip_SK1
  
  # Compute normalization factors
  nf_ctrl <- Q_ctrl_input / Q_ctrl_chip
  nf_test <- Q_test_input / Q_test_chip
  
  # Return reference strain-centric normalization factor
  ref_centric_nf <- nf_test/ nf_ctrl
  
  message('---')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  
  ref_centric_nf
}
