#------------------------------------------------------------------------------#
#                              Analysis helper                                 #
#                                 functions                                    #
#------------------------------------------------------------------------------#

# Load paths to data files
source(here::here('helper_data_file_locations.R'))
# Load IO code
source(here::here('helper_io.R'))

library(GenomicRanges)
library(EnrichedHeatmap)

#' Add genome name to GRanges
#'
#' Adds genome name to \code{GRanges} object.
#' @param gr \code{GRanges} object. No default.
#' @param name Genome name \code{string}. No default.
#' @return Same \code{GRanges} object with a genome name.
#' @examples
#' \dontrun{
#' x <- add_genome_name_to_GR(WT, name='SK1_S288CYue')
#' }
#' @export
add_genome_name_to_GR <- function(gr, name) {
  number_seqs <- length(levels(gr@seqnames))
  gr@seqinfo@genome <- rep(name, number_seqs)
  
  gr
}


#' Get chromosome coordinates
#'
#' Returns chromosome length and centromere coordinates for an hybrid of the SK1
#' and S288C assemblies from Yue et al. 2017 (here named "SK1Yue"). Chromosome
#' lengths were calculated in the bash shell using commands like the following 
#' (SK1 example): \cr
#' \code{cat SK1.genome.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100)
#' "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'}
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1"}
#'   \item \code{"S288C"}
#'   \item \code{"SK1_S288C"}
#' }
#' Defaults to \code{"SK1_S288C"}.
#' @return SK1 coordinates as a \code{GRanges} object.
#' @examples
#' \dontrun{
#' x <- get_chr_coordinates(genome="SK1")
#' }
#' @export

get_chr_coordinates <- function(genome='SK1_S288C') {
  SK1_gff <- rtracklayer::import.gff(SK1_gff_file)
  S288C_gff <- rtracklayer::import.gff(S288C_gff_file)
  
  SK1_start <- SK1_gff[SK1_gff$type == 'centromere']@ranges@start
  SK1_end <- SK1_start + SK1_gff[SK1_gff$type == 'centromere']@ranges@width - 1
  
  S288C_start <- S288C_gff[S288C_gff$type == 'centromere']@ranges@start
  S288C_end <- 
    S288C_start + S288C_gff[S288C_gff$type == 'centromere']@ranges@width - 1
  
  # Coordinate df
  if (genome == 'SK1') {
    coord_table <- data.frame(
      "Chromosome" = paste0('chr', as.roman(1:16)),
      "Start" = SK1_start, "End" = SK1_end,
      "LenChr" =  c(228861, 829469, 340914, 1486921, 589812, 299318,
                    1080440, 542723, 449612, 753937, 690901, 1054145,
                    923535, 791982, 1053869, 946846))
  } else if (genome == 'S288C') {
    coord_table <- data.frame(
      "Chromosome" = paste0('chr', as.roman(1:16)),
      "Start" = S288C_start, "End" = S288C_end,
      "LenChr" = c(219929, 813597, 341580, 1566853, 583092, 271539,
                   1091538, 581049, 440036, 751611, 666862, 1075542,
                   930506, 777615, 1091343, 954457))
  } else if (genome == 'SK1_S288C') {
    coord_table <- data.frame(
      "Chromosome" = c(paste0('chr', as.roman(1:16), '_SK1'),
                       paste0('chr', as.roman(1:16), '_S288C')),
      "Start" = c(SK1_start, S288C_start),
      "End" = c(SK1_end, S288C_end),
      "LenChr" = c(
        c(228861, 829469, 340914, 1486921, 589812, 299318,
          1080440, 542723, 449612, 753937, 690901, 1054145,
          923535, 791982, 1053869, 946846),
        c(219929, 813597, 341580, 1566853, 583092, 271539,
          1091538, 581049, 440036, 751611, 666862, 1075542,
          930506, 777615, 1091343, 954457)))
  } else stop('"genome" argument must be one of ',
              '"SK1", "S288C" or "SK1_S288C".')
  
  # Convert to GRanges
  SK1_S288C <- with(coord_table,
                    GenomicRanges::GRanges(
                      Chromosome, IRanges::IRanges(Start + 1, End),
                      seqlengths=setNames(LenChr, Chromosome)))
  # Add genome name and return
  add_genome_name_to_GR(SK1_S288C, name=genome)
}


#' Normalize by genome average
#'
#' Computes average signal genome-wide and divides every signal value by
#' computed average.
#' @param gr \code{GRanges} object. No default.
#' @return Same \code{GRanges} object with average-normalized signal.
#' @examples
#' \dontrun{
#' x <- norm_by_genome_average(gr=WT)
#' }
#' @export
norm_by_genome_average <- function(gr){
  # Compute genome average
  message('   Genome-wide signal average: ', appendLF=FALSE)
  compute_avrg <- function(x) {
    (sum(GenomicRanges::width(x) * GenomicRanges::score(x)) /
       sum(GenomicRanges::width(x)))
  }
  avrg <- compute_avrg(gr)
  message(round(avrg, 3))
  
  # Divide each score value by genome average signal
  message('   Normalize signal...')
  gr$score <- gr$score / avrg
  
  message('   ---')
  message('   Done!')
  return(gr)
}


#' Convert bedGraph data to binned score data
#'
#' Computes average signal in defined genomic bins genome-wide.
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param mean_normalize Logical specifying whether to normalize signal by
#' genome-wide average signal. Defaults to \code{FALSE}.
#' @param tile_width Integer specifying the length (in bp) of the window
#' to average the signal in. Defaults to \code{100} (i.e. 100 kb).
#' @param coord_genome Character object specifying the genome version
#' (\code{get_chr_coordinates}'s namesake argument); accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1"}
#'   \item \code{"S288C"}
#'   \item \code{"SK1_S288C"}
#' }
#' Defaults to \code{"SK1_S288C"}.
#' @return A \code{GRanges} object contained the calculated inned score.
#' @examples
#' \dontrun{
#' bedGraph_to_binned_score(WT, dot1)
#'
#' bedGraph_to_binned_score(WT, dot1, 1000)
#' 
#' bedGraph_to_binned_score(WT, dot1, 100, coord_genome='SK1')
#' }
#' @export

bedGraph_to_binned_score <- function(gr, mean_normalize=FALSE, tile_width=100,
                                     coord_genome="SK1_S288C") {
  t0  <- proc.time()[3]
  
  if (!is(gr, "GRanges")) {
    stop('"gr" must be a GRanges object.', call. = FALSE)
  }
  
  if (!is(tile_width, "numeric")) {
    stop('"tile_width" must be numeric.', call. = FALSE)
  }
  
  # Make sure it is integer
  tile_width <- floor(tile_width)
  
  if (mean_normalize) {
    # Normalize ChIP signal by genome average
    message('   Normalize by genome-wide average signal...')
    gr <- norm_by_genome_average(gr)
  }
  
  ### Compute binned signal average in selected windows
  message('   Compute binned score (', tile_width, '-bp bins)...')
  # Get "seqlengths"
  genome_info <- get_chr_coordinates(genome=coord_genome)
  # Sort sequences and levels to make sure they match
  genome_info <- sortSeqlevels(genome_info)
  gr <- sort(sortSeqlevels(gr))
  
  # Add info to signal object
  suppressWarnings(
    seqlengths(gr) <- seqlengths(genome_info)
  )
  
  # Compute tiling windows
  bins <- tileGenome(seqlengths(gr), tilewidth=tile_width,
                     cut.last.tile.in.chrom=TRUE)
  
  # Get signal as "RleList"; the signal is stored in the "score" metadata column
  score <- coverage(gr, weight="score")
  
  # Compute average signal per tile
  bins <- binnedAverage(bins, score, "binned_score")
  
  message('---')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  bins
}


#' Compute correlation between signal track data for two samples
#'
#' Returns correlation between ChIP-seq signal of two samples. Useful to check
#' replicates or compare different samples.
#' 
#' @param signal_data_A Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param signal_data_B Second samples' signal track data in the same format as
#' \code{signal_data_A} object. No default.
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1"}
#'   \item \code{"S288C"}
#'   \item \code{"SK1_S288C"}
#' }
#' No default.
#' @param method Character object specifying which correlation coefficient is to
#' be computed. Will be the input to the \code{method} argument of function
#' \code{cor} of the \code{stats} package. Accepts one of (can be abbreviated):
#' \enumerate{
#'   \item \code{"pearson"}
#'   \item \code{"kendall"}
#'   \item \code{"spearman"}
#' }
#' Defaults to \code{"pearson"}.
#' @param window_size Integer specifying the constant width (in bp) of the
#' windows to average signal in. Defaults to \code{1}.
#' @return A floating point value corresponding to the signal correlation.
#' @examples
#' \dontrun{
#' signal_track_correlation(WT, dot1, genome = "SK1")
#' 
#' signal_track_correlation(WT_A, WT_B, genome = "SK1_S288C", method="spearman")
#' }
#' @export

signal_track_correlation <- function(signal_data_A, signal_data_B, genome,
                                     method="pearson", window_size=1) {
  t0  <- proc.time()[3]
  
  if (!is(signal_data_A, "GRanges") | !is(signal_data_B, "GRanges")) {
    stop('"signal_data_A" and "signal_data_B" must be GRanges objects.',
         call. = FALSE)
  }
  
  if (missing(genome)) stop('"genome" is a required argument.\n', call. = FALSE)
  
  # Get "seqlengths"
  message('Get genome coordinates...')
  genome_info <- get_chr_coordinates(genome=genome)
  
  # Sort sequences and levels to make sure they match
  genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
  grA <- sort(GenomeInfoDb::sortSeqlevels(signal_data_A))
  grB <- sort(GenomeInfoDb::sortSeqlevels(signal_data_B))
  
  # Add info to signal object
  suppressWarnings(
    GenomeInfoDb::seqlengths(grA) <- GenomeInfoDb::seqlengths(genome_info)
  )
  suppressWarnings(
    GenomeInfoDb::seqlengths(grB) <- GenomeInfoDb::seqlengths(genome_info)
  )
  
  # Compute 1-bp tiling windows
  message('Compute tiling windows...')
  bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(grA),
                                    tilewidth=window_size,
                                    cut.last.tile.in.chrom=TRUE)
  
  # Get signal as "RleList"; the signal is stored in the "score" metadata column
  message('Compute signal per bp (or tile)...')
  grA <- GenomicRanges::coverage(grA, weight="score")
  grB <- GenomicRanges::coverage(grB, weight="score")
  
  # Get signal per tile
  A_bins <- GenomicRanges::binnedAverage(bins, grA, "binned_score")
  B_bins <- GenomicRanges::binnedAverage(bins, grB, "binned_score")
  
  # Replace zeros by R's representation of missing data
  message('Replace zeros by "NA"s...')
  A_bins[A_bins$binned_score == 0]$binned_score <- NA
  B_bins[B_bins$binned_score == 0]$binned_score <- NA
  
  message('Compute ', method, ' correlation between samples: ',
          deparse(substitute(signal_data_A)), ' and ',
          deparse(substitute(signal_data_B)), '...')
  corr <- cor(x=A_bins$binned_score, y=B_bins$binned_score,
              use="complete.obs", method=method)
  
  message('---')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  
  list(corr,
       data.frame('X'=A_bins$binned_score, 'Y'=B_bins$binned_score))
}


#' Average signal genome-wide and on each chromosome
#'
#' Computes average signal (in the \code{score} \code{GRanges} metadata column)
#' genome-wide and on each chromosome (each individual sequence determined by
#' \code{seqnames}).
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). No default.
#' @return List with two elements:
#' \enumerate{
#'   \item \code{seq_avrg} Average \code{score} on each sequence (dataframe)
#'   \item \code{genome_avrg} Average \code{score} genome-wide (named vector)
#' }
#' @examples
#' \dontrun{
#' average_chr_signal(GRanges_object)
#' }
#' @export
average_chr_signal <- function(gr){
  # IO checks
  if (!is(gr, "GRanges")) stop('input must be a GRanges object.')
  if (!"score" %in% names(GenomicRanges::mcols(gr))) {
    stop(deparse(substitute(gr)), ' does not have a "score" metadata column.')
  }
  
  message('Computing average signal...')
  avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x),
                           na.rm = T)
                       / sum(GenomicRanges::width(x)))
  genome_avrg <- avrg(gr)
  seq_avrg <- sapply(GenomicRanges::split(gr, GenomicRanges::seqnames(gr)),
                     avrg)
  
  # Convert to dataframe
  seq_avrg <- data.frame(chr=names(seq_avrg), avrg_signal=seq_avrg,
                         row.names=NULL, stringsAsFactors = F)
  
  message('Done!')
  return(list("seq_avrg"=seq_avrg, "genome_avrg"=genome_avrg))
}


#' Collect signal from telomeres for all chromosome arms
#'
#' Pulls out the ChIP signal starting at chromosome ends inward up to a
#' specified distance.
#' @param signal_data Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). To load wiggle and bedGraph data
#' run \code{\link{import_wiggle}} and \code{\link{import_bedGraph}},
#' respectively. No default.
#' @param length_to_collect Integer specifying the length (in bp) of the region
#' to collect signal for, starting form the telomeres. Defaults to 100000 (i.e.
#' 100 kb).
#' @param averaging_window Integer specifying the length (in bp) of the window
#' to average the signal in.
#' Defaults to 1 (i.e. 1 bp: no window averaging).
#' @return An R data frame containing 32 rows (one for each chromosome arm) and
#' a number of columns dependent on the specified \code{length_to_collect}. The
#' following columns are included:
#' \enumerate{
#'   \item \code{chr} The chromosome number
#'   \item \code{arm} The chromosome arm; one of "L" and "R"
#'   \item \code{size-cat} The chromosome size category; one of "small" and
#'   "large"
#'   \item \code{t1:tn} Signal columns for each position; number depends on the
#'   the specified collection length (\code{n = length_to_collect})
#' }
#' @examples
#' \dontrun{
#' signal_from_telomeres2(WT)
#' 
#' signal_from_telomeres2(WT, length_to_collect = 50000, averaging_window = 100)
#' }
#' @export


signal_from_telomeres <- function(signal_data, length_to_collect=100000,
                                  averaging_window=1) {
  t0  <- proc.time()[3]
  
  # IO checks
  if (!is(signal_data, "GRanges")) {
    stop('"signal_data" must be a GRanges object.', call. = FALSE)
  }
  
  if (!is(length_to_collect, "numeric")) {
    stop('"length_to_collect" must be numeric.', call. = FALSE)
  }
  
  coord_table <- get_chr_coordinates()
  
  message('Making GRanges object of subtelomeric regions...')
  # Make sure it is integer
  length_to_collect <- floor(length_to_collect)
  
  # Left arms
  left_arm <- coord_table
  
  GenomicRanges::start(left_arm) <- 1
  GenomicRanges::end(left_arm) <- length_to_collect
  GenomicRanges::strand(left_arm) <- '+'
  
  # Right arms
  right_arm <- coord_table
  ends <- GenomeInfoDb::seqlengths(right_arm)
  GenomicRanges::end(right_arm) <- ends
  GenomicRanges::start(right_arm) <- (ends - length_to_collect)
  GenomicRanges::strand(right_arm) <- '-'
  
  # Concatenate arms
  telomeres <- c(left_arm, right_arm)
  
  # Compute signal at each gene using package EnrichedHeatmap
  message('Collecting signal...')
  number_of_windows <- floor(length_to_collect / averaging_window)
  mat <- normalizeToMatrix(signal_data, telomeres, value_column="score",
                           mean_mode="absolute", extend=0, k=number_of_windows,
                           empty_value=NA, smooth=FALSE, target_ratio=1)
  
  # Prepare final data frame
  df <- as.data.frame(mat)
  size_cat <- ifelse(as.character(telomeres@seqnames) %in% c('chrI_SK1',
                                                             'chrIII_SK1',
                                                             'chrVI_SK1',
                                                             'chrI_S288C',
                                                             'chrIII_S288C',
                                                             'chrVI_S288C'),
                     'small', 'large')
  
  df <- cbind(chr=as.character(telomeres@seqnames),
              arm=c(rep('L', 16), rep('R', 16)),
              size_cat=size_cat,
              df)
  
  message('---')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  return(df)
}