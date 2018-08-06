#------------------------------------------------------------------------------#
#                                                                              #
#                                  Figure S5                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(wesanderson)
library(GenomicRanges)
library(tidyverse)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load code to compute spike-in normalization factor
source(here('helper_spikein_normalization_factor.R'))
# Load analysis functions
source(here('helper_analysis.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
wes_colors <- wes_palette(name = "Darjeeling1", n=8, type = "continuous")
wt_color <- 'black'
dot1_color <- wes_colors[1]
set1_color <- wes_colors[4]
set1dot1_color <- wes_colors[8]

#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Spike-in normalization factor                                                #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors

# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah8583', ignore_case=TRUE))]

# Use all WT replicates
WT_chip_rep_counts <- list(
  read_counts_per_chr[['ah119spikeb-062817_trimmed']],
  read_counts_per_chr[['ah119spikeb-100917_trimmed']],
  read_counts_per_chr[['ah119spike-chip_trimmed']])

WT_input_rep_counts <- list(
  read_counts_per_chr[['ah119spikea-062817_trimmed']],
  read_counts_per_chr[['ah119spikea-100917_trimmed']],
  read_counts_per_chr[['ah119spike-inp_trimmed']])

sinfs <- data.frame(
  sample=c('WT', rep('dot1∆', 2), rep('set1∆', 2), rep('dot1∆ set1∆', 2)),
  Red1=c(1,
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8104spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8104spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8104spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8104spike-inp_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8584spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8584spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8584spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8584spike-inp_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8583spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8583spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8583spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8583spike-inp_trimmed']]))))

# Order samples
sinfs$sample <- factor(
  sinfs$sample, levels=c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆'))

# Plot
color_values <- c(wt_color, dot1_color, set1_color, set1dot1_color)

ggplot(sinfs, aes(sample, Red1 * 100, colour=sample, fill=sample)) +
  geom_hline(aes(yintercept = 100), linetype = 3) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_point(size=1.5, alpha=1) +
  scale_colour_manual('', values=color_values, guide=FALSE) +
  scale_fill_manual('', values=color_values, guide=FALSE) +
  scale_x_discrete(labels=c(
    'Wild type', expression(paste(italic('dot1'), Delta)),
    expression(paste(italic('set1'), Delta)),
    expression(paste(italic('dot1'), Delta, italic('set1'), Delta)))) +
  ylim(0, 100) +
  labs(title = '', x = '', y = 'Red1 amount\n(% of wild type)')

#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# Average Red1 signal per chromosome                                           #
#------------------------------------------------------------------------------#
signal_tracks <- list(
  'WT'=import_bedGraph(WT_1st_pileup),
  'dot1∆'=import_bedGraph(dot1_1st_pileup),
  'set1∆'=import_bedGraph(set1_1st_pileup),
  'dot1∆ set1∆'=import_bedGraph(set1dot1_1st_pileup)
)

# Drop spike-in genome
signal_tracks <- lapply(
  signal_tracks, keepSeqlevels, 
  paste0('chr', as.roman(1:16), '_SK1'), pruning.mode="coarse")

avrgs <- lapply(signal_tracks, function(x) average_chr_signal(x)[[2]])

# Compute spike-in normalization factor
# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah8583', ignore_case=TRUE))]

# Use all WT replicates
WT_chip_rep_counts <- list(
  read_counts_per_chr[['ah119spikeb-062817_trimmed']],
  read_counts_per_chr[['ah119spikeb-100917_trimmed']],
  read_counts_per_chr[['ah119spike-chip_trimmed']])

WT_input_rep_counts <- list(
  read_counts_per_chr[['ah119spikea-062817_trimmed']],
  read_counts_per_chr[['ah119spikea-100917_trimmed']],
  read_counts_per_chr[['ah119spike-inp_trimmed']])

sinfs <- data.frame(
  sample=c('WT', rep('dot1∆', 2), rep('set1∆', 2), rep('dot1∆ set1∆', 2)),
  Red1=c(1,
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8104spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8104spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8104spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8104spike-inp_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8584spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8584spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8584spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8584spike-inp_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8583spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8583spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8583spike-chip_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8583spike-inp_trimmed']]))))

# Normalize by average signal and then spike-in factor
for (i in seq_along(signal_tracks)) {
  sample <- names(signal_tracks)[i]
  nf <- mean(sinfs[sinfs$sample == sample, 'Red1'])
  signal_tracks[[i]]$score <- signal_tracks[[i]]$score / avrgs[[sample]] * nf
}

# Compute average signal per chromosome using normalized signal
chr_signal <- lapply(signal_tracks, function(x) average_chr_signal(x)[[1]])

for (i in seq_along(chr_signal)) chr_signal[[i]]$strain <- names(chr_signal)[i]

chr_signal <- do.call('rbind', chr_signal)

# Add chromosome sizes
chr_lengths <- seqlengths(get_chr_coordinates())
chr_lengths <- data.frame(chr=names(chr_lengths), chr_len=chr_lengths)
chr_signal <- dplyr::full_join(
  chr_signal, subset(chr_lengths, str_detect(chr, '_SK1')), by=c('chr'='chr'))

# Order samples
chr_signal$strain <- factor(chr_signal$strain,
                            levels=c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆'))

facet_names <- list('WT'='Wild type',
                    'dot1∆'=expression(paste(italic('dot1'), Delta)),
                    'set1∆'=expression(paste(italic('dot1'), Delta)),
                    'dot1∆ set1∆'=expression(paste(italic('dot1'), Delta,
                                                   italic('set1'), Delta)))
facet_labeller <- function(variable, value){return(facet_names[value])}

ggplot(chr_signal, aes(chr_len / 10^6, avrg_signal, colour=strain)) +
  geom_point(stat="identity", size=2, alpha=0.8) +
  facet_wrap(~strain, ncol=4, labeller=facet_labeller) +
  scale_color_manual(
    '', values=c(wt_color, dot1_color, set1_color, set1dot1_color),
    guide=FALSE) +
  #scale_y_continuous(limits = c(0, 2.4)) +
  labs(title = '', x = 'Chromosome size (Mb)',
       y='Red1 occupancy\n(average per bp)') +
  geom_hline(yintercept=1, lty=3) +
  scale_x_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1.0, 1.5)) +
  ylim(0, 1.7) +
  theme(strip.background = element_blank())

#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Red1 fragment pileup on selected example chromosomes                         #
#------------------------------------------------------------------------------#
WT <- import_bedGraph(WT_1st_pileup)
dot1 <- import_bedGraph(dot1_1st_pileup)
set1 <- import_bedGraph(set1_1st_pileup)
set1dot1 <- import_bedGraph(set1dot1_1st_pileup)

### Normalize signal
# Compute genome-wide signal average
WT_avrg <- average_chr_signal(
  subset(WT, str_detect(WT@seqnames, '_SK1')))[[2]]
dot1_avrg <- average_chr_signal(
  subset(dot1, str_detect(dot1@seqnames, '_SK1')))[[2]]
set1_avrg <- average_chr_signal(
  subset(set1, str_detect(set1@seqnames, '_SK1')))[[2]]
set1dot1_avrg <- average_chr_signal(
  subset(set1dot1, str_detect(set1dot1@seqnames, '_SK1')))[[2]]


# Compute spike-in normalization factor
subset(names(read_counts_per_chr),
       str_detect(names(read_counts_per_chr), '8104'))

# Use all WT replicates
WT_chip_rep_counts <- list(
  read_counts_per_chr[['ah119spikeb-062817_trimmed']],
  read_counts_per_chr[['ah119spikeb-100917_trimmed']],
  read_counts_per_chr[['ah119spike-chip_trimmed']])

WT_input_rep_counts <- list(
  read_counts_per_chr[['ah119spikea-062817_trimmed']],
  read_counts_per_chr[['ah119spikea-100917_trimmed']],
  read_counts_per_chr[['ah119spike-inp_trimmed']])

nfs <- data.frame(
  sample=c('dot1∆', 'set1∆', 'dot1∆ set1∆'),
  nf=c(suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8104spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8104spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8584spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8584spikea-100917_trimmed']])),
         suppressMessages(
           spikein_normalization_factor_from_counts(
             ref_chip_counts=WT_chip_rep_counts,
             ref_input_counts=WT_input_rep_counts,
             test_chip_counts=read_counts_per_chr[[
               'ah8583spikeb-100917_trimmed']],
             test_input_counts=read_counts_per_chr[[
               'ah8583spikea-100917_trimmed']]))))

# Normalize by average signal and then spike-in factor
WT[str_detect(WT@seqnames, '_SK1')]$score <-
  WT[str_detect(WT@seqnames, '_SK1')]$score / WT_avrg

nf <- as.numeric(subset(nfs, sample == 'dot1∆', select='nf'))
dot1[str_detect(dot1@seqnames, '_SK1')]$score <-
  dot1[str_detect(dot1@seqnames, '_SK1')]$score / dot1_avrg * nf

nf <- as.numeric(subset(nfs, sample == 'set1∆', select='nf'))
set1[str_detect(set1@seqnames, '_SK1')]$score <-
  set1[str_detect(set1@seqnames, '_SK1')]$score / set1_avrg * nf

nf <- as.numeric(subset(nfs, sample == 'dot1∆ set1∆', select='nf'))
set1dot1[str_detect(set1dot1@seqnames, '_SK1')]$score <-
  set1dot1[str_detect(set1dot1@seqnames, '_SK1')]$score / set1dot1_avrg * nf


make_df_for_plotting <- function(gr, sample_name='WT',
                                 normalize_by_genome_mean=FALSE, window_size=1,
                                 genome='SK1_S288C', chr='chrI_SK1',
                                 region=c(194000, 208000)) {
  t0 <- proc.time()[3]
  
  message('Compute binned score...')
  gr <- bedGraph_to_binned_score(
    gr, mean_normalize=normalize_by_genome_mean, tile_width=window_size,
    coord_genome=genome)
  
  message('Keep only selected chromosome...')
  gr <- keepSeqlevels(gr, chr, pruning.mode='coarse')
  
  if (!missing(region)) {
    message('Keep only selected region...')
    gr <- subset(gr, start >= region[1] & start <= region[2])
  }
  
  message('Convert to data frame for plotting...')
  gr <- data.frame(seqnames=seqnames(gr),
                   position=start(gr) + floor(width(gr) / 2),
                   signal=gr$binned_score)
  
  message('Add sample name to data frame...')
  gr$strain <- sample_name
  
  message('~~~')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  gr
}

WT_chrI <- make_df_for_plotting(
  WT, sample_name='WT', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')
dot1_chrI <- make_df_for_plotting(
  dot1, sample_name='dot1∆', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')
set1_chrI <- make_df_for_plotting(
  set1, sample_name='set1∆', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')
set1dot1_chrI <- make_df_for_plotting(
  set1dot1, sample_name='dot1∆ set1∆', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')

df_chrI <- rbind(WT_chrI, dot1_chrI, set1_chrI, set1dot1_chrI)

WT_chrII <- make_df_for_plotting(
  WT, sample_name='WT', normalize_by_genome_mean=FALSE,
  window_size=200, chr='chrII_SK1')
dot1_chrII <- make_df_for_plotting(
  dot1, sample_name='dot1∆', normalize_by_genome_mean=FALSE,
  window_size=200, chr='chrII_SK1')
set1_chrII <- make_df_for_plotting(
  set1, sample_name='set1∆', normalize_by_genome_mean=FALSE,
  window_size=200, chr='chrII_SK1')
set1dot1_chrII <- make_df_for_plotting(
  set1dot1, sample_name='dot1∆ set1∆', normalize_by_genome_mean=FALSE,
  window_size=200, chr='chrII_SK1')

df_chrII <- rbind(WT_chrII, dot1_chrII, set1_chrII, set1dot1_chrII)

# Order strains
df_chrI$strain = factor(
  df_chrI$strain, levels=c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆'))
df_chrII$strain = factor(
  df_chrII$strain, levels=c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆'))

# Plot
cen_midpoint <- function(chr, ref_genome='SK1_S288C', y_coord=0) {
  cen <- get_chr_coordinates(genome=ref_genome)
  
  if (!missing(chr)) {
    # Keep only required chr
    cen <- cen[cen@seqnames == chr]
    message('   Drop all chromosomes except ', chr, '...')
    gr <- keepSeqlevels(cen, chr, pruning.mode="coarse")
  }
  
  start <- cen@ranges@start
  half_width <- cen@ranges@width / 2
  cen_mid <- round(start + half_width)
  
  data.frame(seqnames = cen@seqnames, cen_mid = cen_mid, y = y_coord)
}

cen_mid <- cen_midpoint(chr='chrI_SK1', ref_genome='SK1_S288C', y_coord=0)

color_values <- c(wt_color, dot1_color, set1_color, set1dot1_color)

ggplot(df_chrI, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual(values=color_values) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.75, size=0.25) +
  geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
  annotate("text", x=cen_mid$cen_mid / 1000, y=-1.5, label="CEN",
           size=2, colour='black') +
  facet_grid(strain ~ .) +
  scale_fill_manual(values=color_values) +
  labs(title='',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  ylim(-2, 14) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

cen_mid <- cen_midpoint(chr='chrII_SK1', ref_genome='SK1_S288C', y_coord=0)

annotation_text <- data.frame(
  position=750000, signal=10,
  strain=factor(c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆'),
                levels=c('WT', 'dot1∆', 'set1∆', 'dot1∆ set1∆')))

ggplot(df_chrII, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual('', values=color_values) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.75, size=0.25) +
  geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
  annotate("text", x=cen_mid$cen_mid / 1000, y=-1.5, label="CEN",
           size=2, colour='black') +
  facet_grid(strain ~ .) +
  geom_text(data=annotation_text,
            label=c("paste('Wild ', type)",
                    "paste(italic('dot1'), Delta)", 
                    "paste(italic('set1'), Delta)",
                    "paste(italic('dot1'), Delta, italic(' set1'), Delta)"),
            parse=TRUE, size=3, colour=color_values) +
  scale_fill_manual('', values=color_values) +
  labs(title='',
       x=paste0('Position on chr II (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  ylim(-2, 14) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )
