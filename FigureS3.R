#------------------------------------------------------------------------------#
#                                                                              #
#                                  Figure S1                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(wesanderson)
library(tidyverse)
library(GenomicRanges)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load IO code
source(here('helper_io.R'))
# Load analysis functions
source(here('helper_analysis.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
wes_colors <- wes_palette(name = "Darjeeling1", n=8, type = "continuous")
spiked_color <- wes_colors[2]
ctrl_color <- 'black'

#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# In silico spike-in versus control: Red1 fragment pileup example              #
# EDIT: Drop in silico experiment; just use spiked versus non-spiked low-Red1  #
#------------------------------------------------------------------------------#
# In the zoom in add SNPs as vertical red bars
# WT_spike_pileup <- import_bedGraph(WT_rep_pileup)
WT_spike_pileup <- import_bedGraph(WT_1st_pileup)
WT_nonspiked_pileup <- import_bedGraph(nonspiked_rep_pileup[['WT']])


# Discard S288C (spike-in) chromosomes and standardize names
WT_spike_pileup <- keepSeqlevels(WT_spike_pileup,
                                 paste0('chr', as.roman(1:16), '_SK1'),
                                 pruning.mode='coarse')

renaming_map_names <- paste0('chr', as.roman(1:16), '_SK1')
renaming_map <- paste0('chr', as.roman(1:16))
names(renaming_map) <- renaming_map_names
WT_spike_pileup <- renameSeqlevels(WT_spike_pileup, renaming_map)

make_df_for_plotting <- function(gr, sample_name='WT',
                                 normalize_by_genome_mean=FALSE, window_size=1,
                                 genome='SK1_S228C', chr='chrI_SK1',
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

WT_spike_pileup_chrI <- make_df_for_plotting(
  WT_spike_pileup, sample_name='Spiked',
  normalize_by_genome_mean=FALSE, window_size=100, genome='SK1', chr='chrI')
WT_nonspiked_pileup_chrI <- make_df_for_plotting(
  WT_nonspiked_pileup, sample_name='Non-spiked',
  normalize_by_genome_mean=FALSE, window_size=100, genome='SK1', chr='chrI')

WT_spike_chrI_region <- make_df_for_plotting(
  WT_spike_pileup, sample_name='Spiked',
  normalize_by_genome_mean=FALSE, window_size=1, genome='SK1', chr='chrI',
  region=c(53000, 57000))
WT_nonspiked_chrI_region <- make_df_for_plotting(
  WT_nonspiked_pileup, sample_name='Non-spiked',
  normalize_by_genome_mean=FALSE, window_size=1, genome='SK1', chr='chrI',
  region=c(53000, 57000))

df_chrI <- rbind(WT_spike_pileup_chrI, WT_nonspiked_pileup_chrI)
df_chrI_region <- rbind(WT_spike_chrI_region, WT_nonspiked_chrI_region)

# Order strains
df_chrI$strain = factor(
  df_chrI$strain, levels=c('Spiked', 'Non-spiked'))
df_chrI_region$strain = factor(
  df_chrI_region$strain, levels=c('Spiked', 'Non-spiked'))

### Plot
# Whole chromosome
ggplot(df_chrI, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual(values=c(spiked_color, ctrl_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
  annotate('text', x=23, y=7, label=unique(df_chrI$strain),
           size=3, colour=c(spiked_color, ctrl_color)) +
  scale_fill_manual(values=c(spiked_color, ctrl_color)) +
  labs(title='', x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none", strip.background = element_blank()) +
  labs(y=paste0('Red1 occupancy')) +
  ylim(0, 8) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

# Annotate SNPs in region plot
SNPs <- read_tsv(SNPs)
SNPs_chrI <- subset(
  SNPs, chr == 'chrI_SK1' & position >= 53000 & position <= 57000)

# Selected region
ggplot(df_chrI_region, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual(values=c(spiked_color, ctrl_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
  annotate('text', x=53.5, y=6, label=unique(df_chrI_region$strain),
           size=3, colour=c(spiked_color, ctrl_color)) +
  scale_y_continuous(limits=c(-0.7, 8)) +
  annotate("text", x=SNPs_chrI$position/1000, y=-0.5,
           label = "|", colour='red') +
  annotate("text", x=57, y=-0.5, label = "SNPs", colour='red', size=2.5) +
  scale_fill_manual(values=c(spiked_color, ctrl_color)) +
  labs(title='', x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none", strip.background = element_blank()) +
  labs(y=paste0('Red1 occupancy')) +
  # ylim(0, 8) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# In silico spike-in versus control: compare peak number and overlap           #
# EDIT: Drop in silico experiment; just use spiked versus non-spiked low-Red1  #
#------------------------------------------------------------------------------#
spiked_narrow <- import_MACS2_peaks(
  WT_spike_rep_peaks[['narrow']], type='narrow')
nonspiked_narrow <- import_MACS2_peaks(
  nonspiked_rep_peaks[['narrow']], type='narrow')
spiked_broad <- import_MACS2_peaks(WT_spike_rep_peaks[['broad']], type='broad')
nonspiked_broad <- import_MACS2_peaks(
  nonspiked_rep_peaks[['broad']], type='broad')

# Discard S288C (spike-in) chromosomes and standardize names
spiked_narrow <- keepSeqlevels(spiked_narrow,
                               paste0('chr', as.roman(1:16), '_SK1'),
                               pruning.mode='coarse')
spiked_broad <- keepSeqlevels(spiked_broad,
                              paste0('chr', as.roman(1:16), '_SK1'),
                              pruning.mode='coarse')

renaming_map_names <- paste0('chr', as.roman(1:16), '_SK1')
renaming_map <- paste0('chr', as.roman(1:16))
names(renaming_map) <- renaming_map_names
spiked_narrow <- renameSeqlevels(spiked_narrow, renaming_map)
spiked_broad <- renameSeqlevels(spiked_broad, renaming_map)

length(spiked_narrow)
length(nonspiked_narrow)
length(spiked_broad)
length(nonspiked_broad)

# Mind the difference between:
#    1. all overlaps (equal to the reciprocal)
#    2. all peaks with overlaps (NOT equal to the reciprocal)
query <- spiked_narrow
subject <- nonspiked_narrow
message('Overlaps: ')
message('    A in B = B in A: ',
        sum(countOverlaps(query, subject)), ' = ',
        sum(countOverlaps(subject, query)))
message()
message('Peaks with overlaps: ')
message('    A in B != B in A: ',
        length(subsetByOverlaps(query, subject)), ' != ',
        length(subsetByOverlaps(subject, query)))

# Overlaps between pairs of samples
overlaps_between <- function(query, subject){
  message('Peaks in ', deparse(substitute(query)), ' vs ',
          deparse(substitute(subject)), ': ')
  message('   overlapping: ', length(subsetByOverlaps(query, subject)))
  message('   unique: ', length(query[!(query %over% subject)]))
  
  message('Peaks in ', deparse(substitute(subject)), ' vs ',
          deparse(substitute(query)), ': ')
  message('   overlapping: ', length(subsetByOverlaps(subject, query)))
  message('   unique: ', length(subject[!(subject %over% query)]))
  
  get_overlaps <- function(x, y) {
    overlaps <- list(overlap = subsetByOverlaps(x, y),
                     unique = x[!(x %over% y)])
    return(overlaps)
  }
  
  q_on_s <- get_overlaps(query, subject)
  s_on_q <- get_overlaps(subject, query)
  
  return(list(q_on_s=q_on_s, s_on_q=s_on_q))
}

on <- overlaps_between(spiked_narrow, nonspiked_narrow)
ob <- overlaps_between(spiked_broad, nonspiked_broad)

narrow_peak_counts <- as.data.frame(matrix(nrow=6, ncol=3))
colnames(narrow_peak_counts) <- c('Intersection', 'Count', 'Sample')

narrow_peak_counts$Intersection <- c('All', 'Intersected', 'Unique')
narrow_peak_counts$Sample <- c(rep('Non-spiked', 3),
                               rep('Spiked', 3))
narrow_peak_counts$Count <- c(length(on$s_on_q$overlap) + length(on$s_on_q$unique),
                              length(on$s_on_q$overlap),
                              length(on$s_on_q$unique),
                              length(on$q_on_s$overlap) + length(on$q_on_s$unique),
                              length(on$q_on_s$overlap),
                              length(on$q_on_s$unique))

broad_peak_counts <- as.data.frame(matrix(nrow=6, ncol=3))
colnames(broad_peak_counts) <- c('Intersection', 'Count', 'Sample')

broad_peak_counts$Intersection <- c('All', 'Intersected', 'Unique')
broad_peak_counts$Sample <- c(rep('Non-spiked', 3),
                              rep('Spiked', 3))
broad_peak_counts$Count <- c(length(ob$s_on_q$overlap) + length(ob$s_on_q$unique),
                             length(ob$s_on_q$overlap),
                             length(ob$s_on_q$unique),
                             length(ob$q_on_s$overlap) + length(ob$q_on_s$unique),
                             length(ob$q_on_s$overlap),
                             length(ob$q_on_s$unique))

# Order intersection
narrow_peak_counts$Intersection = factor(
  narrow_peak_counts$Intersection, levels=c('Unique', 'Intersected'))
broad_peak_counts$Intersection = factor(
  broad_peak_counts$Intersection, levels=c('Unique', 'Intersected'))

ggplot(narrow_peak_counts[
  narrow_peak_counts$Intersection %in% c('Unique', 'Intersected'), ],
  aes(Sample, Count, alpha=Intersection)) +
  geom_bar(stat='identity', colour='white', fill='Black') +
  scale_alpha_manual('', values=c(0.5, 1)) +
  ylim(0, 2900) +
  # scale_fill_manual('', values=c(insilico_spike_color, ctrl_color)) +
  labs(x='', y='Narrow peak count') +
  theme(axis.text.x = element_text(vjust=0.5))


ggplot(broad_peak_counts[
  broad_peak_counts$Intersection %in% c('Unique', 'Intersected'), ],
  aes(Sample, Count, alpha=Intersection)) +
  geom_bar(stat='identity', colour='white', fill='Black') +
  # ylim(0, 2600) +
  scale_alpha_manual('', values=c(0.5, 1)) +
  ylim(0, 2900) +
  # scale_fill_manual('', values=c(insilico_spike_color, ctrl_color)) +
  labs(x='', y='Broad peak count') +
  theme(axis.text.x = element_text(vjust=0.5))
