#------------------------------------------------------------------------------#
#                                                                              #
#                                   Figure 4                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(GenomicRanges)
library(tidyverse)
library(RColorBrewer)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load code to compute spike-in normalization factor
source(here('helper_spikein_normalization_factor.R'))
# Load IO code
source(here('helper_io.R'))
# Load analysis functions
source(here('helper_analysis.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
strain_colors <- c('black', rev(brewer.pal(9, "YlOrRd")[4:8]))

#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Spike-in normalization factors                                               #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors
# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('gh2ax', ignore_case=TRUE))]


WT_chip_rep_counts <- read_counts_per_chr[['ah119spike-gh2ax-chip_trimmed']]
WT_input_rep_counts <- read_counts_per_chr[['ah119spike-gh2ax-inp_trimmed']]

nfs <- data.frame(
  sample=c('RED1/RED1',
           'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
           'red1_ycs4/red1∆', 'spo11-YF/spo11-YF'),
  nf=c(1,
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=
                read_counts_per_chr[['ah8218spike-gh2ax-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah8218spike-gh2ax-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=
                read_counts_per_chr[['ah8220spike-gh2ax-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah8220spike-gh2ax-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=
                read_counts_per_chr[['ah7011spike-gh2ax-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah7011spike-gh2ax-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=
                read_counts_per_chr[['ah8219spike-gh2ax-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah8219spike-gh2ax-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=
                read_counts_per_chr[['ah4206spike-gh2ax-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah4206spike-gh2ax-inp_trimmed']]))))

# Order samples
nfs$sample <- factor(nfs$sample, levels=c(
    'RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
    'red1_ycs4/red1∆', 'spo11-YF/spo11-YF'))

# Plot
ggplot(nfs, aes(sample, nf * 100, fill=sample)) +
  geom_hline(aes(yintercept=100), linetype=3) +
  geom_bar(stat="identity", width=0.5, alpha=0.5) +
  scale_colour_manual('', values=strain_colors, guide=FALSE) +
  scale_fill_manual('', values=strain_colors, guide=FALSE) +
  scale_x_discrete(
    labels=c(
        expression(
            italic('RED1') * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('RED1'),
            italic('red1') * Delta * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('red1'['ycs4S']),
            italic('red1'['ycs4S']) * '/' * italic('red1') * Delta,
            italic('spo11-YF/spo11-YF')))) +
  ylim(0, 105) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title = '',
       x='\n\nRed1 dosage strain series',
       y=bquote(atop(paste(gamma, '-H2AX amount'), '(% of wild type)')))


#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# Signal onm example chromosome                                                #
#------------------------------------------------------------------------------#
fragment_pileups <- lapply(gH2AX_pileup, import_bedGraph)

make_df_for_plotting <- function(gr, window_size=1, genome='SK1_S288C',
                                 chr='chrI_SK1', region=c(194000, 208000)) {
  t0 <- proc.time()[3]
  
  message('Compute binned score...')
  gr <- bedGraph_to_binned_score(
    gr, mean_normalize=FALSE, tile_width=window_size,
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
  
  message('~~~')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  gr
}

chr <- lapply(fragment_pileups, make_df_for_plotting,
              window_size=200, chr='chrXI_SK1')

# Add sample names to df
for (i in seq_along(chr)) {
  message('>>> ', names(chr)[i])
  chr[[i]]$strain <- names(chr)[i]
}

### Normalize signal
# (also save non-normalized signal)
# Compute genome-wide signal average
compute_avrg <- function(x) {
  # Drop S288C values
  x <- keepSeqlevels(x, paste0('chr', as.roman(1:16), '_SK1'),
                     pruning.mode="coarse")
  sum(width(x) * score(x)) / sum(width(x))
}

fragment_pileup_avrgs <- lapply(fragment_pileups, compute_avrg)

# Rename samples in NF table to match new data
nfs_renamed <- mutate(
    nfs,
    sample=recode(
        sample,
        'RED1/RED1'='WT',
        'red1_ycs4/RED1'='ycs4S/RED1',
        'red1∆/RED1'='red1/RED1',
        'red1_ycs4/red1_ycs4'='ycs4S/ycs4S',
        'red1_ycs4/red1∆'='ycs4S/red1'))

# Normalize by average signal and then spike-in factor
for (i in seq_along(chr)) {
  message('>>> ', names(chr)[i])
  chr[[i]]$non_norm_signal <- chr[[i]]$signal
  chr[[i]]$signal <- chr[[i]]$signal / fragment_pileup_avrgs[[names(chr)[i]]]
  nf <- as.numeric(subset(nfs_renamed, sample == names(chr)[i], select='nf'))
  message('    nf: ', nf)
  chr[[i]]$signal <- chr[[i]]$signal * nf
}

df_chr <- do.call('rbind', chr)

# Rename strains
df_chr <- mutate(
    df_chr,
    strain=recode(
        strain,
        'WT'='RED1/RED1',
        'ycs4S/RED1'='red1_ycs4/RED1',
        'red1/RED1'='red1∆/RED1',
        'ycs4S/ycs4S'='red1_ycs4/red1_ycs4',
        'ycs4S/red1'='red1_ycs4/red1∆'))

# Order strains
df_chr$strain <- factor(df_chr$strain, levels=c(
    'RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
    'red1_ycs4/red1∆', 'spo11-YF/spo11-YF'))

# Plot signal
ggplot(df_chr, aes(position / 1000, non_norm_signal, colour=strain)) +
  scale_color_manual(values=strain_colors) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
  annotate("text", x=560, y=3.2, label=c(
      expression(
          italic('RED1') * '/' * italic('RED1'),
          italic('red1'['ycs4S']) * '/' * italic('RED1'),
          italic('red1') * Delta * '/' * italic('RED1'),
          italic('red1'['ycs4S']) * '/' * italic('red1'['ycs4S']),
          italic('red1'['ycs4S']) * '/' * italic('red1') * Delta,
          italic('spo11-YF/spo11-YF'))),
      size=3.1, colour=strain_colors) +
  scale_fill_manual(values=strain_colors) +
  labs(title='',
       x=paste0('Position on chr XI (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('γ-H2AX occupancy')) +
  ylim(0, 3.6) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

# Plot normalized signal
ggplot(df_chr, aes(position / 1000, signal, colour=strain)) +
    scale_color_manual('', values=strain_colors) +
    scale_fill_manual('', values=strain_colors) +
    geom_area(position='identity', aes(fill=strain), alpha=0.5, size=0.25) +
    facet_grid(strain ~ .) +
    annotate("text", x=600, y=3.9, label=c(
        expression(
            italic('RED1') * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('RED1'),
            italic('red1') * Delta * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('red1'['ycs4S']),
            italic('red1'['ycs4S']) * '/' * italic('red1') * Delta,
            italic('spo11-YF/spo11-YF'))),
        size=3.1, colour=strain_colors) +
    labs(title='',
         x=paste0('Position on chr XI (Kb)'),
         y=expression(paste(gamma, '-H2AX occupancy'))) +
    scale_y_continuous(limits=c(0, 4.9), breaks=c(0, 4)) +
    
    theme(
        legend.position = "none",
        strip.background=element_blank(),
        strip.text = element_blank()
    )

#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Signal from telomeres                                                        #
#------------------------------------------------------------------------------#
chr <- lapply(fragment_pileups, make_df_for_plotting,
              window_size=20, chr='chrXI_SK1')

# Add sample names to df
for (i in seq_along(chr)) {
    message('>>> ', names(chr)[i])
    chr[[i]]$strain <- names(chr)[i]
}

### Normalize signal
# (also save non-normalized signal)
# Compute genome-wide signal average
compute_avrg <- function(x) {
    # Drop S288C values
    x <- keepSeqlevels(x, paste0('chr', as.roman(1:16), '_SK1'),
                       pruning.mode="coarse")
    sum(width(x) * score(x)) / sum(width(x))
}

fragment_pileup_avrgs <- lapply(fragment_pileups, compute_avrg)

# Rename samples in NF table to match new data
nfs_renamed <- mutate(
    nfs,
    sample=recode(
        sample,
        'RED1/RED1'='WT',
        'red1_ycs4/RED1'='ycs4S/RED1',
        'red1∆/RED1'='red1/RED1',
        'red1_ycs4/red1_ycs4'='ycs4S/ycs4S',
        'red1_ycs4/red1∆'='ycs4S/red1'))

# Normalize by average signal and then spike-in factor
for (i in seq_along(chr)) {
    message('>>> ', names(chr)[i])
    chr[[i]]$non_norm_signal <- chr[[i]]$signal
    chr[[i]]$signal <- chr[[i]]$signal / fragment_pileup_avrgs[[names(chr)[i]]]
    nf <- as.numeric(subset(nfs_renamed, sample == names(chr)[i], select='nf'))
    message('    nf: ', nf)
    chr[[i]]$signal <- chr[[i]]$signal * nf
}

df_chr <- do.call('rbind', chr)

# Rename strains
df_chr <- mutate(
    df_chr,
    strain=recode(
        strain,
        'WT'='RED1/RED1',
        'ycs4S/RED1'='red1_ycs4/RED1',
        'red1/RED1'='red1∆/RED1',
        'ycs4S/ycs4S'='red1_ycs4/red1_ycs4',
        'ycs4S/red1'='red1_ycs4/red1∆'))

# Order strains
df_chr$strain <- factor(df_chr$strain, levels=c(
    'RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
    'red1_ycs4/red1∆', 'spo11-YF/spo11-YF'))

ggplot(subset(df_chr, position < 10000), aes(position / 1000, signal, colour=strain)) +
    scale_color_manual('', values=strain_colors) +
    scale_fill_manual('', values=strain_colors) +
    geom_area(position='identity', aes(fill=strain), alpha=0.5, size=0.25) +
    facet_grid(strain ~ .) +
    annotate("text", x=6, y=3.9, label=c(
        expression(
            italic('RED1') * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('RED1'),
            italic('red1') * Delta * '/' * italic('RED1'),
            italic('red1'['ycs4S']) * '/' * italic('red1'['ycs4S']),
            italic('red1'['ycs4S']) * '/' * italic('red1') * Delta,
            italic('spo11-YF/spo11-YF'))),
        size=3.1, colour=strain_colors) +
    labs(title='',
         x=paste0('Position on chr XI (Kb)'),
         y=expression(paste(gamma, '-H2AX occupancy'))) +
    scale_y_continuous(limits=c(0, 4.9), breaks=c(0, 4)) +
    
    theme(
        legend.position = "none",
        strip.background=element_blank(),
        strip.text = element_blank()
    )