#------------------------------------------------------------------------------#
#                                                                              #
#                                   Figure 3                                   #
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
hop1_color <- wes_colors[5]
rec8_color <- wes_colors[3]


#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Spike-in normalization factor                                                #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors

# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah8151', ignore_case=TRUE))]

spinfs <- data.frame(
    sample=c('REC8', rep('rec8∆', 2)),
    spinf=c(1,
            suppressMessages(
                spikein_normalization_factor_from_counts(
                    ref_chip_counts=list(
                        read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                        read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                        read_counts_per_chr[['ah119spike-chip_trimmed']]),
                    ref_input_counts=list(
                        read_counts_per_chr[['ah119spikea-062817_trimmed']],
                        read_counts_per_chr[['ah119spikea-100917_trimmed']],
                        read_counts_per_chr[['ah119spike-inp_trimmed']]),
                    test_chip_counts=
                        read_counts_per_chr[['ah8151spike-chip_trimmed']],
                    test_input_counts=
                        read_counts_per_chr[['ah8151spike-inp_trimmed']])),
            suppressMessages(
                spikein_normalization_factor_from_counts(
                    ref_chip_counts=list(
                        read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                        read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                        read_counts_per_chr[['ah119spike-chip_trimmed']]),
                    ref_input_counts=list(
                        read_counts_per_chr[['ah119spikea-062817_trimmed']],
                        read_counts_per_chr[['ah119spikea-100917_trimmed']],
                        read_counts_per_chr[['ah119spike-inp_trimmed']]),
                    test_chip_counts=
                        read_counts_per_chr[['ah8151spike-red1-chip_trimmed']],
                    test_input_counts=
                        read_counts_per_chr[['ah8151spike-red1-inp_trimmed']]))))

# Order samples
spinfs$sample <- factor(spinfs$sample, levels=c('REC8', 'rec8∆'))

# Plot
message('Red1 amount in rec8∆ relative to wild type:')
message(round(mean(subset(spinfs, sample == 'rec8∆')$spinf), 3),
        '+/-',
        round(sd(subset(spinfs, sample == 'rec8∆')$spinf), 3))

mean_values <- spinfs %>% group_by(sample) %>% summarise(spinf_mean=mean(spinf))

ggplot(spinfs, aes(sample, spinf * 100, colour=sample, fill=sample)) +
    geom_hline(aes(yintercept = 100), linetype = 3) +
    stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
    geom_point(size=1.5, alpha=1) +
    scale_colour_manual('', values=c(wt_color, rec8_color), guide=FALSE) +
    scale_fill_manual('', values=c(wt_color, rec8_color), guide=FALSE) +
    scale_x_discrete(labels=c(expression(italic('REC8')),
                     expression(paste(italic('rec8'), Delta)))) +
    ylim(0, 100) +
    labs(title = '', x = '', y = 'Red1 amount\n(% of wild type)')


#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# Red1 fragment pileup on selected example chromosomes                         #
#------------------------------------------------------------------------------#
WT <- import_bedGraph(WT_1st_pileup)
rec8 <- import_bedGraph(rec8_1st_pileup)

WT_keep <- WT
rec8_keep <- rec8

### Normalize signal
# Compute genome-wide signal average
WT_avrg <- average_chr_signal(
    subset(WT, str_detect(WT@seqnames, '_SK1')))[[2]]
rec8_avrg <- average_chr_signal(
    subset(rec8, str_detect(rec8@seqnames, '_SK1')))[[2]]

# Compute spike-in normalization factor
subset(names(read_counts_per_chr),
       str_detect(names(read_counts_per_chr), '8115'))

nf <- spikein_normalization_factor_from_counts(
    ref_chip_counts=list(read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                         read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                         read_counts_per_chr[['ah119spike-chip_trimmed']]),
    ref_input_counts=list(read_counts_per_chr[['ah119spikea-062817_trimmed']],
                          read_counts_per_chr[['ah119spikea-100917_trimmed']],
                          read_counts_per_chr[['ah119spike-inp_trimmed']]),
    test_chip_counts=read_counts_per_chr[['ah8151spike-red1-chip_trimmed']],
    test_input_counts=read_counts_per_chr[['ah8151spike-red1-inp_trimmed']])

# Normalize by average signal and then spike-in factor
WT[str_detect(WT@seqnames, '_SK1')]$score <-
    WT[str_detect(WT@seqnames, '_SK1')]$score / WT_avrg

rec8[str_detect(rec8@seqnames, '_SK1')]$score <-
    rec8[str_detect(rec8@seqnames, '_SK1')]$score / rec8_avrg * nf


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
    WT, sample_name='REC8', normalize_by_genome_mean=FALSE,
    window_size=100, chr='chrI_SK1')
rec8_chrI <- make_df_for_plotting(
    rec8, sample_name='rec8∆', normalize_by_genome_mean=FALSE,
    window_size=100, chr='chrI_SK1')

df_chrI <- rbind(WT_chrI, rec8_chrI)

WT_chrII <- make_df_for_plotting(
    WT, sample_name='REC8', normalize_by_genome_mean=FALSE,
    window_size=100, chr='chrII_SK1')
rec8_chrII <- make_df_for_plotting(
    rec8, sample_name='rec8∆', normalize_by_genome_mean=FALSE,
    window_size=200, chr='chrII_SK1')

df_chrII <- rbind(WT_chrII, rec8_chrII)

# Order strains
df_chrI$strain = factor(df_chrI$strain, levels=c('REC8', 'rec8∆'))
df_chrII$strain = factor(df_chrII$strain, levels=c('REC8', 'rec8∆'))

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

ggplot(df_chrI, aes(position / 1000, signal, colour=strain)) +
    scale_color_manual(values=c(wt_color, rec8_color)) +
    geom_area(position='identity', aes(fill=strain, colour=strain),
              alpha=0.75, size=0.25) +
    geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
    annotate("text", x=cen_mid$cen_mid / 1000, y=-1, label="CEN",
             size=2, colour='black') +
    scale_fill_manual(values=c(wt_color, rec8_color)) +
    labs(title='',
         x=paste0('Position on chr I (Kb)'), y='') +
    theme(legend.position = "none") +
    labs(y=paste0('Red1 occupancy')) +
    ylim(-1, 12) +
    theme(
        strip.background=element_blank(),
        strip.text=element_blank()
    )

cen_mid <- cen_midpoint(chr='chrII_SK1', ref_genome='SK1_S288C', y_coord=0)

ggplot(df_chrII, aes(position / 1000, signal, colour=strain)) +
    scale_color_manual('', values=c(wt_color, rec8_color)) +
    geom_area(position='identity', aes(fill=strain, colour=strain),
              alpha=0.75, size=0.25) +
    geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
    annotate("text", x=cen_mid$cen_mid / 1000, y=-1, label="CEN",
             size=2, colour='black') +
    scale_fill_manual('', values=c(wt_color, rec8_color)) +
    annotate('text', x=760, y=11, label="italic('REC8')",
             parse=TRUE, size=3, colour=wt_color) +
    annotate('text', x=760, y=9, label="paste(italic('rec8'), Delta)",
             parse=TRUE, size=3, colour=rec8_color) +
    labs(title='',
         x=paste0('Position on chr II (Kb)'), y='') +
    theme(legend.position = "none") +
    labs(y=paste0('Red1 occupancy')) +
    ylim(-1, 12) +
    theme(
        strip.background=element_blank(),
        strip.text=element_blank()
    )


#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Spike-in normalization factor                                                #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors

# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah9120', ignore_case=TRUE))]

spinfs <- data.frame(
  sample=c('HOP1', rep('hop1∆', 2)),
  spinf=c(1,
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=list(
                read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                read_counts_per_chr[['ah119spike-chip_trimmed']]),
              ref_input_counts=list(
                read_counts_per_chr[['ah119spikea-062817_trimmed']],
                read_counts_per_chr[['ah119spikea-100917_trimmed']],
                read_counts_per_chr[['ah119spike-inp_trimmed']]),
              test_chip_counts=
                read_counts_per_chr[['ah9120spikeb-100917_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9120spikea-100917_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=list(
                read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                read_counts_per_chr[['ah119spike-chip_trimmed']]),
              ref_input_counts=list(
                read_counts_per_chr[['ah119spikea-062817_trimmed']],
                read_counts_per_chr[['ah119spikea-100917_trimmed']],
                read_counts_per_chr[['ah119spike-inp_trimmed']]),
              test_chip_counts=
                read_counts_per_chr[['ah9120spike-chip_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9120spike-inp_trimmed']]))))

# Order samples
spinfs$sample <- factor(spinfs$sample, levels=c('HOP1', 'hop1∆'))

# Plot
message('Red1 amount in hop1∆ relative to wild type:')
message(round(mean(subset(spinfs, sample == 'hop1∆')$spinf), 3),
        '+/-',
        round(sd(subset(spinfs, sample == 'hop1∆')$spinf), 3))

mean_values <- spinfs %>% group_by(sample) %>% summarise(spinf_mean=mean(spinf))

ggplot(spinfs, aes(sample, spinf * 100, colour=sample, fill=sample)) +
  geom_hline(aes(yintercept = 100), linetype = 3) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_point(size=1.5, alpha=1) +
  # geom_jitter(width = 0.1) +
  scale_colour_manual('', values=c(wt_color, hop1_color), guide=FALSE) +
  scale_fill_manual('', values=c(wt_color, hop1_color), guide=FALSE) +
  scale_x_discrete(labels=c(expression(italic('HOP1')),
                            expression(paste(italic('hop1'), Delta)))) +
  ylim(0, 100) +
  labs(title = '', x = '', y = 'Red1 amount\n(% of wild type)')


#------------------------------------------------------------------------------#
#                                   Panel D                                    #
# Red1 fragment pileup on selected example chromosomes                         #
#------------------------------------------------------------------------------#
WT <- import_bedGraph(WT_1st_pileup)
hop1 <- import_bedGraph(hop1_1st_pileup)

### Normalize signal
# Compute genome-wide signal average
WT_avrg <- average_chr_signal(
  subset(WT, str_detect(WT@seqnames, '_SK1')))[[2]]
hop1_avrg <- average_chr_signal(
  subset(hop1, str_detect(hop1@seqnames, '_SK1')))[[2]]

# Compute spike-in normalization factor
subset(names(read_counts_per_chr),
       str_detect(names(read_counts_per_chr), '9120'))

nf <- spikein_normalization_factor_from_counts(
  ref_chip_counts=list(read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                       read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                       read_counts_per_chr[['ah119spike-chip_trimmed']]),
  ref_input_counts=list(read_counts_per_chr[['ah119spikea-062817_trimmed']],
                        read_counts_per_chr[['ah119spikea-100917_trimmed']],
                        read_counts_per_chr[['ah119spike-inp_trimmed']]),
  test_chip_counts=read_counts_per_chr[['ah9120spikeb-100917_trimmed']],
  test_input_counts=read_counts_per_chr[['ah9120spikea-100917_trimmed']])

# Normalize by average signal and then spike-in factor
WT[str_detect(WT@seqnames, '_SK1')]$score <-
  WT[str_detect(WT@seqnames, '_SK1')]$score / WT_avrg

hop1[str_detect(hop1@seqnames, '_SK1')]$score <-
  hop1[str_detect(hop1@seqnames, '_SK1')]$score / hop1_avrg * nf


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
  WT, sample_name='HOP1', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')
hop1_chrI <- make_df_for_plotting(
  hop1, sample_name='hop1∆', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrI_SK1')

df_chrI <- rbind(WT_chrI, hop1_chrI)

WT_chrII <- make_df_for_plotting(
  WT, sample_name='HOP1', normalize_by_genome_mean=FALSE,
  window_size=100, chr='chrII_SK1')
hop1_chrII <- make_df_for_plotting(
  hop1, sample_name='hop1∆', normalize_by_genome_mean=FALSE,
  window_size=200, chr='chrII_SK1')

df_chrII <- rbind(WT_chrII, hop1_chrII)

# Order strains
df_chrI$strain = factor(df_chrI$strain, levels=c('HOP1', 'hop1∆'))
df_chrII$strain = factor(df_chrII$strain, levels=c('HOP1', 'hop1∆'))

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

ggplot(df_chrI, aes(position / 1000, signal, colour=strain)) +
  # geom_hline(yintercept=1, lty=3) +
  scale_color_manual(values=c(wt_color, hop1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.75, size=0.25) +
  geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
  annotate("text", x=cen_mid$cen_mid / 1000, y=-1, label="CEN",
           size=2, colour='black') +
  # facet_grid(strain ~ .) +
  scale_fill_manual(values=c(wt_color, hop1_color)) +
  labs(title='',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  ylim(-1, 16) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

cen_mid <- cen_midpoint(chr='chrII_SK1', ref_genome='SK1_S288C', y_coord=0)

ggplot(df_chrII, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual('', values=c(wt_color, hop1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.75, size=0.25) +
  geom_point(data=cen_mid, aes(cen_mid/1000, 0.2), size=1, colour='black') +
  annotate("text", x=cen_mid$cen_mid / 1000, y=-1, label="CEN",
           size=2, colour='black') +
  scale_fill_manual('', values=c(wt_color, hop1_color)) +
  annotate('text', x=760, y=14, label="italic('HOP1')", parse=TRUE,
           size=3, colour=wt_color) +
  annotate('text', x=760, y=11, label="paste(italic('hop1'), Delta)", parse=TRUE,
           size=3, colour=hop1_color) +
  labs(title='',
       x=paste0('Position on chr II (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  ylim(-1, 16) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

#------------------------------------------------------------------------------#
#                                   Panel E                                    #
# Average Red1 signal per chromosome                                           #
#------------------------------------------------------------------------------#
WT <- import_bedGraph(WT_1st_pileup)
hop1 <- import_bedGraph(hop1_1st_pileup)

# Drop spike-in genome
WT_SK1 <- keepSeqlevels(
    WT, paste0('chr', as.roman(1:16), '_SK1'), pruning.mode="coarse")
hop1_SK1 <- keepSeqlevels(
    hop1, paste0('chr', as.roman(1:16), '_SK1'), pruning.mode="coarse")

WT_avrg <- average_chr_signal(WT_SK1)[[2]]
hop1_avrg <- average_chr_signal(hop1_SK1)[[2]]

# Compute spike-in normalization factor
subset(names(read_counts_per_chr),
       str_detect(names(read_counts_per_chr), '9120'))

nf <- spikein_normalization_factor_from_counts(
    ref_chip_counts=list(read_counts_per_chr[['ah119spikeb-062817_trimmed']],
                         read_counts_per_chr[['ah119spikeb-100917_trimmed']],
                         read_counts_per_chr[['ah119spike-chip_trimmed']]),
    ref_input_counts=list(read_counts_per_chr[['ah119spikea-062817_trimmed']],
                          read_counts_per_chr[['ah119spikea-100917_trimmed']],
                          read_counts_per_chr[['ah119spike-inp_trimmed']]),
    test_chip_counts=read_counts_per_chr[['ah9120spikeb-100917_trimmed']],
    test_input_counts=read_counts_per_chr[['ah9120spikea-100917_trimmed']])


# Normalize by average signal and then spike-in factor
WT_SK1$score <- WT_SK1$score / WT_avrg
hop1_SK1$score <- hop1_SK1$score / hop1_avrg * nf

# Compute average signal per chromosome using normalized signal
WT_chrs <- average_chr_signal(WT_SK1)[[1]]
WT_chrs$strain <- 'HOP1'
hop1_chrs <- average_chr_signal(hop1_SK1)[[1]]
hop1_chrs$strain <- 'hop1∆'

chr_signal <- rbind(WT_chrs, hop1_chrs)

# Add chromosome sizes
chr_lengths <- seqlengths(get_chr_coordinates())
chr_lengths <- data.frame(chr=names(chr_lengths), chr_len=chr_lengths)
chr_signal <- dplyr::full_join(
    chr_signal, subset(chr_lengths, str_detect(chr, '_SK1')), by=c('chr'='chr'))

# Order samples
chr_signal$strain <- factor(chr_signal$strain, levels=c('HOP1', 'hop1∆'))

facet_names <- list('HOP1'=expression(italic('HOP1')),
                    'hop1∆'=expression(paste(italic('hop1'), Delta)))
facet_labeller <- function(variable, value){return(facet_names[value])}

ggplot(chr_signal, aes(chr_len / 10^6, avrg_signal, colour=strain)) +
    geom_point(stat="identity", size=2, alpha=0.8) +
    facet_wrap(~strain, ncol=3, labeller=facet_labeller) +
    scale_color_manual('', values=c(wt_color, hop1_color), guide=FALSE) +
    labs(title = '', x = 'Chromosome size (Mb)',
         y='Red1 occupancy\n(average per bp)') +
    geom_hline(yintercept=1, lty=3) +
    scale_x_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1.0, 1.5)) +
    ylim(0, 1.7) +
    theme(strip.background = element_blank())