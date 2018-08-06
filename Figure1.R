#------------------------------------------------------------------------------#
#                                                                              #
#                                   Figure 1                                   #
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
# Load IO code
source(here('helper_io.R'))
# Load analysis functions
source(here('helper_analysis.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
wes_colors <- wes_palette(name = "Darjeeling1", n=8, type = "continuous")
wt_color <- 'black'
low_Red1_color <- wes_colors[7]

#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# red1_ycs4 mutant                                                             #
# EDIT: Use red1_ycs4S instead (because non-spiked sample is also available)   #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors

# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah7011', ignore_case=TRUE))]

Red1_levels <- data.frame(
  sample=c('RED1', rep('red1_ycs4', 2)),
  Red1=c(1,
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
                # read_counts_per_chr[['ah9048spikeb-062817_trimmed']],
                read_counts_per_chr[['ah7011spike-chip_trimmed']],
              test_input_counts=
                # read_counts_per_chr[['ah9048spikea-062817_trimmed']])),
                read_counts_per_chr[['ah7011spike-inp_trimmed']])),
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
                # read_counts_per_chr[['ah9048spikeb-100917_trimmed']],
                read_counts_per_chr[['ah7011spike-red1-chip_trimmed']],
              test_input_counts=
                # read_counts_per_chr[['ah9048spikea-100917_trimmed']]))))
                read_counts_per_chr[['ah7011spike-red1-inp_trimmed']]))),
  se=0)

# Add published Red1 levels (Markowitz et al., PLoS Genet. 2017):
Published_Red1_levels <- data.frame(
  sample=c('RED1', 'red1_ycs4'),
  Red1=c(1, 0.2462982),
  se=c(0, 8), method='Western blot')

Red1_levels$method <- 'SNP-ChIP'

Red1_levels <- rbind(Red1_levels, Published_Red1_levels)

# Order samples
Red1_levels$sample <- factor(
  Red1_levels$sample,
  levels=c('RED1', 'red1_ycs4'))
Red1_levels$method <- factor(
  Red1_levels$method, levels=c('SNP-ChIP', 'Western blot'))

# Plot
ggplot(Red1_levels, aes(sample, Red1 * 100, colour=method, fill=method,
                        ymin=Red1 * 100 - se, ymax=Red1 * 100 + se)) +
    scale_colour_manual('', values=c(low_Red1_color, 'transparent')) +
    scale_fill_manual('', values=c(low_Red1_color, 'black')) +
    geom_hline(aes(yintercept=100), linetype=3) +
    geom_point(size=1.5, position=position_dodge(0.75)) +
    stat_summary(fun.y=mean, geom="bar", position="dodge", width = 0.75,
                 alpha=0.25, colour = NA) +
    geom_errorbar(
        width=.2, position=position_dodge(width=0.75),
        colour=c("transparent", "transparent", "black",
                 "transparent", "transparent")) +
    labs(title='', x='', y='Red1 amount\n(% of wild type)') +
    scale_x_discrete(labels=c(
      expression(italic('RED1'), italic('red1'[italic('ycs4S')])))) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    guides(colour = guide_legend(override.aes = list(shape = NA)))


#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Red1 fragment pileup example before and after spike-in normalization         #
#------------------------------------------------------------------------------#
# Use non-spiked data
WT <- import_bedGraph(nonspiked_rep_pileup[['WT']])
low_Red1 <- import_bedGraph(nonspiked_rep_pileup[['red1_ycs4S']])

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
  WT, sample_name='RED1', normalize_by_genome_mean=FALSE, window_size=100,
  genome='SK1', chr='chrI')
low_Red1_chrI <- make_df_for_plotting(
  low_Red1, sample_name='red1_ycs4', normalize_by_genome_mean=FALSE,
  window_size=100, genome='SK1', chr='chrI')

WT_chrI_region <- make_df_for_plotting(
  WT, sample_name='RED1', normalize_by_genome_mean=FALSE, window_size=1,
  genome='SK1', chr='chrI', region=c(194000, 208000))
low_Red1_chrI_region <- make_df_for_plotting(
  low_Red1, sample_name='red1_ycs4', normalize_by_genome_mean=FALSE,
  window_size=1, genome='SK1', chr='chrI', region=c(194000, 208000))

df_chrI <- rbind(WT_chrI, low_Red1_chrI)
df_chrI_region <- rbind(WT_chrI_region, low_Red1_chrI_region)

### Normalize signal
# Compute genome-wide signal average
compute_avrg <- function(x) {
  (sum(width(x) * score(x)) / sum(width(x)))
}

WT_avrg <- compute_avrg(WT)
low_Red1_avrg <- compute_avrg(low_Red1)

# Compute spike-in normalization factor
sinf <- colMeans(subset(Red1_levels, sample == 'red1_ycs4', select='Red1'))

# Normalize by average signal and then spike-in factor
# Whole chromosome
df_chrI$normalized_signal <- df_chrI$signal
df_chrI[df_chrI$strain == 'red1_ycs4', 'normalized_signal'] <- 
  subset(df_chrI, strain == 'red1_ycs4')$signal / 
  low_Red1_avrg * WT_avrg * sinf

# Selected region
df_chrI_region$normalized_signal <- df_chrI_region$signal
df_chrI_region[df_chrI_region$strain == 'RED1', 'normalized_signal'] <-
  subset(df_chrI_region, strain == 'RED1', signal)
df_chrI_region[df_chrI_region$strain == 'red1_ycs4', 'normalized_signal'] <- 
  subset(df_chrI_region, strain == 'red1_ycs4')$signal / 
  low_Red1_avrg * WT_avrg * sinf

# Order strains
df_chrI$strain = factor(df_chrI$strain, levels=c('RED1', 'red1_ycs4'))
df_chrI_region$strain = factor(df_chrI_region$strain,
                               levels=c('RED1', 'red1_ycs4'))

### Plot before normalization
# Whole chromosome
ggplot(df_chrI, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual(values=c(wt_color, low_Red1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
    annotate('text', x=20, y=6, label=expression(
      italic('RED1'), italic('red1'[italic('ycs4S')])),
      size=3, colour=c(wt_color, low_Red1_color)) +
  scale_fill_manual(values=c(wt_color, low_Red1_color)) +
  labs(title='Before spike-in normalization',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position="none", strip.background=element_blank()) +
  labs(y=paste0('Red1 occupancy')) +
  scale_y_continuous(limits=c(0, 9), breaks=c(0, 8)) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

# Selected region
ggplot(df_chrI_region, aes(position / 1000, signal, colour=strain)) +
  scale_color_manual(values=c(wt_color, low_Red1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
    annotate('text', x=196, y=6, label=expression(
        italic('RED1'), italic('red1'[italic('ycs4S')])),
           size=3, colour=c(wt_color, low_Red1_color)) +
  scale_fill_manual(values=c(wt_color, low_Red1_color)) +
  labs(title='Before spike-in normalization',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none", strip.background = element_blank()) +
  labs(y=paste0('Red1 occupancy')) +
  scale_y_continuous(limits=c(0, 9), breaks=c(0, 8)) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
    )

# Plot after normalization
# Whole chromosome
ggplot(df_chrI, aes(position / 1000, normalized_signal, colour=strain)) +
  scale_color_manual(values=c(wt_color, low_Red1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
  scale_fill_manual(values=c(wt_color, low_Red1_color)) +
  labs(title='After spike-in normalization',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  scale_y_continuous(limits=c(0, 9), breaks=c(0, 8)) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )

# Selected region
ggplot(df_chrI_region, aes(position / 1000, normalized_signal, colour=strain)) +
  scale_color_manual(values=c(wt_color, low_Red1_color)) +
  geom_area(position='identity', aes(fill=strain, colour=strain),
            alpha=0.5, size=0.25) +
  facet_grid(strain ~ .) +
  scale_fill_manual(values=c(wt_color, low_Red1_color)) +
  labs(title='After spike-in normalization',
       x=paste0('Position on chr I (Kb)'), y='') +
  theme(legend.position = "none") +
  labs(y=paste0('Red1 occupancy')) +
  scale_y_continuous(limits=c(0, 9), breaks=c(0, 8)) +
  theme(
    strip.background=element_blank(),
    strip.text=element_blank()
  )


#------------------------------------------------------------------------------#
#                                   Panel D                                    #
# Red1 dosage strain series                                                    #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah8220', ignore_case=TRUE))]

# Use all WT replicates
WT_chip_rep_counts <- list(
  read_counts_per_chr[['ah119spikeb-062817_trimmed']],
  read_counts_per_chr[['ah119spikeb-100917_trimmed']],
  read_counts_per_chr[['ah119spike-chip_trimmed']])

WT_input_rep_counts <- list(
  read_counts_per_chr[['ah119spikea-062817_trimmed']],
  read_counts_per_chr[['ah119spikea-100917_trimmed']],
  read_counts_per_chr[['ah119spike-inp_trimmed']])
  
Red1_levels <- data.frame(
  sample=c('RED1/RED1', rep('red1_ycs4/RED1', 2), rep('red1∆/RED1', 2),
           rep('red1_ycs4/red1_ycs4', 2), rep('red1_ycs4/red1∆', 2)),
  Red1=c(1,
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8218spike-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8218spike-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8218spike-red1-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8218spike-red1-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8220spike-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8220spike-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8220spike-red1-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8220spike-red1-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah7011spike-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah7011spike-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah7011spike-red1-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah7011spike-red1-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8219spike-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8219spike-inp_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=WT_chip_rep_counts,
              ref_input_counts=WT_input_rep_counts,
              test_chip_counts=read_counts_per_chr[[
                'ah8219spike-red1-chip_trimmed']],
              test_input_counts=read_counts_per_chr[[
                'ah8219spike-red1-inp_trimmed']]))),
  se=0)

# Order samples
Red1_levels$sample <- factor(Red1_levels$sample, levels=c(
  'RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
  'red1_ycs4/red1∆'))

# Plot
Published_Red1_levels <- data.frame(
  sample=c('RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4', 
           'red1_ycs4/red1∆'),
  Red1=c(1, 0.4359851, 0.4537080, 0.2462982, 0.1243278),
  se=c(0, 11, 10, 8, 4), method='Western blot')


Red1_levels$method <- 'SNP-ChIP'

Red1_levels <- rbind(Red1_levels, Published_Red1_levels)

# Order samples
Red1_levels$sample <- factor(
  Red1_levels$sample,
  levels=c('RED1/RED1', 'red1_ycs4/RED1', 'red1∆/RED1', 'red1_ycs4/red1_ycs4',
           'red1_ycs4/red1∆'))
Red1_levels$method <- factor(
  Red1_levels$method, levels=c('SNP-ChIP', 'Western blot'))

ggplot(Red1_levels, aes(sample, Red1 * 100, colour=method, fill=method,
                        ymin=Red1 * 100 - se, ymax=Red1 * 100 + se)) +
    scale_colour_manual('', values=c(low_Red1_color, 'transparent')) +
    scale_fill_manual('', values=c(low_Red1_color, 'black')) +
    geom_hline(aes(yintercept=100), linetype=3) +
    geom_point(size=1.5, position=position_dodge(0.75)) +
    stat_summary(fun.y=mean, geom="bar", position="dodge", width = 0.75,
                 alpha=0.25, colour = NA) +
    geom_errorbar(
        width=.2, position=position_dodge(width=0.75),
        colour=c("transparent", "transparent",
                 rep(c("black", "transparent", "transparent"), 4))) +
    labs(title='',
         x='\n\nRed1 dosage strain series', y='Red1 amount\n(% of wild type)') +
    scale_x_discrete(labels=c(
        expression(
            italic('RED1') * '/' * italic('RED1'),
            italic('red1'[italic('ycs4S')]) * '/' * italic('RED1'),
            italic('red1') * Delta * '/' * italic('RED1'),
            italic('red1'[italic('ycs4S')]) * '/' * italic('red1'[italic('ycs4S')]),
            italic('red1'[italic('ycs4S')]) * '/' * italic('red1') * Delta))) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    guides(colour = guide_legend(override.aes = list(shape = NA)))
