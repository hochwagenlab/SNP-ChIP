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
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Histogram of distances between consecutive SNPs                              #
#------------------------------------------------------------------------------#
SNPs <- read_tsv(SNPs)

# Keep only SK1 chromosomes and order table by chromosome and position
SNPs <- SNPs[str_detect(SNPs$chr, '_SK1'), ]
message(nrow(SNPs), ' annotated SNPs')
SNPs <- SNPs[order(SNPs$chr, SNPs$position), ]
SNPs$chr <- factor(SNPs$chr, levels=paste0('chr', as.roman(1:16), '_SK1'))
SNPs <- SNPs[order(SNPs$chr), ]

# Compute distances between consecutive SNPs
SNPs <- SNPs %>% group_by(chr) %>% mutate(Distance = c(NA, diff(position)))

mean_value <- mean(SNPs$Distance, na.rm=TRUE)
median_value <- median(SNPs$Distance, na.rm=TRUE)

p <- ggplot(SNPs, aes(Distance)) +
  geom_histogram(aes(y= ..count..), bins=50, colour='grey30',
                 alpha=0.25, na.rm = T) +
  xlim(-5, 500) +
  #geom_vline(xintercept = mean_value, color = 'red') +
  geom_segment(aes(x=mean_value, y=0, xend=mean_value, yend=3000),
               linetype='dotted', color='red', size=0.4) +
  #geom_vline(xintercept = median_value, color = 'blue') +
  geom_segment(aes(x=median_value, y=0, xend=median_value, yend=5000),
               linetype='dotted', color='red', size=0.4) +
  labs(title='', x='Distance (bp)') +
  annotate('text', x=mean_value, y=4000,
           # fontface='bold',
           label=paste0("mean\n(", round(mean_value, 0), ")")) +
  annotate('text', x=median_value, y=6000,
           # fontface='bold',
           label=paste0("median\n(", round(median_value, 0), ")")) +
  annotate('text', x=350, y=6000,
           # fontface='bold',
           label=paste0('(n = ', nrow(SNPs), ')'))

# Generated pdf is extemely large after importing to vector graphics software
# (large number of data points)
# Alternative: get the histogram summary table and then build plot

d <- ggplot_build(p)$data[[1]]
hist <- data.frame(x=d$x,
                   # xmin=d$xmin, xmax=d$xmax,
                   y=d$y)

ggplot(hist, aes(x, y)) +
  geom_bar(stat='identity', colour='grey30', alpha=0.25, na.rm = T) +
  geom_segment(aes(x=mean_value, y=0, xend=mean_value, yend=2800),
               linetype='dotted', color='red', size=0.4) +
  geom_segment(aes(x=median_value, y=0, xend=median_value, yend=4800),
               linetype='dotted', color='red', size=0.4) +
  labs(title='', x='Distance (bp)', y='Count') +
  annotate('text', x=mean_value, y=4500,
           label=paste0("mean\n(", round(mean_value, 0), ")")) +
  annotate('text', x=median_value, y=6500,
           label=paste0("median\n(", round(median_value, 0), ")")) +
  annotate('text', x=350, y=7000,
           label=paste0('(n = ', nrow(SNPs), ')'))
