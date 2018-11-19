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


#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# Histogram of distances between consecutive SNPs for each chromosome          #
#------------------------------------------------------------------------------#

mean_values <- SNPs %>% group_by(chr) %>%
    summarise(Mean=mean(Distance, na.rm = TRUE))
median_values <- SNPs %>% group_by(chr) %>%
    summarise(Median=median(Distance, na.rm = TRUE))

drop_genome <- function(string) str_replace(string, '_SK1', '')

p <- ggplot(SNPs, aes(Distance)) +
    geom_histogram(aes(y= ..count..), bins=50, colour='grey30',
                   alpha=0.25, na.rm = T) +
    facet_wrap( ~ chr, ncol=4, labeller=as_labeller(drop_genome)) +
    xlim(-5, 500) +
    geom_segment(aes(x=Mean, y=0, xend=Mean, yend=700), data=mean_values,
                 linetype='dotted', color='red', size=0.4) +
    geom_text(aes(x=Mean+80, y=800, label=paste0("mean (", round(Mean, 0), ")")),
              data=mean_values, size=3) +
    geom_segment(aes(x=Median, y=0, xend=Median, yend=1000), data=median_values,
                 linetype='dotted', color='red', size=0.4) +
    geom_text(aes(x=Median+80, y=1100,
                  label=paste0("median (", round(Median, 0), ")")),
                  data=median_values, size=3) +
    labs(title='', x='Distance (bp)', y='Count') +
    theme(
        strip.background=element_blank()
    )

p
