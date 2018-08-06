#------------------------------------------------------------------------------#
#                                                                              #
#                                  Figure S2                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(tidyverse)
library(wesanderson)

# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
wes_colors <- wes_palette(name = 'Zissou1', n=3, type = "continuous")
colors <- c('black', brewer.pal(9, "YlOrRd")[c(3, 5, 9)])
colors <- c('black', wes_colors)


#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Spike-in normalization factors using read counts versus bedtools read pileup #
#------------------------------------------------------------------------------#
# Load data table
sinfs <- read_csv(
  here('data/spikein_normalization_factors_using-different_methods.csv'))
head(sinfs)

data <- sinfs
data <- subset(sinfs, !sample %in% c(
    'AH8583_1', 'AH8583_2', 'AH8584_1', 'AH8584_2', 'AH9120_1', 'AH9120_2',
    'AH8115_1', 'AH8115_2'))

# Rename methods
renamed_data <- mutate(
  data,
  method=recode(method,
                read_counts='Read counts', read_pileups='Read pileups',
                read_pileups_on_SNPs='Read pileups on SNPs',
                read_pileups_on_SNPs_in_peaks='Read pileups on SNPs in peaks'))

# Order methods
renamed_data$method <- factor(
  renamed_data$method, levels=c(
    'Read counts', 'Read pileups', 'Read pileups on SNPs',
    'Read pileups on SNPs in peaks'))

# Order samples by decreasing Red1 amount
renamed_data$sample <- factor(
    renamed_data$sample, levels=c(
        'AH8104_1', 'AH8104_2', 'AH8218_1', 'AH8218_2', 'AH8220_1', 'AH8220_2',
        'AH8151_1', 'AH8151_2', 'AH9048_1', 'AH9048_2', 'AH7011_1', 'AH7011_2',
        'AH8219_1', 'AH8219_2'))

ggplot(renamed_data, aes(sample, sinf, fill=method)) +
  scale_fill_manual('Data type', values=colors) +
  geom_bar(stat="identity", colour='white', position="dodge") +
  ylim(0, 1) +
  labs(title='', x='',
       y='Red1 amount\nrelative to wild type') +
  scale_x_discrete(labels=c(
      expression(
          paste(italic('dot1'), Delta, ' 1'),
          paste(italic('dot1'), Delta, ' 2'),
          paste(
              italic('red1'[italic('ycs4S')]), '/', italic('RED1'), ' 1'),
          paste(
              italic('red1'[italic('ycs4S')]), '/', italic('RED1'), ' 2'),
          paste(italic('red1'), Delta, '/RED1 1'),
          paste(italic('red1'), Delta, '/RED1 2'),
          paste(italic('rec8'), Delta, ' 1'),
          paste(italic('rec8'), Delta, ' 2'),
          paste(italic('red1-pG162A'), ' 1'),
          paste(italic('red1-pG162A'), ' 2'),
          paste(italic('red1'[italic('ycs4S')]), ' 1'),
          paste(italic('red1'[italic('ycs4S')]), ' 2'),
          paste(
              italic('red1'[italic('ycs4S')]), '/', italic('red1'), Delta, ' 1'),
          paste(
              italic('red1'[italic('ycs4S')]), '/', italic('red1'), Delta, ' 2'))
      )) +
  theme(axis.text.x=element_text(face='italic', angle=45, hjust=1),
        legend.text=element_text(size=8))
  
#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# Correlations                                                                 #
#------------------------------------------------------------------------------#
spread_data <- spread(data, key=method, value=sinf)

plot_correlation <- function(data, x, y, x_label, y_label) {
  m <- lm(as.formula(paste(y, '~', x)), data)
  r2 <- round(summary(m)$r.squared, 3)
  
  ggplot(data, aes_string(x, y)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE, colour='red') +
    xlim(0, 1) + ylim(0, 1) +
    geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
    labs(title='', x=x_label, y=y_label) +
    annotate('text', x=0.30, y=0.80, size=4,
             label=paste('italic(R)^{2}==', r2), parse=TRUE)
}

plot_correlation(spread_data, 'read_counts', 'read_pileups',
                 'Read counts', 'Read pileups')
plot_correlation(spread_data, 'read_counts', 'read_pileups_on_SNPs',
                 'Read counts', 'Read pileups\noverlapping SNPs')
plot_correlation(spread_data, 'read_counts', 'read_pileups_on_SNPs_in_peaks',
                 'Read counts', 'Read pileups overlapping\nSNPs in peaks')
