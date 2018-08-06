#------------------------------------------------------------------------------#
#                                                                              #
#                                   Figure 2                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(wesanderson)
library(tidyverse)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load IO code
source(here('helper_io.R'))
# Load code to compute spike-in normalization factor
source(here('helper_spikein_normalization_factor.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Plot colors -----------------------------------------------------------------#
wes_colors <- wes_palette(name = "Darjeeling1", n=20, type = "continuous")
wt_colors <- c('grey30', 'grey50')
low_Red1_color <- rev(wes_colors[17:18])

wt_color <- 'black'


#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Sequencing depth versus aligned reads                                        #
#------------------------------------------------------------------------------#
# Load read alignment counts
list.files(Subsample_aligned_read_counts)

read_all_counts <- function(dir) {
  all_files <- list()
  for (file in list.files(dir)) {
    data <- suppressMessages(
      read_tsv(file.path(dir, file), col_names = c('chr', 'count'), n_max=32)
    )
    data_name <- sub("stats_", "", file)
    data_name <- sub(".txt", "", data_name)
    
    all_files[[data_name]] <- data
  }
  
  all_files
}

all_data <- read_all_counts(dir=Subsample_aligned_read_counts)

df <- data.frame(Genome=character(0), Count=integer(0),
                 Sample=character(0), read_n=integer(0))
for (i in seq_along(all_data)) {
  data_names <- strsplit(names(all_data)[i], "_")[[1]]
  data <- sum_per_genome(all_data[[i]])
  
  df <- rbind(df, data.frame(data.frame(Genome=c('S288C', 'SK1'),
                                        Count=c(data['S288C'], data['SK1'])),
                             Sample=paste(data_names[1],
                                          data_names[2], sep='_'),
                             read_n=as.numeric(sub('Mreads', '',
                                                   data_names[3]))))
}

df <- subset(df, !Sample %in% c('AH9048_input', 'AH9048_ip'))


df <- mutate(
  df,
  Sample=recode(Sample,
                AH119_input='RED1 input', AH119_ip='RED1 IP',
                AH7011_input='red1_ycs4 input', AH7011_ip='red1_ycs4 IP'))

ggplot(df, aes(read_n, Count / 1000000, colour=Sample, shape=Genome)) +
  stat_smooth(method = "lm", size=0.5, se = FALSE) +
  geom_point(size=3) +
  scale_shape_manual('', values = c(1, 16)) +
  scale_color_manual('', values=c(wt_colors, low_Red1_color),
                     labels=c(expression(paste(italic('RED1'), ' input'),
                                         paste(italic('RED1'), ' IP'),
                                         paste(italic('red1'[italic('ycs4S')]),
                                               ' input'),
                                         paste(italic('red1'[italic('ycs4S')]),
                                               ' IP')))) +
  scale_x_discrete(limits=1:10) +
  labs(y=paste0('Number of aligned reads\n(x 1 million)'),
       x=paste0('FASTQ subsample size\n(x 1 million)')) +
  theme(legend.text.align = 0)

#------------------------------------------------------------------------------#
#                                   Panel B                                    #
# All combinations of sequencing depths                                        #
#------------------------------------------------------------------------------#
# Compute SINF for all combinations of FASTQ subsamples
sinf_for_all_combinations <- function() {
  t0 <- proc.time()[3]
  path <- Subsample_aligned_read_counts
  
  depth_unbalance_sinf <- data.frame()
  message('Outer loop iteration (progress tracking): ', appendLF=FALSE)
  for (i in 1:10) {
    message(i, ' ', appendLF=FALSE)
    for (j in 1:10) {
      for (k in 1:10) {
        for (l in 1:10) {
          suppressMessages(
            sinf <- spikein_normalization_factor_from_counts(
              ref_chip_counts=file.path(path, paste0('stats_AH119_ip_', i,
                                                     'Mreads.txt')),
              ref_input_counts=file.path(path, paste0('stats_AH119_input_', j,
                                                      'Mreads.txt')),
              test_chip_counts=file.path(path, paste0('stats_AH7011_ip_', k,
                                                      'Mreads.txt')),
              test_input_counts=file.path(path, paste0('stats_AH7011_input_', l,
                                                       'Mreads.txt')))
          )
          
          depth_unbalance_sinf <- bind_rows(
            depth_unbalance_sinf,
            data.frame(sinf=sinf, n_ref_input=i, n_ref_ChIP=j,
                       n_test_input=k, n_test_ChIP=l))
        }
      }
    }
  }
  
  message()
  message('---')
  message('Completed in ', elapsed_time(t0, proc.time()[3]))
  
  return(depth_unbalance_sinf)
}

depth_unbalance_sinf <- sinf_for_all_combinations()

mean(depth_unbalance_sinf$sinf)
sd(depth_unbalance_sinf$sinf)

ggplot(depth_unbalance_sinf, aes(x = sinf)) +
  geom_histogram(aes(y= ..count..), bins=500,  colour='grey30', alpha=0.25) +
  xlim(0, 1) +
  labs(title='', subtitle='',
       x='Spike-in normalization factor',
       y='Count')

ggplot(depth_unbalance_sinf, aes(x = sinf)) +
  geom_histogram(aes(y= ..count..), bins=50,  colour='black', alpha=0.25) +
  xlim(0.27, 0.30) +
  # geom_vline(xintercept=mean(depth_unbalance_sinf$sinf),
  #            linetype="dashed", color = "red") +
  labs(x='Spike-in normalization factor', y='Count')


#------------------------------------------------------------------------------#
#                                   Panel D                                    #
# Spike-in normalization factor                                                #
#------------------------------------------------------------------------------#
# Make data frame of spike-in normalization factors

# Locate samples
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah119spike-', ignore_case=TRUE))]
names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                      fixed('ah9048spike', ignore_case=TRUE))]

nf <- data.frame(
  sample=c('5%', '10%', '15%', '20%', '25%', '30%'),
  spinf=c(suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-5-95_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-5-95_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-5-95_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-5-95_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-10-90_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-10-90_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-10-90_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-10-90_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-15-85_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-15-85_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-15-85_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-15-85_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-20-80_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-20-80_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-20-80_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-20-80_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-25-75_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-25-75_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-25-75_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-25-75_trimmed']])),
          suppressMessages(
            spikein_normalization_factor_from_counts(
              ref_chip_counts=
                read_counts_per_chr[['ah119spike-chip-30-70_trimmed']],
              ref_input_counts=
                read_counts_per_chr[['ah119spike-input-30-70_trimmed']],
              test_chip_counts=
                read_counts_per_chr[['ah9048spike-chip-30-70_trimmed']],
              test_input_counts=
                read_counts_per_chr[['ah9048spike-input-30-70_trimmed']]))))


# Order samples
nf$sample <- factor(nf$sample,
                    levels = c('5%', '10%', '15%', '20%', '25%', '30%'))

nf_mean <- round(mean(nf$spinf), 3)
nf_sd <- round(sd(nf$spinf), 3)

ggplot(nf, aes(sample, spinf * 100)) +
  geom_bar(stat="identity", width=0.5, alpha=1, colour=NA, fill='black') +
  scale_colour_manual("", values=strain_colors,guide=FALSE) +
  geom_text(aes(label=round(spinf, 3) * 100), vjust=-0.5, size=3) +
  ylim(0, 50) +
  labs(title='',
       x='Spike-in cell percentage',
       y='Red1 amount\n(% of wild type)')


#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# % spike-in cells vs % spike-in reads in Input and ChIP                       #
#------------------------------------------------------------------------------#

get_file_tetrad <- function(pattern) {
  names(read_counts_per_chr)[str_detect(read_counts_per_chr,
                                        fixed(pattern, ignore_case=TRUE))]
}

data <- list('5:95'=get_file_tetrad('5-95'),
             '10:90'=get_file_tetrad('10-90'),
             '15:85'=get_file_tetrad('15-85'),
             '20:80'=get_file_tetrad('20-80'),
             '25:75'=get_file_tetrad('25-75'),
             '30:70'=get_file_tetrad('30-70'))

counts <- lapply(data, function(x) suppressMessages(
  spikein_normalization_factor_from_counts(
    ref_chip_counts=read_counts_per_chr[[x[1]]],
    ref_input_counts=read_counts_per_chr[[x[2]]],
    test_chip_counts=read_counts_per_chr[[x[3]]],
    test_input_counts=read_counts_per_chr[[x[4]]],
    return_counts=TRUE)))

proportions <- data.frame(sample=character(),
                          percent=numeric())
for (i in seq_along(counts)) {
  data <- counts[[i]]
  message('Proportion: ', names(counts)[i])
  new_proportions <- data.frame(sample=c(
    paste0('AH119_chip_', names(counts)[i]),
    paste0('AH119_input_', names(counts)[i]),
    paste0('AH9048_chip_', names(counts)[i]),
    paste0('AH9048_input_', names(counts)[i])),
    percent=c(
      data[['ref_chip']]['S288C'] / (data[['ref_chip']]['S288C'] 
                                     + data[['ref_chip']]['SK1']),
      data[['ref_input']]['S288C'] / (data[['ref_input']]['S288C'] 
                                      + data[['ref_input']]['SK1']),
      data[['test_chip']]['S288C'] / (data[['test_chip']]['S288C'] 
                                      + data[['test_chip']]['SK1']),
      data[['test_input']]['S288C'] / (data[['test_input']]['S288C'] 
                                       + data[['test_input']]['SK1'])))
  
  proportions <- rbind(proportions, new_proportions)
}


# Add additional variables to facet plot
reorganized_table <- data.frame(
  sample=proportions[str_detect(proportions$sample, 'chip'), 'sample'],
  ChIP=proportions[str_detect(proportions$sample, 'chip'), 'percent'] * 100,
  Input=proportions[str_detect(proportions$sample, 'input'), 'percent'] * 100)

reorganized_table$strain <- NA
reorganized_table[
  str_detect(reorganized_table$sample, 'AH119'), 'strain'] <- 'RED1'
reorganized_table[
  str_detect(reorganized_table$sample, 'AH9048'), 'strain'] <- 'red1-pG162A'

reorganized_table$spikein_proportion <- NA
reorganized_table[
  str_detect(reorganized_table$sample, '5:95'), 'spikein_proportion'] <- 5
reorganized_table[
  str_detect(reorganized_table$sample, '10:90'), 'spikein_proportion'] <- 10
reorganized_table[
  str_detect(reorganized_table$sample, '15:85'), 'spikein_proportion'] <- 15
reorganized_table[
  str_detect(reorganized_table$sample, '20:80'), 'spikein_proportion'] <- 20
reorganized_table[
  str_detect(reorganized_table$sample, '25:75'), 'spikein_proportion'] <- 25
reorganized_table[
  str_detect(reorganized_table$sample, '30:70'), 'spikein_proportion'] <- 30

reorganized_table

plot_correlation <- function(data, x, y, plot_title,
                             annotation_position=c(30, 10)) {
  m <- lm(as.formula(paste(y, '~', x)), data)
  r2 <- round(summary(m)$r.squared, 3)
  
  ggplot(data, aes_string(x, y)) +
    geom_smooth(method='lm', se=FALSE, colour='red', fullrange=TRUE) +
    geom_point() +
    scale_y_continuous(breaks=seq(0, 100, 10)) +
    scale_x_continuous(breaks=seq(0, 100, 10)) +
    labs(title=plot_title, x='Input', y='ChIP') +
    expand_limits(x = 0, y = 0) +
    annotate('text', x=annotation_position[1], y=annotation_position[2],
             size=3.5, label=paste("italic(R)^{2}==", r2), parse=TRUE)
}

plot_correlation(subset(reorganized_table, strain == 'RED1'),
                 'Input', 'ChIP', expression(italic('RED1')),
                 annotation_position=c(10, 27.5))
plot_correlation(subset(reorganized_table, strain == 'red1-pG162A'),
                 'Input', 'ChIP', expression(italic('red1-pG162A')),
                 annotation_position=c(10, 55))
