#------------------------------------------------------------------------------#
#                                  ggplot2                                     #
#                                  settings                                    #
#------------------------------------------------------------------------------#

library(ggplot2)

# ggplot2 theme
ggplot2_theme <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=10),
        plot.subtitle=element_text(hjust=0.5, size=10),
        axis.text=element_text(colour='black'),
        axis.ticks=element_line(colour='black'))
theme_set(ggplot2_theme)
