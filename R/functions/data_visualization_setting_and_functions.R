## ggplot setting
library("ggplot2")

theme_classic_2 <- function(base_size = 16){
    theme_classic() %+replace%
        theme(
            axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black")
        )
}

theme_set(theme_classic_2())
grDevices::pdf.options(useDingbats = FALSE)

## knitr setting
library("knitr")
opts_chunk$set(dev = c("png", "pdf"))
