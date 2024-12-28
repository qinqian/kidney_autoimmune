library(ggplot2)
library(patchwork)
library(tidyverse)
library(tidyplots)

custom_colors <- c(
  "minisvl+tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+gs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:gs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+g" = rgb(249, 134, 130, maxColorValue = 255),
  "msv:g" = rgb(249, 134, 130, maxColorValue = 255),
  "msv:tg" = rgb(249, 134, 130, maxColorValue = 255),
  "nanomonsv" = rgb(147, 203, 118, maxColorValue = 255),
  "savana" = rgb(22, 183, 139, maxColorValue = 255),
  "severus" = rgb(16, 174, 228, maxColorValue = 255),
  "severus_lowaf" = rgb(16, 174, 228, maxColorValue = 255),
  "sniffles" = rgb(154, 130, 251, maxColorValue = 255),
  "svision" = rgb(248, 89, 206, maxColorValue = 255)
)

get_theme <- function(size=12, angle=0) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data = read_tsv(input[['tabs']]) 
    celltype <- data %>% select(celltype_hint2) %>% distinct() %>% pull()

    pdf(out_path[["fig"]], width=12, height=5)
    for (cell in celltype) {
	    print(cell)
	    print(data %>% filter(celltype_hint2==cell) )
	    print(dim(data %>% filter(celltype_hint2==cell) ))
	    print(ggplot(data=data %>% filter(celltype_hint2==cell)) +
		          geom_line(aes(x=distance, y=probability, color=group, shape=sample)) +
		          geom_point(aes(x=distance, y=probability, color=group, shape=sample)) + theme_classic(base_size=15)+ggtitle(cell))
	    print(data %>% filter(celltype_hint2==cell) %>% tidyplot(x=distance, y=probability, color=group) %>%
                  add_mean_line() %>%
                  add_mean_dot() %>%
                  add_sem_ribbon()+ggtitle(cell))
	    print(data %>% filter(celltype_hint2==cell) %>% tidyplot(x = group, y=probability, color=group) %>% 
		  add_boxplot() %>% add_test_pvalue(ref.group = 1)+ggtitle(cell))
    }
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)

