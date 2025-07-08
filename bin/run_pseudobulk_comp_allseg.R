set.seed(999)
library("fastknn")

suppressPackageStartupMessages({
    library(glue)
    library(tessera)
    library(scCustomize)
    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(sf)
    ## Plotting functions 
    ## Not imported by Tessera
    library(ggplot2)
    library(ggthemes)
    library(viridis)
    library(patchwork)
    library(harmony)
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(presto)
    library(tibble)
    library(scDotPlot)
    library(SingleR)    
    library(Seurat)
    library(presto)
    library(dplyr)
    library(ggplot2)
    library(argparse)
    library(tidyverse)
})


parser <- ArgumentParser(prog="run_annotation.r", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()

output = args$output
input_data = args$data
print(input_data)

xen.int <- readRDS(input_data[1])
xen.int2 <- readRDS(input_data[2])
singleR.int = readRDS(input_data[3])

nn1 = AggregateExpression(object = xen.int, group.by = c('tech', 'cell_type'), slot='data')$RNA
nn2 = AggregateExpression(object = xen.int2, group.by = c('tech', 'cell_type'), slot='data')$RNA
singleR = AggregateExpression(object = singleR.int, group.by = c('tech', 'fine_ids'), slot='data')$RNA


print(dim(nn1))

nn1 = log1p(as.matrix(nn1))
nn2 = log1p(as.matrix(nn2))
singleR = log1p(as.matrix(singleR))
print(head(singleR))
nn1.cor = c()
nn2.cor = c()
singleR.cor = c()

pdf(glue("output/{args$output}_Rplots.pdf"), width=18, height=6)
for (cel in unique(xen.int$cell_type)) {
     n1.cor <- cor(nn1[, glue("sc_{cel}")], nn1[, glue("xenium_{cel}")], method='spearman')
     n2.cor <- cor(nn2[, glue("sc_{cel}")], nn2[, glue("xenium_{cel}")], method='spearman')
     s.cor <- cor(singleR[, glue("sc_{cel}")], singleR[, glue("xenium_{cel}")], method='spearman')
     nn1.cor <- c(nn1.cor, n1.cor)
     nn2.cor <- c(nn2.cor, n2.cor)
     singleR.cor <- c(singleR.cor, s.cor)

     n1.cor2 <- cor(nn1[, glue("sc_{cel}")], nn1[, glue("xenium_{cel}")], method='pearson')
     n2.cor2 <- cor(nn2[, glue("sc_{cel}")], nn2[, glue("xenium_{cel}")], method='pearson')
     s.cor2 <- cor(singleR[, glue("sc_{cel}")], singleR[, glue("xenium_{cel}")], method='pearson')

corr_label1 <- paste0("r = ", round(n1.cor, 3))
corr_label2 <- paste0("r = ", round(n2.cor, 3))
corr_label3 <- paste0("r = ", round(s.cor, 3))

corr2_label1 <- paste0("r = ", round(n1.cor2, 3))
corr2_label2 <- paste0("r = ", round(n2.cor2, 3))
corr2_label3 <- paste0("r = ", round(s.cor2, 3))


df = data.frame(x=nn1[, glue("sc_{cel}")], y=nn1[, glue("xenium_{cel}")])
df2 = data.frame(x=nn2[, glue("sc_{cel}")], y=nn2[, glue("xenium_{cel}")])
df3 = data.frame(x=singleR[, glue("sc_{cel}")], y=singleR[, glue("xenium_{cel}")])

p0=ggplot(df3, aes(x = x, y = y)) +
  geom_point() +  # scatterplot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # regression line with confidence interval
  annotate("text", 
           x = max(df$x), 
           y = max(df$y), 
           label = glue("{corr_label1}\n{corr2_label1}"), 
           hjust = 1.1, vjust = 1.1, 
           size = 5) + 
  theme_minimal() +
  labs(title = glue("{cel} fgraph"),
       x = "Single cell KPMP",
       y = "Xenium")
# Scatterplot with linear regression line
p1=ggplot(df, aes(x = x, y = y)) +
  geom_point() +  # scatterplot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # regression line with confidence interval
  annotate("text", 
           x = max(df$x), 
           y = max(df$y), 
           label = glue("{corr_label2}\n{corr2_label2}"), 
           hjust = 1.1, vjust = 1.1, 
           size = 5) + 
  theme_minimal() +
  labs(title = glue("{cel} UMAP NN"),
       x = "Single cell KPMP",
       y = "Xenium")
p2=ggplot(df2, aes(x = x, y = y)) +
  geom_point() +  # scatterplot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # regression line with confidence interval
  annotate("text", 
           x = max(df2$x), 
           y = max(df2$y), 
           label = glue("{corr_label3}\n{corr2_label3}"), 
           hjust = 1.1, vjust = 1.1, 
           size = 5)+ 
  theme_minimal() +
  labs(title = glue("{cel} SingleR"),
       x = "Single cell KPMP",
       y = "Xenium")
print(p0+p1+p2)
}
dev.off()


