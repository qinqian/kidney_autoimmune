devtools::install_github("jinworks/CellChat")
install.packages("sf")

library(CellChat)
library(patchwork)
library(CellChat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(sf)
library(RColorBrewer)

ptm = Sys.time()

data <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")


#Part I: Data input & processing and initialization of CellChat object

#Laod data
tile_obj <- data$tile_obj #Extract object
obj <- data$obj

sample_to_analyze <- "BS21-N65682A2__20241025" #Subset sample
obj_sub <- subset(obj, subset = sample_id == sample_to_analyze)

tile_meta <- tile_obj@meta.data #Metadata
cell_meta <- obj_sub@meta.data

## merged_meta <- merge(
##   cell_meta,
##   tile_meta[, c("id", "niche_label")],
##   by.x = "tile_id",   
##   by.y = "id",         
##   all.x = TRUE
## )
## obj_sub@meta.data$niche_label <- merged_meta$niche_label
## Idents(obj_sub) <- "niche_labels"  


data.input <- GetAssayData(obj_sub, layer = "data", assay = "RNA")

meta <- data.frame(
  labels = factor(Idents(obj_sub)),  
  row.names = colnames(obj_sub)
)

meta <- meta[!is.na(meta$labels), , drop = FALSE]
meta$labels <- droplevels(meta$labels)

coords <- obj_sub@meta.data[, c("X", "Y")]
rownames(coords) <- rownames(obj_sub@meta.data)

common_cells <- Reduce(intersect, list(colnames(data.input), rownames(meta), rownames(coords)))
data.input <- data.input[, common_cells]
meta       <- meta[common_cells, , drop = FALSE]
coords     <- coords[common_cells, , drop = FALSE]

spatial.factors <- data.frame(
  ratio = rep(1, nrow(coords)),
  tol = rep(1, nrow(coords))
)
rownames(spatial.factors) <- rownames(coords)


#Create a CellChat object
cellchat <- createCellChat(
  object = data.input,      # your normalized expression matrix
  meta = obj_sub@meta.data,              # includes "labels"
  group.by = "celltype",      # what to cluster/group by
  datatype = "spatial",     # critical: tells CellChat to use spatial modeling
  coordinates = as.matrix(coords) ,     # your custom X/Y
  spatial.factors = spatial.factors
)

#Set the ligand-receptor interaction database
library(CellChat)
data(CellChatDB.human)
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
devtools::install_github("immunogenomics/presto")

cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

color_vec <- c(
  "TAL" = "#E41A1C",                           
  "Immune infiltrated interstitial" = "#FF7F00",  
  "Proximal Tubule" = "#4DAF4A",               
  "Fibrogenic interstitial" = "#984EA3",       
  "Collecting Duct" = "#377EB8",               
  "Glomeruli" = "#FFD700",                     
  "Epithelial+Muscle(Vessel)" = "#FF69B4",     
  "Muscle Outlier" = "#A65628"                 
)

ggplot(tile_sf) +
  geom_sf(aes(fill = niche_label), color = "black", size = 0.05) +  # adds visible thin borders
  scale_fill_manual(values = color_vec) +
  theme_void() +
  ggtitle(paste("Tissue Niches:", sample_to_analyze)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )

#Or
plot(tile_sf["niche_label"])




#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
ptm = Sys.time()

library(BiocNeighbors)

cellchat <- subsetData(cellchat)

options(future.globals.maxSize = 2 * 1024^3)
cellchat <- computeCommunProb(
  cellchat, 
  raw.use = TRUE,  
  ## type = "truncatedMean", 
  trim = 0.1,
  distance.use = FALSE,           
  interaction.range = 250,      
  scale.distance = NULL,
  contact.dependent = TRUE,      
  contact.range = 100            
)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 1))
par(mar = c(2, 2, 2, 2))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = rowSums(cellchat@net$count),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)

par(mfrow = c(1, 1))
par(mar = c(2, 2, 2, 2))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = rowSums(cellchat@net$weight),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

# Signaling Visualization
pathways.show <- c("TGFb")

par(mfrow = c(1.5, 1), xpd = TRUE)   #Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Generate centrality heatmap 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
par(mfrow = c(1, 1))
netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways.show,
  width = 23,
  height = 15,
  font.size = 10
)

netAnalysis_contribution(cellchat, signaling = "LT", return.data = FALSE)

#Inspect MON PTS1/S2 interaction
# View significant interactions from MON to PT-S1/S2
net <- subsetCommunication(cellchat, sources.use = "MON", targets.use = "PT-S1/S2")
head(net)
