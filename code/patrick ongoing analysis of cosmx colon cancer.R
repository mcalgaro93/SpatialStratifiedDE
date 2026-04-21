# for other sections & data formats, see 
# https://zenodo.org/records/15574384

# necessary libraries
library(data.table) # for more memory efficient data frames
library(Matrix) # for sparse matrices like our counts matrix
library(ggplot2)
library(InSituCor) # devtools::install_github("https://github.com/Nanostring-Biostats/InSituCor")
library(HieraType) # remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", subdir = "_code/HieraType", ref = "Main")
library(Seurat)
library(SpatialExperiment)
library(DelayedArray)
library(Matrix)


if (FALSE) {
  # retrieval
  options(timeout=1e3)
  url <- "https://zenodo.org/records/15574384/files/221.zip"
  tf <- tempfile(fileext=".zip")
  download.file(url, tf)
  # decompression
  pa <- dirname(tf)
  unzip(tf, exdir=pa)
  # importing
  library(alabaster.sce)
  nm <- file.path(pa, "221")
  sce <- readObject(nm)
  # coersion
  xy <- sprintf("Center%s_global_mm", c("X", "Y"))
  xy <- as.matrix(colData(sce)[xy])
  (spe <- toSpatialExperiment(sce, spatialCoords=xy))
  
  #saveRDS(spe, file = "data/spe.RDS")
  
  
  ## counts matrix (genes x cells)
  counts <- counts(spe)
  counts <- DelayedArray(counts)
  counts <- as(counts, "dgCMatrix")
  
  ## cell (column) metadata
  metadata <- as.data.frame(colData(spe))
  
  ## spatial coordinates (cells x dims, usually 2 columns)
  xy <- as.matrix(spatialCoords(spe))
  
  ## optional: make sure dimensions align
  stopifnot(ncol(counts) == nrow(metadata),
            ncol(counts) == nrow(xy))
  
  save(counts, metadata, xy, file = "basicdata.RData")
}

load("basicdata.RData")
counts <- Matrix::t(counts)

# normalize:
scale_row <- mean(metadata$nCount_RNA) / metadata$nCount_RNA
norm <- counts
norm@x <- norm@x * scale_row[norm@i + 1L]

##### explore:
granular <- as.vector(metadata$lv2)
granular[is.na(granular)] <- "unknown"

# big clusts:
clust <- granular
clust[grepl("epi", clust)] <- "epi"
clust[grepl("fib", clust)] <- "fib"
clust[grepl("CAF", clust)] <- "CAF"
clust[grepl("B\\.", clust)] <- "B"
clust[grepl("PC", clust)] <- "PC"
clust[grepl("SMC", clust)] <- "SMC"
clust[grepl("DC", clust)] <- "DC"

cellcols <- rep(c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3",
  "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
  "#B15928", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072"
), 3)[1:length(unique(clust))]
names(cellcols) <- unique(clust)

par(mar = c(0,0,0,0))
plot(xy, asp = 1, cex = 0.1, pch = 16,
     col = cellcols[clust])
legend("topleft", pch = 16, col = cellcols, legend = names(cellcols))
#points(xy[is.element(clust, c("BC", "TC")), ], col = 1, cex = 0.2, pch = 16)


#### define "Z" / spatial context: ----------------------------

### NO: too complicated:
# embed cells via scPearson:
# get neighbors:
# random sketch over neighbors:

## simple Z: just count cell types:

# get neighbors:
if (FALSE) {
  neighbors <- InSituCor:::nearestNeighborGraph(x = xy[,1], y = xy[,2], N = 50, subset = 1)
  saveRDS(neighbors, file = "data/neighbors.RDS")
} else {
  neighbors <- readRDS(file = "data/neighbors.RDS")
}

# count neighboring cells:
if (FALSE) {
  env <- InSituCor:::neighbor_tabulate(clust, neighbors)
  saveRDS(env, file = "data/env.RDS")
} else {
  env <- readRDS(file = "data/env.RDS")
}


#### define design variable: ----------------------------------------

inds.tumor <- clust == "epi"
inds.t     <- clust == "Tc"
stopifnot(any(inds.tumor), any(inds.t))
dist2tcell <- FNN::get.knnx(
  data  = as.matrix(xy[inds.t, , drop=FALSE]),
  query = as.matrix(xy),
  k = 1
)$nn.dist[, 1]

#### define patches ------------------------------------------------------

#source(url("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/refs/heads/Main/_code/PatchDE/DEutils.R"))
source(url("https://raw.githubusercontent.com/mcalgaro93/SpatialStratifiedDE/refs/heads/debug-getPatches/code/patchDE.R"))

# only use tumor cells kinda near T-cells:
use <- inds.tumor & (dist2tcell < 0.1)

patches <- getPatches(xy = xy[use, ], 
                      X = dist2tcell[use], 
                      npatches = 1500,
                      bitesize = 0.025,
                      maxradius = 0.2,
                      roundness = 0.25,
                      initwithhotspots = TRUE,
                      n_iters = 25,
                      plotprogress = FALSE
) # produces a vector of patch assigments, potentially including NA's for cells not given a patch
saveRDS(patches, file = "processed/patches1500.RDS")


plot(xy[use, ], pch = 16, cex = 0.1, col = "grey80")
points(xy[use, ], pch = 16, cex = 0.1, col = rep(cellcols,100)[as.numeric(as.factor(patches))])


####  run DE  -----------------------------------------------------

res <- patchDE(y = norm[use, ],
               df = data.frame("dist2tcell" = dist2tcell[use]),
               patch = patches)
str(res)