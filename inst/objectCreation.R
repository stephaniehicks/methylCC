#' @title A reduced size of the FlowSorted.Blood.450k dataset 
#' 
#' @description A reduced size of the FlowSorted.Blood.450k dataset 
#' 
#' @docType data 
#' @format A RGset object with 2e5 rows (probes)
#' and 6 columns (whole blood samples).
#' 
library(FlowSorted.Blood.450k)
data(FlowSorted.Blood.450k)
# take a random sample to make object size in build smaller
set.seed(12345)
cpg_ids <- sample(seq_len(nrow(FlowSorted.Blood.450k)), 1.5e5)
FlowSorted.Blood.450k.sub <- FlowSorted.Blood.450k[cpg_ids,
            pData(FlowSorted.Blood.450k)$CellTypeLong %in% "Whole blood"]

## Save FlowSorted.Blood.450k.sub
save(FlowSorted.Blood.450k.sub, 
     file = "data/FlowSorted.Blood.450k.sub.RData", 
     compress='xz') 

#' @title Regions with low methylation values in six pure blood cell types
#'
#' @description
#' This script is used to create a ADD ME
#'
#' @docType data
#'
library(FlowSorted.Blood.450k)
library(minfi)
data(FlowSorted.Blood.450k)

Mset <- preprocessIllumina(updateObject(FlowSorted.Blood.450k))
Mset <- mapToGenome(Mset, mergeManifest = FALSE)

# Only consider autosomes
rowKeep <- which(!seqnames(granges(Mset)) %in% c("chrX","chrY"))
Mset <- Mset[rowKeep,]

# Based on MDS plots, we see sample 13 is an oulier; we remove it
IDs = c("CD8T","CD4T", "NK","Bcell","Mono","Gran")
colKeep <- which(pData(Mset)$CellType %in% IDs)
colKeep = setdiff(colKeep, 13)

# order samples by cell types (IDs)
colKeep = colKeep[order(factor(pData(Mset[, colKeep])$CellType, levels = IDs))]
Mset <- Mset[, colKeep]

pd = as.data.frame(pData(Mset))
pd$CellType <- factor(pd$CellType, levels = IDs)
p = getBeta(Mset, type = "Illumina")    # beta values
colnames(p) = pd$Sample_Name = rownames(pd) = gsub("\\+","",pd$Sample_Name)

## Create Genomic Ranges object
gr_train450k <- granges(Mset)
chr <- as.character(seqnames(gr_train450k))
pos <- start(gr_train450k)


## Calculate posterior probablity of being in group 1 (groups: 0, 1)
##      using Mclust and find d0 and d1 regions using p0 data
library(mclust)
probZ1 = sapply(seq_len(ncol(p)), function(x){
   out = Mclust(p[,x], G = 2)
   out$z[,2]
})
row.names(probZ1) <- NULL
colnames(probZ1) <- colnames(p)

Z0dat = regionFinder(1 / rowMax(probZ1), chr = chr, pos = pos,
                     cluster = cl, cutoff = 1 / 0.05 )
offMethRegions = Z0dat[Z0dat$L > 3,] # DMR regions with at least 3 probes

## Save offMethRegions
save(offMethRegions, file = "data/offMethRegions.RData", compress='xz') 


#' @title Regions with low methylation values in six pure blood cell types
#'
#' @description
#' This script is used to create a ADD ME
#'
#' @docType data
#' @format A MethylSet object with 1e4 rows (probes)
#' and 58 columns (samples).
#' 
Z1dat = regionFinder(1 / (rowMax(1-probZ1)+1e-10), chr = chr, pos = pos,
                     cluster = cl, cutoff = 1 / 0.05)
onMethRegions = Z1dat[Z1dat$L>3,] # DMR regions with at least 3 probes

## Save onMethRegions
save(onMethRegions, file = "data/onMethRegions.RData", compress='xz') 


