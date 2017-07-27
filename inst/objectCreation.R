#' @name offMethRegions
#' @title Regions with low methylation values in six pure blood cell types
#'
#' @description
#' This script is used to create a ADD ME
#'
#' @docType data
#' @keywords datasets
#' @format A MethylSet object with 1e4 rows (probes)
#' and 58 columns (samples).
#'
#' @rdname offMethRegions
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

## Save offMethRegions to be used in HousemanSeq()
save(offMethRegions, file = "data/offMethRegions.RData")



#' @name onMethRegions
#' @title Regions with low methylation values in six pure blood cell types
#'
#' @description
#' This script is used to create a ADD ME
#'
#' @docType data
#' @keywords datasets
#' @format A MethylSet object with 1e4 rows (probes)
#' and 58 columns (samples).
#'
#' @rdname onMethRegions
#'
Z1dat = regionFinder(1 / (rowMax(1-probZ1)+1e-10), chr = chr, pos = pos,
                     cluster = cl, cutoff = 1 / 0.05)
onMethRegions = Z1dat[Z1dat$L>3,] # DMR regions with at least 3 probes

## Save onMethRegions to be used in HousemanSeq()
save(onMethRegions, file = "data/onMethRegions.RData")


#' @name celltypeSpecificDMRegions
#' @title Cell type specific regions in four blood cell types
#'
#' @description
#' This script is used to create a ADD ME
#'
#' @docType data
#' @keywords datasets
#' @format
#'
#' @rdname celltypeSpecificDMRegions
#' @export
#'
library(rafalib)
# First: Group cell types (CD4T, CD8T and NK) into Tcell.
#        See SVD plot which shows 3 cell types are very similar
oldcelltypes <- pd$CellType
colcelltypesID <- unique(pd$CellType)
s <- svd(p-rowMeans(p))
mypar()
plot(s$v[,1:2],col=as.fumeric(as.character(oldcelltypes)),pch=16)
legend("bottomright",levels(oldcelltypes),col=1:6,pch=16)

# define new cell types
cell <- as.character(pd$CellType)
cell[cell%in%c("CD4T","CD8T","NK")] <- "Tcell"
cells <- unique(cell)

# Second: Find celltype specific DMRs
library(limma)
library(bumphunter)
cl <- clusterMaker(chr, pos)

allPoss = as.matrix(expand.grid(c(0,1), c(0,1), c(0,1), c(0,1)))
Z = allPoss = allPoss[-c(1, 16),]
allPoss <- (allPoss == TRUE)


Xmat = cbind(rep(1, length(cell)), model.matrix(~factor(cell, levels = cells) - 1))
colnames(Xmat) = c("Intercept", cells)


regions = zmat <- c() # regions will contain all celltype-specific DMRs
for(ind in seq_len(nrow(allPoss))){ # length 14

    # subset object and Xmat for pairwise design matrix
    subXmat = cbind("Intercept" = Xmat[, "Intercept"],
                    "cellTypes" = rowSums(as.matrix(Xmat[, cells[allPoss[ind,]] ],
                                                    ncols = length(cells[allPoss[ind,]]))))
    # find bumps
    bumps = bumphunter(object = p, design = subXmat, chr = chr,
                       pos = pos, cluster = cl, cutoff = 0.2, B = 0,
                       smooth = FALSE, smoothFunction = loessByCluster)

    # y are the beta values collapsed by regions from bumphunter
    y <- t(apply(bumps$table[,7:8], 1, function(z){
                    colMeans(p[(z[1]):(z[2]),,drop=FALSE]) } ))

    # Only include region with more than 1 CpG (L > 1)
    #       OR only 1 CpG in region if no other larger regions possible
    keepInd <- which( (bumps$table$L > 1 | (bumps$table$L==1 & bumps$table$clusterL == 1)) &
                        apply(abs(sweep(y, 2, subXmat[,2], FUN = "-")) < 0.4, 1, all) )
    cat(min(length(keepInd), 25), "regions found.\n")
    if(nrow(bumps$table[keepInd,]) > 25){
            regions <- rbind(regions, head(bumps$table[keepInd,], n = 25))
            zmat <- rbind(zmat, t(replicate(25, Z[ind,])))
    } else {
        regions <- rbind(regions, bumps$table[keepInd,])
        if(nrow(bumps$table[keepInd,]) > 0){
            zmat <- rbind(zmat, t(replicate(nrow(bumps$table[keepInd,]), Z[ind,])))
        }
    }
}
colnames(zmat) <- cells
zmat

y <- t(apply(regions[,7:8],1,function(ind){
    colMeans(p[(ind[1]):(ind[2]),,drop=FALSE])
}))

# quick check of regions
mypar()
for(i in sample(nrow(y),25)){
    plot(y[i,],col=as.fumeric(cell),main=i,ylim=c(0,1))
    abline(h=0.5)
}

###more visual inspection
heatmap(y,labCol=cell,scale="none", Rowv = NA)
Z[,c(3,4,2,1)]

###create a profile for each cell type
profiles <- sapply(splitit(cell),function(ind)
    rowMeans(y[,ind]))


celltypeSpecificDMRegions <- makeGRangesFromDataFrame(regions)
mcols(celltypeSpecificDMRegions) <- zmat
save(celltypeSpecificDMRegions, file = "data/celltypeSpecificDMRegions.RData")




#### Creating WGBS object

library(bsseq)
library(readr)
load(file.path("/n/irizarryfs01_backed_up/kkorthauer/WGBS/MCGILL/DATA", "WB_BSSeq.RData")) # loads WB object
dat <- read_delim("/n/irizarryfs01_backed_up/kkorthauer/WGBS/MCGILL/sampleInfo.txt",
                  delim = "\t")
dat <- dat[dat[,9] %in% "blood",]


load("data/WB_BSSeq.RData")

