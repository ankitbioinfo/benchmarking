Sys.setenv(LANG="en")

suppressMessages(library(MERINGUE))

data(BCL) # object contains counts and positions for all 4 sections. Will only need sections 1-3
head(BCL$pos)
BCL$counts[1:5,1:5]

## parse out the 3 sections
BCL1.pos <- BCL$pos[BCL$pos$slice==1,][,c("x", "y")]
BCL2.pos <- BCL$pos[BCL$pos$slice==2,][,c("x", "y")]
BCL3.pos <- BCL$pos[BCL$pos$slice==3,][,c("x", "y")]
BCL1.counts <- BCL$counts[,BCL$pos$slice==1]
BCL2.counts <- BCL$counts[,BCL$pos$slice==2]
BCL3.counts <- BCL$counts[,BCL$pos$slice==3]

# Get common set of genes for the sections
genes.have <- Reduce(intersect, list(
  rownames(BCL1.counts),
  rownames(BCL2.counts),
  rownames(BCL3.counts)
))

# Combine into large counts matrix. Counts of genes for all the sections being assessed
counts <- cbind(
  BCL1.counts[genes.have,],
  BCL2.counts[genes.have,],
  BCL3.counts[genes.have,]
)

# section factor
section <-  c(
  rep('L1', ncol(BCL1.counts)),
  rep('L2', ncol(BCL2.counts)),
  rep('L3', ncol(BCL3.counts))
)
names(section) <- colnames(counts)

# List of positions
posList <- list(
  BCL1.pos[colnames(BCL1.counts),],
  BCL2.pos[colnames(BCL2.counts),],
  BCL3.pos[colnames(BCL3.counts),]
)

cc <- cleanCounts(counts, min.reads = 100, min.lib.size = 100, plot=TRUE)
mat <- normalizeCounts(cc, log=FALSE, verbose=TRUE)

## Normalizing matrix with 767 cells and 4985 genes.

## normFactor not provided. Normalizing by library size.

## Using depthScale 1e+06

posList[[1]] <- posList[[1]][intersect(rownames(posList[[1]]), colnames(mat)),]
posList[[2]] <- posList[[2]][intersect(rownames(posList[[2]]), colnames(mat)),]
posList[[3]] <- posList[[3]][intersect(rownames(posList[[3]]), colnames(mat)),]

# Plot
par(mfrow=c(1,3), mar=rep(2,4))
plotEmbedding(posList[[1]], groups=section, main='section 1', cex=2)
plotEmbedding(posList[[2]], groups=section, main='section 2', cex=2)
plotEmbedding(posList[[3]], groups=section, main='section 3', cex=2)

write.csv(posList,"breast_cancer_pos_data.csv",quote=FALSE,row.names=TRUE)
