Sys.setenv(LANG="en")
library(spacexr)
library(Matrix)

x  <- readRDS("myRCTDde_updated2.rds")

#slotNames(x)
#[1] "spatialRNA"         "originalSpatialRNA" "reference"
#[4] "config"             "cell_type_info"     "internal_vars"
#[7] "results"            "de_results"         "internal_vars_de"

#> slotNames(x@spatialRNA)
#[1] "coords" "counts" "nUMI"

#> dim(x@spatialRNA@counts)
#[1]  5160 44091

#> dim(x@spatialRNA@coords)
#[1] 44091     2


#E1 = as.data.frame(as.matrix(x@spatialRNA@counts))
E2 = as.data.frame(as.matrix(x@spatialRNA@coords))

gname <- rownames(x@spatialRNA@counts)
cname <- colnames(x@spatialRNA@counts)


writeMM(x@spatialRNA@counts,'matrix.mtx')
write.csv(x@results$results_df, 'rctd_doublet.csv',quote=FALSE,row.names=TRUE)
write.csv(E2,"coords.csv",quote=FALSE,row.names=TRUE)
write.csv(gname,"gname.csv",quote=FALSE,row.names=TRUE)
write.csv(cname,"cname.csv",quote=FALSE,row.names=TRUE)

#write.csv(E1,"counts.csv",quote=FALSE,row.names=TRUE)
