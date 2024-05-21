library(MERINGUE)
Sys.setenv(LANG="en")

data("drosophila")
pos <- drosophila$pos
gexp <- drosophila$mat

dim(gexp)

par(mfrow=c(1,3), mar=rep(2,4))
g <- 'apt'
plotEmbedding(pos[,c(1,2)], col=gexp[g,], main='X-Y')
plotEmbedding(pos[,c(3,2)], col=gexp[g,], main='Z-Y')
plotEmbedding(pos[,c(1,3)], col=gexp[g,], main='X-Z')

write.csv(pos,"drosophila_pos_data.csv",quote=FALSE,row.names=TRUE)
write.csv(gexp,"drosophila_counts_data.csv",quote=FALSE,row.names=TRUE)
