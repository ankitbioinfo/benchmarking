
```
install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
conda install conda-forge::r-rcppziggurat
```
to quickly check that 

```
load("predictive_variable_slideseq_cerebellum.RData")
c<-load("predictive_variable_slideseq_cerebellum.RData")
#[1] "anterior_rep1" "anterior_rep2" "anterior_rep3" "nodular_rep1"
#[5] "nodular_rep2"  "nodular_rep3"
write.csv(anterior_rep1,"a1.csv",quote=FALSE,row.names=TRUE)
```

find reindex columns 
```
df3=pd.read_csv(path+'coords.csv',index_col=0)
coord=df3.to_numpy()
df4= pd.read_csv(path+'rctd_doublet.csv',index_col=0)
newdf= df4.reindex(df3.index)
```
