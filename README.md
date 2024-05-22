install.packages("devtools")

options(timeout = 600000000) ### set this to avoid timeout error

devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

conda install conda-forge::r-rcppziggurat

```
load("predictive_variable_slideseq_cerebellum.RData")
c<-load("predictive_variable_slideseq_cerebellum.RData")
#[1] "anterior_rep1" "anterior_rep2" "anterior_rep3" "nodular_rep1"
#[5] "nodular_rep2"  "nodular_rep3"
write.csv(anterior_rep1,"a1.csv",quote=FALSE,row.names=TRUE)
```
