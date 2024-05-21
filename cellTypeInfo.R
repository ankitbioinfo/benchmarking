
Sys.setenv(LANG="en")


load("joint_clean_com_conservative.RData")

#"com.main"           "com.final.all"      "com.malig.final"    "com.nonmalig.final"
#cat("nrow = B1 ", nrow(ct), "ncol = ", ncol(ct), "\n");

write.csv(com.main,"JF_com_main.csv",quote=FALSE,row.names=TRUE)
write.csv(com.final.all,"JF_comFinalAll.csv",quote=FALSE,row.names=TRUE)
write.csv(com.malig.final,"JF_comMaligFinal.csv",quote=FALSE,row.names=TRUE)
write.csv(com.nonmalig.final,"JF_comNonMaligFinal.csv",quote=FALSE,row.names=TRUE)
