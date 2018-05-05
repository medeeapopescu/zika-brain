# Performs differential expression between pre-defined populations

.libPaths(c("/datastore/wanxinw/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))
library(scde)

# NEED USER DEFINITION
# Load necessary data
root_dir<-"/local10G/wanxinw/scde_diff_exp/ERA_run1.2.3.4/"
counts <- read.csv(paste0(root_dir,"counts_epi_250.csv"),row.names=1,check.names=F,stringsAsFactors=F)
cluster.id <- as.matrix(read.csv(paste0(root_dir,"groups_epi_250.csv"),row.names=1,check.names=F,stringsAsFactors=F))


for (i in 1:2){
    for ( j in (i+1):3){
    
    # define two clusters of interest
    cells<-rownames(cluster.id)[union(which(cluster.id[,1]==i),which(cluster.id[,1]==j))]
    print(length(cells))
    # slice grping factor and counts accordingly
    grping<-cluster.id[,1][union(which(cluster.id[,1]==i),which(cluster.id[,1]==j))]
    print(length(grping))
    counts_curr<-counts[,which(colnames(counts)%in%cells)]
    print(dim(counts_curr))
    
    # output name
    name_o.ifm<-paste0(root_dir, "epi/", paste0("epi.o.ifm.",i,j,".csv"))
    name_ediff<-paste0(root_dir, "epi/", paste0("epi.ediff.clust.",i,j,".ordered.csv"))
    print(name_o.ifm)
    print(name_ediff)

    # generate "current group" object associating cell names with grp id
    grping.curr <- factor(grping,levels=c(i,j))
    
    # calculate o.ifm and o.prior
    o.ifm <- scde.error.models(count=counts_curr, groups=grping.curr, n.cores=8, threshold.segmentation=T, save.model.plots=F, verbose=1)
    write.csv(o.ifm,name_o.ifm)
    o.prior <- scde.expression.prior(models=o.ifm, counts=counts_curr, length.out=400, show.plot=F)
    
    # diff exp analysis and order
    ediff.temp <- scde.expression.difference(o.ifm, counts=counts_curr, o.prior, groups  =  grping.curr, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
    ediff.ordered <- ediff.temp[order(ediff.temp$Z, decreasing = TRUE),]
    
    print("Calculations done.")
    # output
    write.csv(ediff.ordered, file=name_ediff)
    print("Output done.")
    }
}


