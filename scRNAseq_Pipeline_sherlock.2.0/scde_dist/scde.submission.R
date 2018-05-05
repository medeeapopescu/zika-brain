# Specify path to look for scde and load libraries
.libPaths(c("/datastore/wanxinw/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))
library(scde)
library(parallel)
library(boot)

# Get name for input count file (has to be .csv)
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    stopifnot(length(args)>0)
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

print(input) #.csv file name

# declare output file names
o.ifm_output<-gsub(".csv","_o.ifm.csv",input)
mode.fail_output<-gsub(".csv","_mode.fail.csv",input)

print(o.ifm_output)
print(mode.fail_output)

# import count files
counts <- read.csv(input,row.names=1,stringsAsFactors=F,check.names=F)
o.ifm <- scde.error.models(counts=counts,n.cores=12,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1)
write.csv(o.ifm,o.ifm_output)
# o.ifm <- read.csv(o.ifm_output, row.names=1,stringsAsFactors=F, check.names=F)

print(dim(o.ifm))
valid.cells <- o.ifm$corr.a>0

print(length(valid.cells))
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = T)
o.fpm <- scde.expression.magnitude(o.ifm,counts=counts)

o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
par(mfrow=c(1,1),mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1);
plot(c(),c(),xlim=range(o.prior$x),ylim=c(0,1),xlab="expression magnitude (log10)",ylab="drop-out probability")
invisible(apply(o.fail.curves,2,function(y) lines(x=o.prior$x,y=y,col="orange")))

p.self.fail <- scde.failure.probability(models=o.ifm,counts=counts)
cell.names <- colnames(counts); names(cell.names) <- cell.names;

# reclculate posteriors with the individual posterior modes 
jp <- scde.posteriors(models=o.ifm,counts,o.prior,return.individual.posterior.modes=T,n.cores=8)
# find joint posterior modes for each gene - a measure of MLE of group-average expression
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes here)
mat <- log10(exp(jp$modes)+1);
# weighted distance
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
  unlist(lapply(cell.names,function(nam2) {
    corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
  }))
},mc.cores=8)),upper=F);

mode.fail.matrix <- as.matrix(mode.fail.dist)
write.csv(mode.fail.matrix,mode.fail_output)