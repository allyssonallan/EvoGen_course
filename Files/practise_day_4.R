
## This script shows an example of how to detecting natural positive selection in SNP data

# Data set is comprised of ~ 200,000 SNPs from the metaboChip in HapMap populations: YRI (Africa), CEU (Europe), CHB (East Asian)

# Data is in the widely-used plink format (http://pngu.mgh.harvard.edu/~purcell/plink/), specifically in allele frequency per cluster format

# if you want to have a look at the input file:
# less -S ../Data/hapmap.frq.strat

# ----

# one of the most powerful approaches to detect selection is to compare genetic variation between 2 or more populations, assuming that the control populations are neutrally evolving for that specific gene

# FST is a very commonly used metric to identify changes in allele (or haplotype) frequencies between populations

# we can compute per-site FST values using a method-of-moments estimator

Rscript ../Scripts/plink2fst.R

# this generates a file Results/hapmap.fst
# it may take a while so you can copy it from ../Data

# less -S Results/hapmap.fst


# Manhattan plots
Rscript ../Scripts/plotManFST.R Results/hapmap.fst Plots/hapmap.fst.pdf

# calculate PBS
fst=read.table("Results/hapmap.fst", stringsAsFact=F, header=T)

# YRI.CEU + CEU.CHB - YRI.CHB
pbs= ( (-log(1-fst$FST.YRI.CEU)) + (-log(1-fst$FST.CEU.CHB)) - (-log(1-fst$FST.YRI.CHB)) ) / 2
pbs[which(pbs<0)]=0

# plot
cols=rep("grey", nrow(fst)); cols[which( (fst$chrom %% 2) == 1)]="lightgrey"

plot(x=fst$cpos, y=pbs, col=cols, frame=F, xlab="", xaxt="n", ylab="PBS", main="CEU", pch=16)

# check top hits
fst[which.max(pbs),]
# https://genome-euro.ucsc.edu
# http://haplotter.uchicago.edu/
# http://hgdp.uchicago.edu/cgi-bin/gbrowse/HGDP/

# second hit
fst[which(pbs>1.4 & pbs<2),]

# check which gene
# rs482000
# SLC35F3

# compute p-value
source("../Scripts/functions.R")

ms.command <- "Software/msdir/ms 326 100000 -s 1 -I 3 118 120 88 -n 1 1.68 -n 2 1.12 -n 3 1.12 -eg 0 2 72 -eg 0 3 96 -ma x 2.42 1.52 2.42 x 7.73 1.52 7.73 x -ej 0.029 3 2 -en 0.029 2 0.29 -ej 0.19 2 1 -en 0.30 1 1 | gzip > Results/ms.txt.gz"

system(ms.command, intern=F)

sim.chroms=readMs("Results/ms.txt.gz" , 326)$hap

nreps=length(sim.chroms)

sim.pbs<-rep(NA, nreps)
for (i in 1:nreps) {

	pops=list()
        pops[[1]]=sim.chroms[[i]][1:118] # YRI
        pops[[2]]=sim.chroms[[i]][119:238] # CEU
        pops[[3]]=sim.chroms[[i]][239:326] # CHB

        sim.fst=chroms2fst(pops)

        sim.pbs[i]=( (-log(1-sim.fst[1])) + (-log(1-sim.fst[2])) - (-log(1-sim.fst[3])) ) / 2

        if ((i %% 100)==0) cat(i, "\t", length(which(sim.pbs>=1.42))/i, "\n")

}

sim.pbs[which(sim.pbs<0)]=0

hist(sim.pbs, main="Simulations under neutrality", xlab="PBS", breaks=20)
abline(v=1.42, lty=2)





