
fin=commandArgs(T)[1]

source("Software/treemix-1.12/src/plotting_funcs.R");
pdf(paste(fin,".pdf",sep="",collapse=""))
par(mfrow=c(1,2))
plot_tree(fin)
plot_resid(fin, "Data/poporder.txt")
dev.off()




