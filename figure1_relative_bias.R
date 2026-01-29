rm(list=ls())
library(ggplot2)
library(ggpubr)

N <- c(500,1000,5000,10000,50000)
runs <- 1000
ols.res <- matrix(,nrow=runs,ncol=length(N))
tsls.res <- matrix(,nrow=runs,ncol=length(N))
beta <- 1
rho <- 0.5

for (N0 in 1:length(N)){
  
fname = paste('relativebias_',N[N0],'.RData', sep="")
load(fname)
# collect simulation results
ols.res[,N0] <- ols
tsls.res[,N0] <- tsls
}

hist.dat.ols <- list(); hist.plot.ols <- list(); est.plot.ols <- list(); hist.dat.tsls <- list(); hist.plot.tsls <- list(); est.plot.tsls <- list()
for (i in 1:length(N)){
hist.dat.ols[[i]]<-as.data.frame(data.frame(method = factor(rep(c("OLS"), each=runs)), std.est = c(ols.res[,i])))
hist.dat.tsls[[i]]<-as.data.frame(data.frame(method = factor(rep(c("2SLS"), each=runs)), std.est = c(tsls.res[,i])))
hist.plot.ols[[i]]<-ggplot(hist.dat.ols[[i]], aes(x = std.est))
hist.plot.tsls[[i]]<-ggplot(hist.dat.tsls[[i]], aes(x = std.est))
minq <- min(c(quantile(ols.res,0.001)[[1]],quantile(tsls.res,0.001)[[1]]))
maxq <- min(c(quantile(ols.res,0.999)[[1]],quantile(tsls.res,0.999)[[1]]))
bin <- seq(minq,maxq,(maxq-minq)/28)
est.plot.ols[[i]] <- hist.plot.ols[[i]] + geom_density(aes(color = method, fill = method) , linewidth = 0,  alpha = 1) + coord_cartesian(xlim = c(minq, maxq)) + theme_bw() + scale_color_manual(values=rev(c('forestgreen'))) + scale_fill_manual(values=rev(c('forestgreen'))) + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="right", legend.text = element_text(size=12)) + xlab("OLS estimates") + ggtitle(paste0("n = ",N[i])) + geom_vline(xintercept = 1, colour = "firebrick1") + theme(plot.title = element_text(hjust = 0.5, size=11), plot.subtitle = element_text(hjust = 0.5), legend.key.size=unit(1, 'cm'), legend.text = element_text(size=12), legend.title = element_text(size=13), axis.title=element_text(size=12)) + theme(plot.title=element_text(size=13))
est.plot.tsls[[i]] <- hist.plot.tsls[[i]] + geom_density(aes(color = method, fill = method) , linewidth = 0,  alpha = 1) + coord_cartesian(xlim = c(minq, maxq)) + theme_bw() + scale_color_manual(values=rev(c('gold1'))) + scale_fill_manual(values=rev(c('gold1'))) + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="right", legend.text = element_text(size=12)) + xlab("2SLS estimates") + ggtitle(paste0("n = ",N[i])) + geom_vline(xintercept = 1, colour = "firebrick1") + theme(plot.title = element_text(hjust = 0.5, size=11), plot.subtitle = element_text(hjust = 0.5), legend.key.size=unit(1, 'cm'), legend.text = element_text(size=12), legend.title = element_text(size=13), axis.title=element_text(size=12)) + theme(plot.title=element_text(size=13))
}

ols.plots <- ggarrange(est.plot.ols[[1]],est.plot.ols[[2]],est.plot.ols[[3]],est.plot.ols[[4]],est.plot.ols[[5]], common.legend = TRUE, legend = "right", ncol = 5, nrow = 1)
tsls.plots <- ggarrange(est.plot.tsls[[1]],est.plot.tsls[[2]],est.plot.tsls[[3]],est.plot.tsls[[4]],est.plot.tsls[[5]], common.legend = TRUE, legend = "right", ncol = 5, nrow = 1)

estimation.plot <- ggarrange(ols.plots,tsls.plots,ncol=1,nrow=2)

ggsave(filename="figure1.png", plot=estimation.plot, device="png",width=12.5,height=4.25,dpi=300, bg="white")
