### Load Data 
#setwd("~/Documents/classes/1st year/4th term/Methods IV/Final") 
setwd("~/GitHub/methods_project/")
data <- read.delim("challenge1_data.txt")

z<-c(0.5,0.5,-0.5,-0.5)
k<-kronecker(diag(10),t(z))
y<-data[,6:ncol(data)]
y<-as.matrix(y[,as.vector(rbind(matrix(1:20,2,10),matrix(21:40,2,10)))])

y<-y^0.25

diffs <-y %*% t(k)

cov = t(diffs) %*% diffs/nrow(data)
invcov = solve(cov)

quad.form <- function(row,sigma.inv) {
  qform = t(row) %*% sigma.inv %*% row
  return(qform)
}

f.stats <- apply(diffs,1,quad.form,sigma.inv=invcov)

data = cbind(data,f.stats)

#Plot F statistics, merge with data
#ggplot of F statistic by chromosome with facet grid

require(ggplot2)
require(dplyr)
require(reshape2)
require(longitudinal)

x = data

factnames = seq(1,22,1)
factnames = c(factnames,"X")

for (i in 1:23) {
  factnames[i] = paste("chr",factnames[i],sep="")
}

x$chr= factor(x$chr, levels = factnames)

pdf(file = paste("f.statistic_by_location",".pdf", sep=""),
    width=18, height=6, onefile=T)
print(
  ggplot(x)+
    theme_bw() +
    theme(axis.text.x=element_blank())+
    theme(
      plot.background = element_blank()
      #,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
    ) +
    theme(axis.line = element_line(color = 'black'))+
    theme(legend.position="none") + 
    labs(x = "Locus",
      y = "F statistic")+
    geom_point(aes(x=start, y=f.stats, color=chr, size=f.stats), alpha=0.25) + 
    facet_grid(.~chr, scales="free_x", space="free_x")+
    ggtitle("F statistic for each locus, sorted by chromosome")
)
dev.off()