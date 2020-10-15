setwd("F:/Fruitfly stuff/Model YV Sensitivity/Baseline_small")


alltrees = read.csv('Trees.csv')
alltrees = alltrees[,2:5]

#### record the 400 closest for each tree
allcloseones=list()


for (i in 1:nrow(alltrees)){
  print(i)
  dists = with(alltrees,sqrt( (X-X[i])^2 + (Y-Y[i])^2 )*15)
  dists[i] = 9999999999
  closeones = data.frame(id=order(dists)[1:400],dist=dists[order(dists)][1:400])
  allcloseones[[i]] = closeones 
}

save(allcloseones,file='allcloseones')
load('allcloseones')
##check###
i=4900
i=900
aco=allcloseones[[i]]
plot(alltrees[aco$id,]$X, alltrees[aco$id,]$Y,cex=0.01)
text(alltrees[aco$id,]$X, alltrees[aco$id,]$Y,1:400,cex=0.6)
text(alltrees[aco$id,]$X[1:200], alltrees[aco$id,]$Y[1:200],1:200,cex=0.6)
points(alltrees[i,]$X, alltrees[i,]$Y,col='red', cex =1,pch = 16)
points(alltrees[aco$id,]$X[1:3], alltrees[aco$id,]$Y[1:3],col='green')
































