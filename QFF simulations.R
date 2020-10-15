nruns=10000
setwd("F:/Fruitfly stuff/Model YV Sensitivity/Baseline_small")
#### model parameters
male.max.age = 42
luredist = 405		# metres
mmdd = 600  ## mean male daily distance

### notes on distance units... 
## the units on the map are 'pixels', which are 15m by 15m
## units in alltrees and the design are also in pixels, so multiply by 15m to get metres
## simulation is run on pixels, so units in metres are divided by 10 in the code

# class      : RasterLayer 
# dimensions : 3853, 2974, 11458822  (nrow, ncol, ncell)
# resolution : 15, 15  (x, y)
# extent     : 347229, 391839, 5794668, 5852463  (xmin, xmax, ymin, ymax)
# crs        : +proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
# source     : F:/Fruitfly stuff/Model YV Sensitivity/Baseline_small/Risk Map/Intro_risk.tif 
# names      : Intro_risk 
# values     : 0.05873247, 1  (min, max)

# This is the information of the raster layer that needs to be converted to a matrix

col_names = seq(347229, 391839-1, 15)
row_names = seq(5794668, 5852463-1, 15)

#################  surveillance design... 


design=read.csv('design_grid.csv')
design$x = abs((design$x-347229) /15)
design$y = abs((design$y-5794668) /15)
tinds = read.csv('Trees.csv')
tinds$size = ceiling(tinds$InclProb*10)

tsize = tinds$size

alltrees = tinds[,c(3:6)]
names(alltrees) = c('x','y','id','size')
alltrees$id = seq(1, nrow(alltrees),1)

alltrees$x = abs((alltrees$x-347229) /15)
alltrees$y = abs((alltrees$y-5794668)/15)

load('allcloseones')
head(allcloseones)

load('allpleave')  #allpleave
head(allpleave)

load('introrisk')
introrisk = mat_risk_t


##########################################################################
### Run these lines to get a plot of the female dispersal kernel. 
distwtfunc = function(dist){
	2^(-dist/200)
}
plot(1:1000,distwtfunc(1:1000),col='red',t='l',lwd=2,xlab="Distance (m)",ylab="Weighted probability")
##########################################################################


introrisk2=round(introrisk,5)
unvals=unique(as.vector(introrisk2))
forplotrisk = introrisk2
forplotrisk[introrisk2==unvals[1]]=1
forplotrisk[introrisk2==unvals[4]]=0.75
forplotrisk[introrisk2==unvals[2]]=0.5
forplotrisk[introrisk2==unvals[3]]=0.4
forplotrisk[introrisk2==unvals[5]]=0

dev.new()

image(x = 0:2974, y = 0:3853, forplotrisk, useRaster = TRUE)
alph = 0.2*tsize
alph[alph>1]=1
points(alltrees,cex=0.1,col=rgb(0,0,0,alph),pch=16)

points(design$x,design$y,col='purple',pch=3, cex = 3)
ntraps = dim(design)[1]




##############################################


allres=NULL
for (run in 1:nruns){
# dev.new()
print(run)
#image(x=1:1000,y=1:1000,introrisk)
#initincurpt = sample(length(introrisk),1000,prob=introrisk)
#introrisk[initincurpt]=2
#for (i in 1:100){
# image(x = 0:4648, y = 0:6200, introrisk, useRaster = TRUE)  
initincurpt = sample(length(introrisk),1,prob=introrisk)
rr = initincurpt%%nrow(introrisk) 
rr[rr==0]=nrow(introrisk)
cc = (initincurpt-1+nrow(introrisk))%/%nrow(introrisk)
#introrisk[cbind(rr,cc)]=10
dists = with(alltrees,sqrt( (x-rr)^2 + (y-cc)^2 )*15)	##multiply by 10 to covert to metres
wts = sapply(dists,distwtfunc )
firsttree = alltrees[sample(length(wts),1,prob=wts),]
# points(firsttree$x,firsttree$y,col='blue',pch=16,cex=0.5)
#}
# dev.off()

malelist = data.frame()
infectedlist = data.frame(id=firsttree$id,day=1)
for (t in 2:600){
  # dev.new()
	print(t)
	szs=alltrees$size[infectedlist$id]
	days=t-infectedlist$day
	nleaversp = unlist(sapply(1:length(szs) , function(ii) allpleave[[szs[ii]]][days[ii]]))
	nleaversp[is.na(nleaversp)]=0
	### males
	nleaversm = rpois(length(nleaversp), nleaversp)
	prodids = infectedlist$id[nleaversm>0]
	nleaversm = nleaversm[nleaversm>0]
	ntreesprod = length(prodids)
	if (ntreesprod>0){
		newmales = data.frame(x=rep(alltrees$x[prodids],nleaversm) , y=rep(alltrees$y[prodids],nleaversm), 
			age=0) 
		malelist = rbind(malelist , newmales )
		}
	malelist$age=malelist$age+1
	malelist=subset(malelist,age<=male.max.age)
	nmales=dim(malelist)[1]
	mdist = rnorm(nmales)*mmdd/0.7977/15
	a = runif(nmales)*pi*2
	malelist$x=malelist$x+mdist*cos(a)
	malelist$y=malelist$y+mdist*sin(a)
	### females
	nleavers = rpois(length(nleaversp), nleaversp)
	prodids = infectedlist$id[nleavers>0]
	#producers = alltrees[,prodids]
	nleavers = nleavers[nleavers>0]
	ntreesprod = length(prodids)
	if (ntreesprod>0) for (ii in 1:ntreesprod ){
		wts=sapply( allcloseones[[prodids[ii]]]$dist , distwtfunc )	## this distance already in metres
		nexttreesids = allcloseones[[prodids[ii]]]$id[sample(length(wts),nleavers[ii],replace=TRUE,prob=wts)]
		#points(nexttrees$x,nexttrees$y,col='green',pch=16,cex=0.9)
		nexttreesids = setdiff(nexttreesids,infectedlist$id )
		if (length(nexttreesids)>0) infectedlist = rbind( infectedlist , data.frame(id=nexttreesids ,day=t))
	}
	infectedtrees = alltrees[infectedlist$id,]
	
	
	### plot
	# 0:4648, y = 0:6200
	# plot(malelist$x,malelist$y,col='darkgreen',pch=16,cex=0.05,xlim=c(0,4648),ylim=c(0,6200))
	# points(malelist$x,malelist$y,col='darkgreen',pch=16,cex=0.05)
	# points(infectedtrees$x,infectedtrees$y,col='green',pch=16,cex=0.3)
	### plot to file
	# if (TRUE & t%%10==0){
	#   # dev.new()
	#   png(paste('outfig',run,+t+1000,'.png'))
	#   image(x = 0:4648, y = 0:6200,forplotrisk, useRaster = TRUE)
	#   alph = 0.2*tsize
	#   alph[alph>1]=1
	#   points(tinds,cex=0.1,col=rgb(0,0,0,alph),pch=16)
	#   points(design$x,design$y,col='purple',pch=3)
	#   points(malelist$x,malelist$y,col='darkgreen',pch=16,cex=0.05)
	#   points(infectedtrees$x,infectedtrees$y,col='green',pch=16,cex=0.3)
	#   dev.off()
	# }
	
	
	###### surveillance
	if (!is.null(design) & nmales>0){
		dzs=sapply(1:ntraps, function(i) (design$x[i]-malelist$x)^2 + (design$y[i]-malelist$y)^2 )
		if (min(dzs) < (luredist/15)^2) {		## dividing by 10 converts to pixels
			dzs=matrix(dzs,ncol=ntraps)
			trapi=which(dzs==min(dzs),arr.ind=TRUE)[2]
			print(paste("run",run,"DETECTED at time:",t))
			#points(design$x[trapi],design$y[trapi],col='cyan',pch=16)
			allres=rbind(allres,c(t,trapi,nmales,dim(infectedlist)[1],rr, cc))
			
			png(paste('outfig',run,+t+1000,'.png'))
			image(x = 0:2974, y = 0:3853,forplotrisk, useRaster = TRUE)
			points(tinds,cex=0.1,col=rgb(0,0,0,alph),pch=16)
			points(design$x,design$y,col='purple',pch=3)
			points(malelist$x,malelist$y,col='darkgreen',pch=16,cex=0.05)
			points(infectedtrees$x,infectedtrees$y,col='green',pch=16,cex=0.3)
			dev.off()

			break
		}

		
	}
	
	else if(t == 600){
	  trapi='NA'
	  allres=rbind(allres,c(t,trapi,nmales,dim(infectedlist)[1],rr, cc))
	  print(paste("No Detection, nmales:", nmales))
	  
	}
	
	# dev.off()
}}

allres = allres[1:1000,]
allres = data.frame(allres)
allres = data.frame(apply(allres, 2, as.numeric))
allres_sub = allres[allres$time != 600,]

base = base [base $trap.id != 'NA', ]
names(allres)=c('time','trap.id','nmales','ntrees', 'x', 'y')


write.csv(allres, "results_baseline_small_kmeans_traps.csv")
str(allres)




