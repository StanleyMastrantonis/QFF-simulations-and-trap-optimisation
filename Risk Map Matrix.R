setwd("F:/Fruitfly stuff/Model YV Sensitivity/Baseline_small/Risk Map")
library(raster)
intro = raster("Intro_risk.tif")
plot(intro)

mat_risk = raster::as.matrix(intro)

col_names = seq(347229, 391839-1, 15)
row_names = seq(5794668, 5852463-1, 15)
mat_risk = matrix(mat_risk,3853, 2974)
dimnames(mat_risk) = list(as.numeric(row_names), as.numeric(col_names))

mat_risk_df = which(mat_risk >=0, arr.ind = T)
image(mat_risk, useRaster = TRUE)

rotate = function(x) t(apply(x, 2, rev))
mat_risk_t = rotate(mat_risk)
mat_risk_t[is.na(mat_risk_t)] = 0 
image(mat_risk_t, useRaster = TRUE)

save(mat_risk_t,file = 'introrisk')
