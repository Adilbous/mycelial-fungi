rm(list=ls())
library(compiler)
library(deSolve)
library(pracma)

#setwd(getwd())
#source("solver v1.R")
source("/Users/adil/Desktop/OBT/Projet 25 Scripts/solver v1.R")

U = solver()
m = U[a1:b1,]
mp = U[a2:b2,]
p = U[a3:b3,]
si = U[a4:b4,]
se = U[a5:b5,]


X = seq(0,L-h, by = h)

plot(X, m[,1]+mp[,1], lwd=2, type="p", col="green",)
plot(X, m[,2]+mp[,2], lwd=2, type="p", col="orange")
plot(X, m[,3]+mp[,3], lwd=2, type="p", col="red")
plot(X, m[,4]+mp[,4], lwd=2, type="p", col="blue")
plot(X, m[,5]+mp[,5], lwd=2, type="p", col="blue")
plot(X, m[,6]+mp[,6], lwd=2, type="p", col="blue")
plot(X, m[,7]+mp[,7], lwd=2, type="p", col="blue")
plot(X, m[,8]+mp[,8], lwd=2, type="p", col="blue")
plot(X, m[,9]+mp[,9], lwd=2, type="p", col="blue")

# plot(X, se[,1], lwd=2, type="p", col="green",)
# plot(X, se[,2], lwd=2, type="p", col="orange")
# plot(X, se[,3], lwd=2, type="p", col="red")
# plot(X, se[,4], lwd=2, type="p", col="blue")
# plot(X, se[,5], lwd=2, type="p", col="blue")
# plot(X, se[,6], lwd=2, type="p", col="blue")
# plot(X, se[,7], lwd=2, type="p", col="blue")
# plot(X, se[,8], lwd=2, type="p", col="blue")
# plot(X, se[,9], lwd=2, type="p", col="blue")

plot(X, p[,1], lwd=2, type="p", col="green",)
plot(X, p[,2], lwd=2, type="p", col="orange")
plot(X, p[,3], lwd=2, type="p", col="red")
plot(X, p[,4], lwd=2, type="p", col="blue")
plot(X, p[,5], lwd=2, type="p", col="blue")
plot(X, p[,6], lwd=2, type="p", col="blue")
plot(X, p[,7], lwd=2, type="p", col="blue")
plot(X, p[,8], lwd=2, type="p", col="blue")
plot(X, p[,9], lwd=2, type="p", col="blue")