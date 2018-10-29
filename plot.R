rm(list=ls())
library(compiler)
library(deSolve)
library(pracma)
#install.packages("sensitivity")
#library(sensitivity)

setwd("/Users/adil/Desktop/untitled folder")
source("solver.R")

c1=9e1

U = solver()
si = U[a4:b4,]
X = seq(0,L-h, by = h)


# c2=1
# c3=1e2
# c4=1e-8

C1 = seq(9e1/10, 9e1*7, length.out = 10)
n = length(C1)
mae = matrix(NA, nrow = 5, ncol = n)

for (i in 1:n){
  
  c1 = C1[i]
  U_iter = solver()
  
  si_iter = U_iter[a4:b4,]
  
  mae[1,i] = mean(si_iter[,1000] - si[,1000])
  mae[2,i] = mean(si_iter[,2000] - si[,2000])
  mae[3,i] = mean(si_iter[,3000] - si[,3000])
  mae[4,i] = mean(si_iter[,4000] - si[,4000])
  mae[5,i] = mean(si_iter[,5000] - si[,5000])
  
}

plot_mmp <-function(U){
  m = U[a1:b1,]
  mp = U[a2:b2,]
  
  plot(X, m[,5000]+mp[,5000], main = "Total hypal density", ylab = "M + Mp", xlab = "Distance from center (cm)", 
       lwd=2, type="l", col="blue", ylim = c(0,max( m[,1000:5000] + mp[,1000:5000])))
  points(X, m[,1000]+mp[,1000], lwd=2, type="l", col="orange")
  points(X, m[,2000]+mp[,2000], lwd=2, type="l", col="red")
  points(X, m[,3000]+mp[,3000], lwd=2, type="l", col="black")
  points(X, m[,4000]+mp[,4000], lwd=2, type="l", col="yellow")
}


plot_m <-function(U){
  
  m = U[a1:b1,]
  
  plot(X, m[,1000], main = "Active hypal density", lwd=2, ylab = "M", xlab = "Distance from center (cm)",
       type="l", col="orange", ylim=c(0,max(m[,1000:5000])))
  points(X, m[,2000], lwd=2, type="l", col="red")
  points(X, m[,3000], lwd=2, type="l", col="black")
  points(X, m[,4000], lwd=2, type="l", col="yellow")
  points(X, m[,5000], lwd=2, type="l", col="blue")
  
}


plot_si <- function(U){
  
  si = U[a4:b4,]
  
  plot(X, si[,1000], main = "Internal substrate concentration", lwd=2, ylab = "Si", xlab = "Distance from center (cm)",
       type="l", col="orange", ylim=c(0,max(si[,1000:5000])))
  points(X, si[,2000], lwd=2, type="l", col="red")
  points(X, si[,3000], lwd=2, type="l", col="black")
  points(X, si[,4000], lwd=2, type="l", col="yellow")
  points(X, si[,5000], lwd=2, type="l", col="blue")

}

plot_se <- function(U){
  
  se = U[a5:b5,]
  
  plot(X, se[,1000], main = "External substrate concentration", lwd=2,ylab = "Si", xlab = "Distance from center (cm)",
       type="l", col="orange", ylim=c(0,max(se[,1000:5000])))
  points(X, se[,2000], lwd=2, type="l", col="red")
  points(X, se[,3000], lwd=2, type="l", col="black")
  points(X, se[,4000], lwd=2, type="l", col="yellow")
  points(X, se[,5000], lwd=2, type="l", col="blue")
}

 plot_p <-function(U){
   
  p = U[a3:b3,]
   

  plot(X, p[,1000], main = "Tip density", lwd=2, type="l", ylab = "P", xlab = "Distance from center (cm)",
       col="orange", ylim=c(0,max(p[,1000:5000])))
  points(X, p[,2000], lwd=2, type="l", col="red")
  points(X, p[,3000], lwd=2, type="l", col="black")
  points(X, p[,4000], lwd=2, type="l", col="yellow")
  points(X, p[,5000], lwd=2, type="l", col="blue")
}
 
plot_mp <- function(U){
  
  mp = U[a2:b2,]
  
  plot(X, mp[,1000], lwd=2, type="l", col="orange", ylim=c(0,max(mp[,1000:5000])))
  points(X, mp[,2000], lwd=2, type="l", col="red")
  points(X, mp[,3000], lwd=2, type="l", col="black")
  points(X, mp[,4000], lwd=2, type="l", col="yellow")
  points(X, mp[,5000], lwd=2, type="l", col="blue")
}
