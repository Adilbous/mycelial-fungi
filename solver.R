rm(list=ls())
library(compiler)
library(deSolve)
library(pracma)

# Import des paramètres d'intêret 
#source("/Users/adil/Desktop/OBT/Mod Vivant/R/Projet 25/parameters.R")

# Import des fonctions (dérivées etc) et des paramètres d'intêret

#setwd(getwd())
#source("function.R")
source("/Users/adil/Desktop/OBT/Projet 25 Scripts/function.R")

# Fonction des dérivées stiff

fs <-function(u){ 
  
  du_s=rep(NA,5*NS)
  
  mj  = u[a1:b1]
  mpj = u[a2:b2]
  pj  = u[a3:b3]
  sij = u[a4:b4]
  sej = u[a5:b5]
  
  du_s[a1:b1]  = dmdt_s(mj, mpj, pj, sij, sej) 
  du_s[a2:b2] = dmpdt_s(mj, mpj, pj, sij, sej)
  du_s[a3:b3] = dpdt_s(mj, mpj, pj, sij, sej)
  du_s[a4:b4] = dsidt_s(mj, mpj, pj, sij, sej)
  du_s[a5:b5] = dsedt_s(mj, mpj, pj, sij, sej)
  
  return(du_s)
}

fns <-function(u){ 
  
  du_ns=rep(NA,5*NS)
  
  mj  = u[a1:b1]
  mpj = u[a2:b2]
  pj  = u[a3:b3]
  sij = u[a4:b4]
  sej = u[a5:b5]
  
  du_ns[a1:b1]  = dmdt_ns(mj, mpj, pj, sij, sej) 
  du_ns[a2:b2] = dmpdt_ns(mj, mpj, pj, sij, sej)
  du_ns[a3:b3] = dpdt_ns(mj, mpj, pj, sij, sej)
  du_ns[a4:b4] = dsidt_ns(mj, mpj, pj, sij, sej)
  du_ns[a5:b5] = dsedt_ns(mj, mpj, pj, sij, sej)
  
  return(du_ns)
}

# Initialisation de U

init_U <-function(){
  
  U = matrix(NA,5*NS,NT)
  
  # on utilise U, le vecteur (mj m'j pj sij sej)
  # on a alors pour les index :
  # de   a1 = 1      à  b1 = NS : mj
  # de   a2 = NS+1   à  b2 = 2*NS : mpj
  # de   a3 = 2*NS+1 à  b3 = 3*NS : pj
  # de   a4 = 3*Ns+1 à  b4 = 4*NS : sij
  # de   a5 = 4*NS+1 à  b5 = 5*NS : sej
  
  U[a1:b1,1]  = c(rep(m0,lambda%/%h),rep(0,NS-lambda%/%h))
  U[a2:b2,1] = rep(0,NS)
  U[a3:b3,1] = c(rep(p0,lambda%/%h),rep(0,NS-lambda%/%h))
  U[a4:b4,1] = c(rep(si0,lambda%/%h),rep(0,NS-lambda%/%h))
  U[a5:b5,1] = rep(se0,NS)
  
  # for (i in a1:5){
  #   U[i,1] = (5-i)*m0/5
  # }
  # U[6:b1,1] = 0
  #  
  #  for (i in a3:(a3+5)){
  #    U[i,1] = (5-i)*p0/5
  #  }
  #  U[(a3+6):b3,1] = 0
  # 
  #  for (i in a4:(a4+5)){
  #    U[i,1] = (5-i)*si0/5
  #  }
  #  U[(a4+6):b4,1] = 0
  #  
  #  U[a2:b2,1] = rep(0,NS)
  #  U[a5:b5,1] = rep(se0,NS)
  
  return(U)
}

# Execution du schéma numérique

# Fonction de calcul d'erreur

rho <-function(Utau, Udemitau){
  
  Sum = 0
  for (i in 1:5*NS){
    Sum = Sum + ((Udemitau[i] - Utau[i]) / ( tol + tol*abs(Utau[i])))**2
  }
  return( 1/3 * sqrt( 1/(5*NS) * Sum ) )
  
}

# Solver pour chaque pas de temps i 

solver_step <-function(U, tau, i){
  #z1 : Runge Kutta d'ordre 2 pour la premi?re partie de la d?composition
  
  k1 = fns(U[,i-1])
  k2 = fns(U[,i-1]+tau/2*1/2*k1)
  k3 = fns(U[,i-1]+tau/2*(1/2*k1+1/2*k2))
  z1 = U[,i-1]+tau/2*(1/3*k1+1/3*k2+1/3*k3)
  
  #z2 : calcul avec le Jacobien
  # Essayer rampe (condition initiale), si c'est mieux => Schéma Upwind pour l'espace
  # Sinon, Euler implicite 
  
  v0 = z1
  v1 = v0+tau/2*fs(v0)
  A = diag(5*NS)-tau/2*jacobian(fs,v1)
  # print(kappa(A, exact = TRUE))
  
  v2 = v1+tau/2*solve(A)%*%fs(v1)
  z2 = v2
  
  #z3 : Runge Kutta d'ordre 2 pour la premi?re partie de la d?composition
  # Crank Nicholson pour les 2 ? 
  
  k1 = fns(z2)
  k2 = fns(z2+tau/2*1/2*k1)
  k3 = fns(z2+tau/2*(1/2*k1+1/2*k2))
  z3 = z2+tau/2*(1/3*k1+1/3*k2+1/3*k3)
  
  return(z3)
}

solver <-function(){
  
  U = init_U()
  tau = to
  
  for (i in 2:NT){
    print(i)
    print(tau)
    
    Utau = solver_step(U, tau, i)
    
    Uinterm = U 
    Uinterm[,i] = solver_step(U, tau/2, i)
    
    Udemitau = solver_step( Uinterm, tau/2, i)
    
    print(rho(Utau, Udemitau))
    
    if ( rho(Utau, Udemitau) < 1 ){
      U[,i] = Udemitau
    }
    
    else {
      
      Uinterm1 = U
      Uinterm1[,i] = solver_step(Uinterm1, tau/3, i)
      
      Uinterm2 = Uinterm1
      Uinterm2[,i] = solver_step(Uinterm2, tau/3, i)
      
      U[,i] = solver_step( Uinterm2, tau/3, i)
    }
    
    tau = tau * min(2, max( (4/5) * rho(Utau, Udemitau)**(-1/3)  , 1/4 ) )
    
    #U[,i]=z3
  }
  
  return(U)
}