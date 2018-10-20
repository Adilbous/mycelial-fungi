rm(list=ls())
library(compiler)
library(deSolve)
library(pracma)

# Paramètres du système mis à l'échelle

nu=10
d=1e-1
b=1e3
f=1e4
c1=9e1
c2=1
c3=1e2
c4=1e-8
Di=1e1
Da=1e1
De=1e-1
lambda=0.2

m0=0.1
p0=0.1
si0=0.4
se0=0.3

# Paramètres de la résolution

# Discrétisation spatiale

L=2 #Etendue
h=0.1 #grain
epsilon=1e-30
NS=L%/%h 

# Discrétisation temporelle

Ttot=4 #Etendue
tau=1/24/10 #grain
NT=9

a1 = 1
b1 = NS
a2=1*NS+1
b2=2*NS
a3=2*NS+1
b3=3*NS
a4=3*NS+1
b4=4*NS
a5=4*NS+1
b5=5*NS


# Fonction de limitation de Van Leer

Phi = function(r){
  return (r+abs(r))/(1+abs(r))
}

# Dérivées temporelles stiff des variables d'état 

dmdt_s <-function(mj, mpj, pj, sij, sej){
  dmjs=rep(NA,NS)
  for(j in 1:NS){
    dmjs[j]=nu*sij[j]*pj[j]-d*mj[j]
  }
  return(dmjs)
}

dmpdt_s <-function(mj, mpj, pj, sij, sej){
  dmpjs=rep(NA,NS)
  for(j in 1:NS){
    dmpjs[j]=d*mj[j]
  }
  return(dmpjs)
}

dpdt_s <-function(mj, mpj, pj, sij, sej){
  dpjs=rep(NA,NS)
  for(j in 1:NS){
    dpjs[j]=b*sij[j]*mj[j]-f*mj[j]*pj[j]
  }
  return(dpjs)
}

dsidt_s <-function(mj, mpj, pj, sij, sej){
  dsijs=rep(NA,NS)
  for(j in 1:NS){
    if(j==1){
      pxj=(-3*pj[1]+4*pj[2]-pj[3])/(2*h)
    }
    else if(j==NS){
      pxj=(3*pj[NS]-4*pj[NS-1]+pj[NS-2])/(2*h)
    }
    else{
      pxj=(pj[j+1]-pj[j-1])/(2*h)
    }
    dsijs[j]=c1*sij[j]*mj[j]*sej[j]-c2*nu*sij[j]*pj[j]-c4*Da*mj[j]*sij[j]*abs(pxj)
  }
  
  return(dsijs)
}

dsedt_s <-function(mj, mpj, pj, sij, sej){
  dsejs=rep(NA,NS)
  for(j in 1:NS){
    if(j==1){
      A=De*(sej[j+1]-2*sej[j]+sej[j+1])/(h*h)
    }
    else if(j==NS){
      A=De*(sej[j-1]-2*sej[j]+sej[j-1])/(h*h)
    }
    else{
      A=De*(sej[j+1]-2*sej[j]+sej[j-1])/(h*h)
    }
    dsejs[j]=A-c3*sij[j]*sej[j]*mj[j]
  }
  
  return(dsejs)
}

# Dérivées temporelles non-stiff des variables d'état 

dmdt_ns <-function(mj, mpj, pj, sij, sej){
  dmjns=rep(NA,NS)
  for(j in 1:NS){
    dmjns[j]=0
  }
  
  return(dmjns)
}

dmpdt_ns <-function(mj, mpj, pj, sij, sej){
  dmpjns=rep(NA,NS)
  for(j in 1:NS){
    dmpjns[j]=0
  }
  return(dmpjns)
}

dpdt_ns <-function(mj, mpj, pj, sij, sej){
  dpjns=rep(NA,NS)
  for(j in 1:NS){
    fp=function(k){
      return(nu*sij[k]*pj[k])
    }
    rp=function(k){
      if(k==1){ #fp-1=0
        return((fp(k+1)-fp(k)+epsilon)/(fp(k)+epsilon))
      }
      else if(k==NS){ #fpNS+1=0
        return((-fp(k)+epsilon)/(fp(k)-fp(k-1)+epsilon))
      }
      else{
        return((fp(k+1)-fp(k)+epsilon)/(fp(k)-fp(k-1)+epsilon))
      }
    }
    if(j==1){
      dpjns[j]=-1/h*(fp(j)+1/2*Phi(rp(j))*(fp(j)))
    }
    else if(j==2){
      dpjns[j]=-1/h*(fp(j)+1/2*Phi(rp(j))*(fp(j)-fp(j-1))-
                      (fp(j-1)+1/2*Phi(rp(j-1))*(fp(j-1))))
    }
    else{
      dpjns[j]=-1/h*(fp(j)+1/2*Phi(rp(j))*(fp(j)-fp(j-1))-
                      (fp(j-1)+1/2*Phi(rp(j-1))*(fp(j-1)-fp(j-2))))
    }
  }
  return(dpjns)
}

dsidt_ns <-function(mj, mpj, pj, sij, sej){
  dsijns=rep(NA,NS)
  #on d?finit quelques fonctions interm?diaires pour all?ger par la suite
  px=function(k){
    if(k==1){
      pxk=(-3*pj[1]+4*pj[2]-pj[3])/(2*h)
    }
    else if(k==NS){
      pxk=(3*pj[NS]-4*pj[NS-1]+pj[NS-2])/(2*h)
    }
    else{
      pxk=(pj[k+1]-pj[k-1])/(2*h)
    }
    return (pxk)
  }
  w=function(k){
    return(Da*mj[k]*px(k))
  }
  fa=function(k){
    return(w(k)*sij[k])
  }
  ra=function(k){
    if(k==1){ #fa-1=0
      return((fa(k+1)-fa(k)+epsilon)/(fa(k)+epsilon))
    }
    else if(k==NS){ #faNS+1=0
      return((-fa(k)+epsilon)/(fa(k)-fa(k-1)+epsilon))
    }
    else{
      return((fa(k+1)-fa(k)+epsilon)/(fa(k)-fa(k-1)+epsilon))
    }
  }
  Fd=function(k){
    if(k==0){ 
      return(0)
    }
    else if(k==NS){
      return(0)
    }
    else{
      return(Di*(mj[k+1]-mj[k])/2*(sij[k+1]-sij[k])/h)
    }
  }
  Fa = function(k){
    #print(w(k))
    if( w(k)>=0 | is.na(w(k)) ) {
        if(k==0){ 
          return(0)
        }
        else if(k==1){
          return(fa(k)+1/2*Phi(ra(k)*(fa(k))))
        }
        else{
          return(fa(k)+1/2*Phi(ra(k)*(fa(k)-fa(k-1))))
        }        
    }
    else{
      if(k==NS-1){ 
        return(fa(k+1)+1/2*Phi(1/ra(k+1)*(fa(k+1))))
      }
      else if(k==NS){
        return(0)
      }
      else{
        return(fa(k+1)+1/2*Phi(1/ra(k+1)*(fa(k+1)-fa(k+2))))
      }
    }
     
  }
  #----fin des fonctions interm?diaires
  dsijns[1]=1/h*(Fd(1))-1/h*(Fa(1))
  for(j in 2:NS){
    dsijns[j]= 1/h * (Fd(j) - Fd(j-1)) -1/h * (Fa(j) - Fa(j-1))
  }
  
  return(dsijns)
  
}

dsedt_ns <-function(mj, mpj, pj, sij, sej){
  dsejns=rep(NA,NS)
  for(j in 1:NS){
    dsejns[j]=0
  }
  return(dsejns)
}
  



