"sim" <-
    function(model,times, modulator=NULL,X0=NULL, ...)  # this is a wrapper for lsoda
{
  mi=summary(model) # mi = model info
# print(mi);next
#   my.atol <- rep(1e-4,mi$nStates)
  attach(model$globalParameters)  # in Morrison this makes Keq globally available 
  mod=0
  if (class(modulator)=="numeric") mod=1
  if (class(modulator)=="list") mod=2
#print(modulator)
  
  fderiv<-function(t, X, p) { # state derivative function sent to ODEsolve
    v=rep(0,mi$nReactions)
    xp=rep(0,mi$nStates)
    St=mi$S0
    St[mi$BC==FALSE]=X

# TR version of this    
    if (mi$nRules>0) 
      for (j in 1:mi$nRules)
        St[model$rules[[j]]$idOutput]=model$rules[[j]]$law(St[model$rule[[j]]$inputs]) 
    
    if (p["mod"]==1) {
      if (t<0) 
        m=rep(1,length(modulator)) else 
        m=modulator
      names(m)<-mi$rIDs
    }
    
    for (j in 1:mi$nReactions) {
        mrj=model$reactions[[j]]
        rm=c(mrj$reactants,mrj$modifiers)
        P=mrj$parameters
        if (p["mod"]==0)   v[j]=mrj$law(St[rm],P)
        if (p["mod"]==1)   v[j]=m[mi$rIDs[j]]*mrj$law(St[rm],P)
        if (p["mod"]==2)   v[j]=modulator[[mi$rIDs[j]]](t)*mrj$law(St[rm],P)
      }
    xp=mi$incid%*%v
    names(xp)<-names(mi$y0)
    names(v)<-mi$rIDs
    aux=c(v,St[mi$BC==TRUE])
    list(xp,aux)
  }    
  # ******************  END fderiv function definition
  
  if (is.null(X0)) X0=mi$y0
    out=ode(y=X0,times=times,fderiv,  parms=c(mod=mod), ...)
  detach(model$globalParameters)
  out
} 
