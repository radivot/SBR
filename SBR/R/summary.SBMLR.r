"summary.SBMLR"<-function(object, ...)
{
	model=object
	sIDs=names(model$species)
	rIDs=names(model$reactions)
	ruleIDs=names(model$rules)
	nReactions=length(model$reactions);nSpecies=length(model$species);nRules=length(model$rules) 
	
# Species
	S0=NULL;BC=NULL # initialize 
	for (i in 1:nSpecies){
		BC[i]=model$species[[i]]$bc; 
		S0[i]=model$species[[i]]$ic
	}
	names(S0)<-sIDs 
	names(BC)<-sIDs 
	y0=S0[BC==FALSE]
	nStates=length(y0)

# these 4 lines replace the block commented after it  
  rLaws=sapply(model$reactions,function(x) x$strLaw)
  globals=model$globalParameters
	attach(globals)  # e.g. for global coordination of k5 in SOD2012
	V0=sapply(model$reactions,function(x) x$law(S0[c(x$reactants,x$modifiers)],x$parameters))
	names(V0)<-rIDs
	detach(globals)  
  
# Incidence Matrix
	incid=matrix(rep(0,nStates*nReactions),nrow=nStates)
	indx=(1:nSpecies)[BC==FALSE]
	for (i in 1:nStates)
		for (j in 1:nReactions)
		{if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$products)) 
				incid[i,j] = summary(factor(model$reactions[[j]]$products))[[model$species[[indx[i]]]$id]]
			if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$reactants)) 
				incid[i,j] = incid[i,j]-summary(factor(model$reactions[[j]]$reactants))[[model$species[[indx[i]]]$id]]  }     
	rownames(incid)<-names(y0)
	
# return a list of model information
	options(stringsAsFactors=FALSE)
	DFs=data.frame(index=1:nSpecies,initialConcentrations=S0,boundaryConditions=BC);row.names(DFs)<-sIDs
	DFr=data.frame(index=1:nReactions,Laws=rLaws,initialFluxes=V0);   row.names(DFr)<-rIDs
	list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,nStates=nStates,
			y0=y0,nReactions=nReactions,rIDs=rIDs,rLaws=rLaws, V0=V0,globalVec=unlist(globals), # P0, VP,
			incid=incid,nRules=nRules,ruleIDs=ruleIDs,species=DFs,reactions=DFr) 
}



