"readSBMLR"<-function(filename)
{# SBMLR file input, SBML model object output
  
  
  source(filename,local=TRUE) # loads in model as a list of lists
  
  ## initialize names for indexing by names/IDs of reactions, species and rules
#nSpecies=length(model$species)
#if (nSpecies>0){
#sIDs=NULL;  # initialize before assignments
#for (i in 1:nSpecies) sIDs[i]<-model$species[[i]]$id
#names(model$species)<-sIDs
#}
# ********* Commented block above is much simpler with simplified lapply = sapply
  if (length(model$species)>0) 
    names(model$species)<-sapply(model$species,function(x) x$id)
  
  
# # not clear why I wanted this next block  ... could probably get rid of it
# # OK, see readSBML, it was free coming in from xmlTreeParse so I kept it, 
# # but I must have forgotten that this meant I would need to create it coming in from R source. 
#   notes=model$notes
#   nNotes=length(notes)
#   if (nNotes>0){
#     con <- xmlOutputDOM()
#     con$addTag("notes", close=FALSE)
#     con$addTag("body",attrs=c(xmlns="http://www.w3.org/1999/xhtml"), close=FALSE  )
#     for (i in 1:nNotes)  
#     {
#       con$addTag("p",close=FALSE  )
#       con$addNode(xmlTextNode(notes[i]))
#       con$closeTag()
#     }
#     con$closeTag();con$closeTag()
#     dom=con$value()[[1]]
#     model$htmlNotes=dom
#   }
#   
#   
#   
  
# Make rules in MathML, expressions and functions 
  nRules=length(model$rules)
  if (nRules>0){
    ruleIDs=NULL
    for (i in 1:nRules)
    {
      #cat(i,"     ",class(model$rules[[i]]$strLaw),"\n")
      model$rules[[i]]$exprLaw=parse(text=model$rules[[i]]$strLaw)[[1]]
      r=model$rules[[i]]$inputs
      e=model$rules[[i]]$exprLaw
#print(r)
      model$rules[[i]]$law=makeLaw(r,p=NULL,e)
#print("here")
      # model$rules[[i]]$mathmlLaw=R2MathML(e)$value()[[1]]
      ruleIDs[i]<-model$rules[[i]]$idOutput
    }
    names(model$rules)<-ruleIDs
  }
  
# Make reactions in MathML, expressions and functions 
  nReactions=length(model$reactions)
  if (nReactions>0){
    rIDs=NULL
    for (i in 1:nReactions)
    {
      if (is.null(model$reactions[[i]]$reversible)) model$reactions[[i]]$reversible=FALSE 
#cat(i,"     ",class(model$reactions[[i]]$strLaw),"\n")
      model$reactions[[i]]$exprLaw=parse(text=model$reactions[[i]]$strLaw)[[1]]
      r=model$reactions[[i]]$reactants
      m=model$reactions[[i]]$modifiers
#      r=c(r,m)
      p=names(model$reactions[[i]]$parameters)
      e=model$reactions[[i]]$exprLaw
      model$reactions[[i]]$law=makeLaw(c(r,m),p,e)
#      model$reactions[[i]]$law=makeLaw(r,p,e)
      # model$reactions[[i]]$mathmlLaw=R2MathML(e)$value()[[1]]
      rIDs[i]<-model$reactions[[i]]$id
      
    }
    names(model$reactions)<-rIDs
  }
  
  
  class(model)<-"SBMLR"
  return(model)
}

# 
# R2MathML <-function(e)
# {  # takes R expressions into mathML XMLNode
#   a <- xmlOutputDOM()
#   
#   r2ml<- function(type)   # maps R operator symbols into MathML 
#     switch(type,
#            "*" = "times",
#            "/" = "divide",
#            "+" = "plus",
#            "-" = "minus",
#            "^" = "power",
#            "exp" = "exp",
#            "log" = "ln",
#            "not found")    # end of r2ml sub-function definition. 
#   
#   
#   recurs<-function(e)  
#   { if(e[[1]]=="(") {e=e[[2]]} # remove parentheses  
#     a$addTag("apply", close=FALSE  )
#     a$addTag(r2ml(as.character(e[[1]])))
#     for (j in 2:length(e))
#       if ((class(e[[j]])=="name")|(class(e[[j]])=="numeric"))  {
#         op<-switch(class(e[[j]]),name="ci",numeric="cn")
#         a$addTag(op,as.character(e[[j]]))
#       } else  Recall(e[[j]]) 
#     
#     a$closeTag()
#   }
#   recurs(e)
#   return(a)
# }

# NOTE: because I decided to keep the mathML when reading SBML, to make it easier to save later, 
# I thought I should make the mathML upfront from the SBMLR representation, so that the model in R
# is the same regardless of the source being SBML or SBMLR

# the following is called by both readSBML and readSBMLR so it outside where both can reach it.
# Note that keeing it here instead of in a separate file => no need to document it
"makeLaw"<-function(r,p,e){
  # takes reactant list r, parameter list p and rate law R expression e 
  # and makes a reaction rate law function out of them.
  lawTempl=function(r,p=NULL){ }
  i=2
  for (j in seq(along=p)){
    #		if(!is.null(p))
    #		for (j in 1:length(p)){
    body(lawTempl)[[i]]<-call("=",as.name(p[j]),call("[",as.name("p"),p[j]))
    i=i+1}
  #   for (j in 1:length(r)){ 
  for (j in seq(along=r)){
    body(lawTempl)[[i]]<-call("=",as.name(r[j]),call("[",as.name("r"),r[j]))
    i=i+1}
  body(lawTempl)[[i]]<-e
  lawTempl
}

