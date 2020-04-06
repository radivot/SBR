# This script models the remdesivir (RDV) PK data on the last page (page 15) of 
# https://www.who.int/ebola/drc-2018/summaries-of-evidence-experimental-therapeutics.pdf?ua=1
#rdv plasma
auc_1=c(1840,3261,1255)  #h*ng/ml for 75, 150, 75 mg over 2, 2, and 0.5 hrs
tp5_1=c(0.84,1.1,1)  # t1/2 in hours
#rdv-TP (GS-443902) in PBM cells
auc_2=c(176,295,394)
tp5_2=c(43,36,49)
1/(1-exp(-0.693*24/tp5_2))  #recreats ratio=c(3.1,2.7,3.46), check

library(mrgsolve)# setwd("~/ccf/hobbs/methods/mrgsolve")
options(mrgsolve.soloc="~/ccf/hobbs/methods/mrgsolve/C")
mod1=mread("pk1",modlib())#modlib()=>"/Library/Frameworks/R.framework/Versions/3.6/Resources/library/mrgsolve/models"
str(mod1)
example("param")
param(mod1)
V=20 # assuming V=20L
(CL=0.693*V/tp5_1[1])  # 0.693 × Vd /tp5
mod1=param(mod1,list(CL=16.5,KA=100)) #fast absorption mimics injection 
mev=mod1 %>% ev(amt = 75, ii = 24, addl = 0) 
mev%>%mrgsim(end = 15, delta = 0.1)%>%plot(CP~time)
d=mev%>%mrgsim_df(end = 15, delta = 0.1)
head(d)
auc1=0.1*sum(d$CP) # in h*mg/L  
auc1*1000 # to get to h*ng/ml, 4396.77 vs 1840 measured is a bit off (lost faster)
#smaller step size only increases it to 4542, so not helping. Need extra compartment sucking it off


mod2=mread("pk2",modlib())#modlib()=>"/Library/Frameworks/R.framework/Versions/3.6/Resources/library/mrgsolve/models"
param(mod2)
V=20 # assuming V=20L
(CL=0.693*V/tp5_1[1])  # 0.693 × Vd /tp5
mod2=param(mod2,list(CL=16.5,KA=100)) #fast absorption mimics injection 
mev=mod2 %>% ev(amt = 75, ii = 24, addl = 0) 
mev%>%mrgsim(end = 15, delta = 0.1)%>%plot(CP~time)
d=mev%>%mrgsim_df(end = 15, delta = 0.1)
head(d)
(auc1=0.1*sum(d$CP)) # in h*mg/L  
auc1*1000 # to get to h*ng/ml, 4342 is smaller (Q=2)

mod2=param(mod2,list(CL=16.5,KA=100,Q=0.01)) #fast absorption mimics injection 
mod2=param(mod2,list(CL=16.5,KA=100,Q=10)) #fast absorption mimics injection 
mod2=param(mod2,list(CL=16.5,KA=100,Q=1)) #fast absorption mimics injection 
mev=mod2 %>% ev(amt = 75, ii = 24, addl = 0) 
mev%>%mrgsim(end = 15, delta = 0.1)%>%plot(CP~time)
d=mev%>%mrgsim_df(end = 15, delta = 0.1)
head(d)
(auc1=0.1*sum(d$CP)) # in h*mg/L  
auc1*1000 # Q=0.01 => 4394 (worse), Q=10 4393, Q-1 4317 
#need some permanently lost, either stored in cells or metabolized (or excreted by other routes).


