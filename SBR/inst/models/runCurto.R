library(SBR) 
(curto=readSBMLR(file.path(system.file(package="SBR"), "models/curto.r")) )
summary(curto)
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=sim(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))
# out is a dataframe with class set to deSolve to overload plot, see ?plot.deSolve
