#Coding the ABW decomposition model. Allison SD, Bradford MA & Wallenstein MD. 2010. Nature Geoscience. Soil-carbon response to warming dependent on microbial physiology. 3: 336-340.

#initial parameter values, based on spinup in Supplementary Table 2
#time-note I am dropping the interval from original ABW
time<- 200000 #end model run after this many model hours. 
temp<-20 #temperature, degrees C.
t<-0 #start at time 0

#initial pool sizes
SOC<- 100 #mg/cm3
DOC<- 0.5 #mg/cm3
MIC<- 0.5 #mg/cm3
Enz<- 0.01 #mg/cm3

#inputfluxes
inputSOC<- 0.0005 #mg/cm3/h
inputDOC<- 0.0005 #mg/cm3/h

#constant turnover rates and parameters
r.death<- 0.0002 # biomass turnover rate (1/h)
r.enz.prod<-0.000005 #fraction allocated to enzyme production (1/h)
r.enz.loss<-0.001 #fraction enzyme pool that turns over per hour (1/h)
MICtoSOC<- 0.5 #fraction of microbial turnover that enters SOC pool, rather than DOC

#temperature sensitive parameter values
CUE.0<- 0.63 #intercept CUE
CUE.slope<- -0.016 #change in CUE per degree temperature
Vmax.0<- 100000000 #intercept Vmax term. mg SOM / cm3 / h
Vmax.uptake.0 <- 100000000 #intercept uptake Vmax term
Km.0 <- 500 #intercept of temperature sensitive Km value
Km.uptake.0 <- 0.1 #Intercept of temperature sensitive Km.uptake
Km.slope<- 5 #change in Km value per degree C
Km.uptake.slope <- 0.01 #change in Km.uptake parameter per degree C
Ea<- 47 #kj/mol Activation energy for SOC degrading enzymes
Ea.uptake<- 47 #kj/mol Activation energy for uptake transporters 
gas.const <- 0.008314 #the universal gas constant, bro. 

#calculating temperature sensitive parameters
Vmax.uptake<- Vmax.uptake.0 * exp(-Ea.uptake/(gas.const * (temp + 273)))
Km.uptake<- Km.uptake.slope*temp + Km.uptake.0
CUE<- CUE.slope*temp + CUE.0
Vmax<- Vmax.0 * exp(-Ea/(gas.const*(temp+273)))
Km<- Km.slope*temp + Km.0

#create output matrix for saving outputs
out<-matrix(rep(0,time*6),nrow=time,dimnames=list(NULL,c('hours','SOC','DOC','mbc','Enz','CO2')))

for(i in 1:time){
  #fluxes!
  ASSIM = Vmax.uptake * MIC * (DOC / ( Km.uptake + DOC)) #microbial uptake
  DEATH = r.death * MIC #microbial death
  EPROD = r.enz.prod * MIC #microbial enzyme production
  ELOSS = r.enz.loss * Enz #microbial enzyme turnover
  DECOMP = Vmax * Enz * (SOC / (Km + SOC))
  #DECOMP flux not alowed to be greater than total SOC pool size
  DECOMP<- ifelse(DECOMP>SOC,SOC,DECOMP)

  #changes in pools per unit time!
  dMICdt<- ASSIM*CUE - DEATH - EPROD #change in microbial biomass per unit time
  dEnzdt<- EPROD-ELOSS #change in enzyme pool size per unit time
  dSOCdt<- inputSOC + DEATH*MICtoSOC - DECOMP #change in SOC pool per unit time
  dDOCdt<- inputDOC + DEATH*(1-MICtoSOC) + DECOMP + ELOSS - ASSIM #change in DOC pool per unit time.
  
  #CO2 respiration
  CO2<- ASSIM*(1-CUE)
  
  #update pools
  MIC<- MIC + dMICdt
  Enz<- Enz + dEnzdt
  SOC<- SOC + dSOCdt
  DOC<- DOC + dDOCdt
  
  #update output matrix
  out[i,]<-c(i,SOC,DOC,MIC,Enz,CO2)
#end model function
}
out<-data.frame(out)

plot(out$Enz)

