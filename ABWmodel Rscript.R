#Coding the ABW decomposition model. Allison SD, Bradford MA, Wallenstein MD. 2010. Nature Geoscience.Soil-carbon response to warming dependent on microbial physiology. 3: 336-340.

#initial parameter values.
#time
endTime<- 26280000 #end model run after this many model hours. 
interval<- 5000
temp<-20
t<-0 #start at time 0

#initial pool sizes
SOC<- 111.8765 #mg/cm3
DOC<- 0.00144928 #mg/cm3
MIC<- 2.19159 #mg/cm3
Enz<- 0.0109579 #mg/cm3

#inputfluxes
inputSOC<- 0.0005 #mg/cm3/h
inputDOC<- 0.0005 #mg/cm3/h

#constant turnover rates and parameters
r.death<- 0.0002 # biomass turnoevr rate 1/h
r.enz.prod<-0.000005 #fraction allocated to enzyme production (1/h)
r.enz.loss<-0.001 #fraction enzyme pool that turns over per hour (1/h)
MICtoSOC<- 0.5 #fraction of microbial turnover that enters SOC pool, rather than DOC

#temperature sensitive parameter values
CUE.0<- 0.63 #intercept CUE
CUE.slope<- -0.016 #change in CUE per degree temperature
Vmax.0<- 100000000 #intercept Vmax term. mg SOM / cm3 / h
Vmax.uptake.0 <- 100000000 #intercept uptake Vmax term
Km.0 <- 500 #intercept Km value
Km.uptake.0 <- 0.1 #Intercept Km of uptake value
Km.slope<- 5
Km.uptake.slope <- 0.01
Ea<- 47 #kj/mol
Ea.uptake<- 47
gas.const <- 0.008314

#create empty vectors to record model outputs
output.t<- c()
output.MIC<- c()
output.SOC<- c()
output.CO2<- c()

while(t<endTime){
  #calculating temperature sensitive parameters
  Vmax.uptake<- Vmax.uptake.0 * exp(-Ea.uptake/(gas.const * (temp + 273)))
  Km.uptake<- Km.uptake.slope*temp + Km.uptake.0
  CUE<- CUE.slope*temp + CUE.0
  Vmax<- Vmax.0 * exp(-Ea/(gas.const*(temp+273)))
  Km<- Km.slope*temp + Km.0
  
  #fluxes!
  ASSIM = Vmax.uptake * MIC * (DOC / ( Km.uptake + DOC)) #microbial uptake
  DEATH = r.death * MIC #microbial death
  EPROD = r.enz.prod*Enz #microbial enzyme production
  ELOSS = r.enz.loss*Enz #microbial enzyme turnover
  DECOMP = Vmax*Enz*(SOC/(Km + SOC))
  #DECOMP flux not alowed to be greater than total SOC pool size
  DECOMP<- ifelse(DECOMP>SOC,SOC,DECOMP)

  #changes in pools per unit time!
  dMICdt<- ASSIM*CUE - DEATH - EPROD #change in microbial biomass per unit time
  dEnzdt<- EPROD-ELOSS #change in enzyme pool size per unit time
  dSOCdt<- inputSOC + DEATH*MICtoSOC - DECOMP
  dDOCdt<- inputDOC + DEATH*(1-MICtoSOC) + DECOMP + ELOSS - ASSIM
  
  #CO2 respiration
  CO2<- ASSIM*(1-CUE)
  
  #update pools
  MIC<- MIC + dMICdt
  Enz<- Enz + dEnzdt
  SOC<- SOC + dSOCdt
  DOC<- DOC + dDOCdt
  
  #advance time step
  t=t+interval
  
  #save model outputs every timestep
  output.t<- c(output.t,t)
  output.MIC<- c(output.MIC,MIC)
  output.SOC<- c(output.SOC,SOC)
  output.CO2<- c(output.CO2,CO2)
}

#problems- getting negative CO2 fluxes. SOC pool increases linearly with time. Biomass oscillates, but thats to be expected.