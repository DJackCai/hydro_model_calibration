#### ==== This is the R version of the matlab codes for the 
#### Modello Idrologico Semi-distributio in Continuo (MISDc) model in Brocca et al. (2012)

### Original codes and Data can be found HERE:
## https://github.com/IRPIhydrology/MISDc

#### Main script 
rm(list = ls())
graph.theme.beta<-
  theme(axis.title=element_text(size=14),axis.text = element_text(size=11),plot.title = element_text(hjust=0.5,size=14),legend.text = element_text(size=12),legend.title = element_text(size=14))

library(lubridate)
library(pracma)
library(ggplot2)

## Load the data 
data.mat<-read.table("migi_0205_daily.txt")
date<-data.mat[,1]
#### Important: Change the datenum format in matlab to UTC date format in R
date.Rformat<-as.POSIXct((date - 719529)*86400, origin = "1970-01-01", tz = "UTC")

data.mat.new<-data.mat
data.mat.new[,1]<-date.Rformat
datamat<-data.mat.new
X_INIT_cali<-c(0.4292,242.4420,10,1.4417,1.5316,1.7985,4.4634,2.6975,15.4304,7.4725)

MISDc.simu<-function(datamat,PARS,Ab,name=NA) {
# datamat:  a) date (in numeric Matlab date numbers format,need to translate to R data) 
# b) rainfall depth (in mm)
# c) air temperature (°C) # d) observed discharge data (m^3/s)
 M<-dim(datamat)[1] # total number of time steps 
 date<-datamat[,1]
 RAIN<-datamat[,2]; TEMP<-datamat[,3]; Qobs<-datamat[,4]
 #### define the time step
 delta_T<-round(mean(diff(date,lag = 1),na.rm=T) *24 *10000)/10000
 delta_T<-as.numeric(delta_T)
 Month<-month(date)
 
 #### PARS: model parameters input values
 Frac.sz.init<-PARS[1]  # initial fraction of surface wetness
 Wmax.rz<-PARS[2]  # total water capacity in root-zone
 m1<-PARS[3]  # exponent of drainage of surface layer
 KS1<-PARS[4]  # hydraulic conductivity of the surface layer 
 gamma1<-PARS[5]  # lag time coefficient in determining the IUH
 Kc<-PARS[6]  # parameter of potential evapotranspiration
 alpha<-PARS[7]  # exponent of runoff
 Cm<-PARS[8] #snow module parameter degree day (Bayesian GLUE in Barry)
 m2<-PARS[9] # exponent of drainge for root-zone
 KS2<-PARS[10] # hydraulic conductivity in root-zone
 Wmax.sz<-150   # maximum water capacity in surface layer treated as fixed
 
 dt=0.2 # computation time step in hour 
 # hydraulic conductivity: mm/h --> mm/delta_T
 KS1<-KS1*(delta_T)
 KS2<-KS2*(delta_T)
 
 # Snow module, with snow module function. Return a time series
 Rain_out<-snow_model(as.matrix(RAIN),as.matrix(TEMP),-0.5,0.5,Cm)[[1]];
 SWE.out<-snow_model(as.matrix(RAIN),as.matrix(TEMP),-0.5,0.5,Cm)[[2]];
 # Potential Evapotranspiration parameter over year: here treated as fixed 
 # L=c(0.21,0,22,0.23,0.28,0.3,0.31,0.3,0.29,0.27,0.25,0.22,0.2)
 L=rep(0.25,12)
 Ka<-1.26
 ET<-(TEMP>0)*(Kc*(Ka*L[Month]*(0.46*TEMP+8)-2))/(24/delta_T) # vector of ET

 # initialisation for storgae
 # WW_sz & WW_rz: relative wetness; BF: baseflow; RO: surface runoff
 BF<-RO<-WW_sz<-WW_rz<-Perc_all<-DeepPrec_all<-rep(0,M)  
 SZ.all<-RZ.all<-rep(0,M)  
 ### main routine 
 # Starting conditions of soil wetness index  
 W.sz<-Frac.sz.init*Wmax.sz
 W.rz<-Frac.sz.init*Wmax.rz  
 S<-NaN; Pcum<-IE<-0
 
 ## loop through the time step
 for (t in 1:M) {
   # Direct runoff on first day (A general soil wetness-effective runoff relation)
   IE<-Rain_out[t]*((W.sz/Wmax.sz) ^ alpha)
   Evap<-ET[t]*(W.sz/Wmax.sz)  # only ET from surface?
   if (W.rz<Wmax.rz) {Perc=KS1* ((W.sz/Wmax.sz) ^m1)  }
   else {Perc=0}
   
   # deep percolation
   DeepPrec<-KS2*((W.rz/Wmax.rz) ^m2)
   # Additional constraint on the percolation
   Perc[Perc>(0.6*Wmax.sz)]=0.6*Wmax.sz
   Perc_all[t]<-Perc;  DeepPrec_all[t]<-DeepPrec
   # DeepPrec[DeepPrec>0.6*Wmax.rz] =0.6*Wmax.rz
   
   #### Soil water balance accounting for both layers
   W.sz<-max(0,W.sz+Rain_out[t]-IE-Perc-Evap+SWE.out[t])  # correct
   W.rz<-max(0,W.rz+Perc-DeepPrec)
 
   # saturation excess
   if (W.sz>=Wmax.sz) {Excess.sz<-W.sz-Wmax.sz;W.sz=Wmax.sz }
   else {Excess.sz<-0}
   if (W.rz>=Wmax.rz) {Excess.rz<-W.rz-Wmax.rz; W.rz=Wmax.rz }
   else {Excess.rz<-0}
   if (W.sz<0) {W.sz<-0}
   if (W.rz<0) {W.rz<-0}
   SZ.all[t]<-W.sz
   RZ.all[t]<-W.rz  
   # store the results for soil wetness index in time step t 
   WW_sz[t]<-W.sz/Wmax.sz
   WW_rz[t]<-W.rz/Wmax.rz
   # baseflow 
   BF[t]<-KS2*( ((W.sz+W.rz)/(Wmax.sz+Wmax.rz))^(m2)) # based on drainge component of root-zone
   # Calculate the total runoff contribution 
   RO[t]<-IE+Excess.sz+Excess.rz  # surface runoff + saturation excess
 }

 # Convolution over the Instantaneous Unit Hydrograph
  
  IUH1<-IUH_comp(gamma1,Ab,dt,delta_T)*dt;
  IUH1<-IUH1/sum(IUH1)
  # assume 1 storage only
  IUH2<-IUH_NASH(1,0.5*gamma1,Ab,dt,delta_T)*dt; 
  IUH2<-IUH2/sum(IUH2)
  # interpolate the time series of the input runoff 
  RO_int<-(interp1(1:M,RO,seq(1,M,dt)))
  BF_int<-(interp1(1:M,BF,seq(1,M,dt)))

  ## Direct runoff hydrograph to determine the volume 
  temp1<-conv(IUH1,RO_int)  # temporary state 
  temp2<-conv(IUH2,BF_int)
 # Qsim1<-temp2[seq(1,M*round(1/dt),round(1/dt))]*(Ab*1000/delta_T/3600);  # baseflow component 
  
 Qsim<-( temp1[seq(1,M*round(1/dt),round(1/dt))] + temp2[seq(1,M*round(1/dt),round(1/dt))]) *(Ab*1000/delta_T/3600)
  # calculate streamflow prediction performance 
  RMSE<-( mean((Qsim-Qobs)**2,na.rm=T) ) ^0.5
  NS<-1-(sum((Qsim-Qobs)**2,na.rm=T)/sum((Qobs-mean(Qobs,na.rm = T))**2,na.rm=T))
  # ANSE<-1-sum ( ((Qobs+mean(Qobs,na.rm = T))*(Qsim-Qobs))**2,na.rm = T )/
   # sum ( ((Qobs+mean(Qobs,na.rm = T))*(Qobs-mean(Qobs,na.rm = T)))**2,na.rm = T )
  KGE<-klingupta(Qsim,Qobs)
  df.stream<-data.frame(date,Qobs,Qsim)
  assign("streamdf",df.stream,envir = globalenv())
  assign(paste0("streamplot_",name),ggplot(df.stream,aes(x=date,y=Qobs))+geom_line(col="black",linetype="dashed")+
   geom_line(col="green",linetype="solid",size=2), envir = globalenv())
  return(list(RMSE=RMSE,NS=NS,KGE=KGE,Qsim=Qsim))
}
### Visualise the results 
streamplot_Try+geom_line(data=streamdf,aes(x=date,y=Qsim),col="red",linetype="dashed")+theme_classic()+
  annotate("text", x=streamdf$date[500], y=50, label= paste("NSE =",round(Testrun.out[[2]],4),"KGE =",round(Testrun.out[[3]],4)),size=5)+
  graph.theme.beta

### GIUH for surface runoff
IUH_comp<-function(gamma1,Ab,dt,delta_T){
  Lag=(gamma1*1.19*(Ab^0.33))/delta_T
  hp=0.8/Lag
  IUHdata=read.table("IUH.txt",header=F)
  t<-IUHdata[,1]*Lag # time lag 
  IUH_0<-IUHdata[,2] * hp #  the dimensionless IUH coordinates
  ti<-seq(0,max(t),dt)
  IUH.GIUH<-interp1(t,IUH_0,ti)  
  return(IUH.GIUH)  # x,y, points to interpolate  
  }
###  Nash IUH for baseflow component 
IUH_NASH<-function(n, gamma1, Ab, dt, delta_T) {
  K=(gamma1*1.19*(Ab^0.33))/delta_T
  time<-seq(0,100,dt) # return a vector 
  IUH.linear<-(time/K)^(n-1) * (exp(-time/K))* (1/(K*factorial(n-1)))
  return(IUH.linear)
}

### Model to account for snow water equivalent
snow_model<-function(precipitation,temperature,temp_min,temp_max,Cm) {
  rainfall<-snowfall<-SWE_snowpack<-SWE_melting<-matrix(nrow=length(precipitation),ncol=1)
  if ( is.nan(precipitation[1,1]) == T | is.nan(temperature[1,1])==T) {snowfall[1,1]=precipitation[1,1]=NaN
  }
  else if (temperature[1,1]<=temp_min) {
    snowfall[1,1]<-precipitation[1,1]
    rainfall[1,1]<-0
    SWE_snowpack[1,1]<-snowfall[1,1]
    SWE_melting[1,1]<-0
  }
  else if (temperature[1,1]>=temp_max)  {
    snowfall[1,1]<-0
    rainfall[1,1]<-precipitation[1,1]
    SWE_snowpack[1,1]<-SWE_melting[1,1]<-0
  }
  else {
    rainfall[1,1] = precipitation[1,1] * ((temperature[1,1]-temp_min)/(temp_max-temp_min));
    snowfall[1,1] = precipitation[1,1] - rainfall[1,1];
    SWE_snowpack[1,1] = snowfall[1,1];
    SWE_melting[1,1] = 0
  }
  # later steps 
  for (i in 2:length(precipitation)) {
  if (is.nan(precipitation[i,1] )== T | is.nan(temperature[i,1])== T) {
  rainfall[i,1] = NaN;
  snowfall[i,1] = NaN; }
  else if (temperature[i,1] <= temp_min) {
  # if the temperature is less than the low threshold, the precipitation is entirely snowfall
  rainfall[i,1] = 0;
  snowfall[i,1] = precipitation[i,1];
  SWE_snowpack[i,1] = SWE_snowpack[i-1,1] + snowfall[i,1];
  SWE_melting[i,1] = 0 }
  else if (temperature[i,1] > temp_max ) {
  # if the temperature is more than the high threshold,the precipitation is entirely rainfall
  rainfall[i,1] = precipitation[i,1];
  snowfall[i,1] = 0;
  SWE_melting[i,1] = Cm * (temperature[i,1] - temp_max);
# h_melting[i,1] = rho_water * SWE_melting[i,1] / rho_snow;
  # Check the snowpack SWE
  if (SWE_snowpack[i-1,1] >= SWE_melting[i,1])
  {SWE_snowpack[i,1] = SWE_snowpack[i-1,1] - SWE_melting[i,1]}
  else {SWE_melting[i,1]= SWE_snowpack[i-1,1]
  SWE_snowpack[i,1] = 0 } }
  
  else { rainfall[i,1] = precipitation[i,1] * ((temperature[i,1]-temp_min)/(temp_max-temp_min));
  snowfall[i,1] = precipitation[i,1] - rainfall[i,1];
  SWE_snowpack[i,1] = SWE_snowpack[i-1,1] + snowfall[i,1];
  SWE_melting[i,1] = 0;
  } } 
  return(list(rainfall,SWE_melting,SWE_snowpack))
}

# KGE for model metric
klingupta<-function(mod, obs) {
  mod[is.na(mod)==T]<-NaN
  flows<-cbind(mod,obs)
  flows<-na.omit(flows)
  sd.mod<-sd(mod,na.rm = T)
  sd.obs<-sd(obs,na.rm = T)
  mean.mod<-mean(mod,na.rm = T)
  mean.sd<-mean(obs,na.rm = T)
  rmat<-cor(data.frame(mod,obs))
  r<-rmat[1,2]
  relvar<-sd.mod/sd.obs
  bias<-mean.mod/mean.sd
  kge<-1-sqrt ( (r-1)**2 + (relvar-1)**2 + (bias-1)**2) 
  return(kge)
}

##### Calibration of model parameters 
convert_adim<-function(X0){
  LOW=c(0.1,100,2,0.1,0.5,0.4,1,0.1/24,5,0.01)
  UP<-c(0.9,3000,10,40,3.5,3.0,15,3,35,65)
  X=(LOW)+(UP-LOW)*X0
  return(X)
} 
# Idea: optimise the weighting parameter x0 between the parameter range
calibOK<-function(x0,data,Ab){
  X<-convert_adim(x0)  
  output<-MISDc.simu(datamat = data,PARS = X,Ab = Ab)
  KGE<-output$KGE
  err<-1-KGE
  return(1-KGE)  # error term to be minimised  
}
### Need to comment out the fixed the parameter value after debugging

cal_MISDc<-function(data,X_init,Ab) {
  NPAR<-10
  Res<-fmincon(x0 = X_init,lb = rep(0,10),ub=rep(1,10),fn = calibOK,data=data,Ab=Ab,maxiter = 300,maxfeval = 1000,tol = 1e-5)
  #Res<-optim(par=as.matrix(X_init),fn=calibOK,method = "L-BFGS-B",lower = rep(0,10),upper=rep(1,10),data=data,Ab=Ab,control = list(maxit=500,trace=1))
  X<-convert_adim(Res$par) # calibrated parameters, transformed to the true parameter range
  outputs<-MISDc.simu(datamat = data,PARS = X,Ab = Ab)  # final run of the model using calibrated parameters 
  return(outputs)
}

### Check the consistency in results: direct calibration & output + input the calibrated values 
cali.modelrun<-cal_MISDc(data = data.mat.new,X_init =X_INIT, Ab = 137)
Testrun.out<-MISDc.simu(data.mat.new,PARS = X_INIT_cali, Ab = 137,name = "Try")
print(paste("NSE =",round(Testrun.out[[2]],4),"KGE =",round(Testrun.out[[3]],4)))
print(paste("NSE =",round(cali.modelrun[[2]],4),"KGE =",round(cali.modelrun[[3]],4)))
### Visualise the results 
streamplot_Try+geom_line(data=streamdf,aes(x=date,y=Qsim),col="red",linetype="dashed")+theme_classic()+
  annotate("text", x=streamdf$date[500], y=50, label= paste("NSE =",round(Testrun.out[[2]],4),"KGE =",round(Testrun.out[[3]],4)),size=5)+
  graph.theme.beta
