library(hydromad)
library(lubridate)
klingupta<-function(mod, obs) {
  flows<-cbind(mod,obs)
  flows<-flows[which(is.na(obs)==F),]
  sd.mod<-sd(mod,na.rm = T)
  sd.obs<-sd(obs,na.rm = T)
  mean.mod<-mean(mod,na.rm = T)
  mean.sd<-mean(obs,na.rm = T)
  rmat<-cor(flows)  # after removing the missing values 
  r<-rmat[1,2]
  relvar<-sd.mod/sd.obs
  bias<-mean.mod/mean.sd
  kge<-1-sqrt ( (r-1)**2 + (relvar-1)**2 + (bias-1)**2) 
  return(kge)    }
hydromad.stats("KGE" = function(Q, X, ...) {
  klingupta(mod = X,obs = Q)
})

### Set up the MISDc options 
hydromad.options(MISDc=list(
  Frac.sz.init=c(0.1,0.9),
  Wmax.rz=c(100,3000),  # total water capacity in root-zone
  m1=c(2,10),      # exponent of drainage of surface layer
  KS1=c(0.1,40),   # hydraulic conductivity of the surface layer
  tau_q=c(0,3),   # recession rate for overland runoff
  tau_s=c(2,100),   # recession rate for baseflow
  Kc=c(0.4,3),     # parameter of potential evapotranspiration
  alpha=c(1,15),   # exponent of runoff
  m2=c(5,35),     # exponent of drainge for root-zone
  KS2=c(0.01,60),   # hydraulic conductivity in root-zone
  Wmax.sz=c(20,300),
  Xs_0=0,Xq_0=0
  ))

### SMA module, but here put the routing module in as well
# Argument "Routeornot" governs whether we want the routing module to execute
### DATA: a zoo object data for the target catchment 


MISDc.sim<-function(DATA,Frac.sz.init,Wmax.rz,m1,KS1,tau_q,tau_s,Kc,alpha,m2,KS2,Wmax.sz,return_state=T,Routeornot=T,Xs_0 = 0, Xq_0 = 0) {

 ### Step 1: check if the input is correct
  stopifnot(c("P", "E") %in% colnames(DATA))
  stopifnot(tau_q>=0)
  stopifnot(tau_s>=0)
  stopifnot(length(which(c(Frac.sz.init,Wmax.rz,m1,KS1,Kc,alpha,m2,KS2,Wmax.sz)<=0))==0)
  xpar=c(Frac.sz.init,Wmax.rz,m1,KS1,Kc,alpha,m2,KS2,Wmax.sz)
  
  ### Step 2: Load Data 
  date=as.Date(index(DATA),"%Y-%m-%d")
  P <- DATA[, "P"]
  E<- DATA[, "E"]
  bad <- is.na(P) | is.na(E)
  P[bad] <- 0
  E[bad] <- 0  # remove all missing data if any present 

  RAIN<-P
  TEMP<-E
  M<-length(date)
  
  #### define the time step: here delta_T = should be 24 hours
  delta_T<-round(mean(diff(date,lag = 1),na.rm=T) *24 *10000)/10000
  delta_T<-as.numeric(delta_T)
  Month<-month(date)
  
  # hydraulic conductivity: mm/h --> mm/delta_T
  KS1<-KS1*(delta_T)
  KS2<-KS2*(delta_T)
  
  # Potential Evapotranspiration parameter over year
  L=rep(0.25,12)
  Ka<-1.26
  PET<-(TEMP>0)*(Kc*(Ka*L[Month]*(0.46*TEMP+8)-2))/(24/delta_T)   # time series of PET
 
   ####  Store the states 
  # WW_sz & WW_rz: relative wetness; BF: baseflow; RO: surface runoff
  # Perc: percolation from surface zone to root-zone (0-1 m)
  # DeepPrec: percolation from root-zone to deeper layers
  # SZ/RZ: absolute wetness of soil (mm)
  # IE: infiltration excess;  SE: saturation excess 
  
  BF<-RO<-WW_sz<-WW_rz<-Perc_all<-DeepPrec_all<-ET_all<-rep(0,M)  
  SZ.all<-RZ.all<-IE_all<-SE_all<-rep(0,M);
  
  ##### ============  SMA module start   ==========
  # Starting conditions of soil wetness index  
  W.sz<-Frac.sz.init*Wmax.sz
  W.rz<-Frac.sz.init*Wmax.rz  
  S<-NaN; IE<-0
  
  ## loop through the time step
  for (t in 1:M) {
    # Direct runoff on first day (A general soil wetness-effective runoff relation)
    IE<-RAIN[t]*((W.sz/Wmax.sz) ^ alpha)      # rainfall is involved here 
    Evap<-PET[t]*(W.sz/Wmax.sz)       
    ET_all[t]<-Evap
    
    if (W.rz<Wmax.rz) {Perc=KS1* ((W.sz/Wmax.sz) ^m1)  }
    else {Perc=0}
    
    ###  Deep percolation
    DeepPrec<-KS2*((W.rz/Wmax.rz) ^m2)
    Perc[Perc>(0.6*Wmax.sz)]=0.6*Wmax.sz     # Additional constraint on the percolation
    Perc_all[t]<-Perc;  
    DeepPrec[DeepPrec>0.6*Wmax.rz] =0.6*Wmax.rz
    DeepPrec_all[t]<-DeepPrec
    
    #### Moisture for both layers (W.sz for surface; W.rz for root-zone)
    W.sz<-max(0,W.sz+RAIN[t]-IE-Perc-Evap)
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
    
    ####  store the results for soil wetness index in time step t 
    WW_sz[t]<-W.sz/Wmax.sz
    WW_rz[t]<-W.rz/Wmax.rz
    
    #### Baseflow 
    BF[t]<-KS2*( ((W.sz+W.rz)/(Wmax.sz+Wmax.rz))^(m2)) # based on drainge component of root-zone
   
   #### Total overland runoff  
    SE<-Excess.sz+Excess.rz
    IE_all[t]<-IE; SE_all[t]<-SE
    RO[t]<-IE+Excess.sz+Excess.rz  # surface runoff + saturation excess
  }
  
  ####### ========== Routing module, by linear transfer function ================
  ##### One store for each flow component, so one parameter for each 
  alpha_s <- exp(-1 / tau_s)
  alpha_q <- exp(-1 / tau_q)
  if (Routeornot) {
  Xq_ts<-stats::filter(RO*(1-alpha_q),alpha_q,method="recursive",init=Xq_0)
  Xs_ts<-stats::filter(BF*(1-alpha_s),alpha_s,method="recursive",init=Xs_0)
  ####    total streamflow
  Qsim<- (Xq_ts+Xs_ts) } else { 
    Qsim<- RO+BF }

  ##########   Store results 
  if  (return_state) {
    if (Routeornot) {
      ## Try to artificially add a U column
     ans<-list(ET=ET_all,Perc=Perc_all,Deep=DeepPrec_all,SZ=SZ.all,RZ=RZ.all,IE=IE_all,SE=SE_all,ET=ET_all,Xq=Xq_ts,Xs=Xs_ts,IE=IE_all,SE=SE_all,
           BF=BF,RO=RO,U=Qsim)
             
       for (i in 1:length(ans) ) { attributes(ans[[i]]) <-attributes(P) }
     ans <- do.call(cbind,ans) 
     return(ans)} else {
      ans.noroute<-list(ET=ET_all,Perc=Perc_all,Deep=DeepPrec_all,WW_sz=WW_sz,WW_rz=WW_rz,SZ=SZ.all,RZ=RZ.all,IE=IE_all,SE=SE_all,
      BF=BF,U=RO,U=Qsim)
      for (i in 1:length(ans.noroute)) { attributes(ans.noroute[[i]]) <-attributes(P) }
     ans <- do.call(cbind,ans.noroute )
     return(ans)
    } }   else{
      U<-Qsim
      attributes(U)<-attributes(P)
      return(U)  }
}

### Read the data for running the model. Data also uploaded in the repo. 

# Try fitByOptim
data.136111.Cali.HMD1<-read.zoo("data.136111.Cali.HMD.csv",header=T,index.column = 1,sep=",")
MISDc.form.136111=hydromad(DATA=data.136111.Cali.HMD,sma="MISDc",warmup = 365,
                           return_state=F,Xs_0=0,Xq_0=0,Routeornot=T)
set.seed(1)
hmd.MISDc.136111.Cali<-fitByOptim(MISDc.form.136111,hmadstat("r.squared"))

######### 10.27 Update: Conduct Sobol' method to assess parameter sensitivity in this model 
library(sensitivity)
library(hydromad)
## Define the KGE star objective function 
hydromad.stats("KGEstar" = function(Q, X,...) {
  KGE<-klingupta(mod = X,obs = Q)
  KGE/(2-KGE)
})

#### Run Sobol' method for global sensitivity analysis
n <- 1000
set.seed(1)
ss.403217.long.KGEstar<- sobol2002(model = evalPars,
                                   ## Draw two random samples of parameters
                                   X1 =parameterSets(getFreeParsRanges(modx),n),
                                   X2 = parameterSets(getFreeParsRanges(modx),n),
                                   ## Number of bootstrap replicates
                                   nboot = 4,
                                   ## Arguments to be passed to evalPars
                                   object=modx,
                                   
                                   objective=~hmadstat("KGEstar")(Q,X)
                                   
)
## Store the results of FSI and TSI 
FSI<-ss.403217.long$V$original
TSI<-ss.403217.long$T$original


