##0.Loading packages and setting directories
require(Rmpfr)
library(RColorBrewer)
library(lattice)
library(grid)
library(gridExtra)
library(ggplot2)
library(plot.matrix)
library(fields)
library(nleqslv)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding")
source("AP11.Plotting_multiple_images.R")

options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=256
ncontour=5
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato"))
palet<-jet.colors(ncol)
pal<-list(palet)

#1.Adaptive dynamics Low toxicity and variations in pathway yield (10:1, 1:1, 1:10)
###1.Influence of degradation rate and fitness contributions on the allocation between subpathways (no influence of transporters)

##Function splitting the pathway with regards to parameters
P1s_func<-function(P1s,Pts,ratio1_2){
  if(missing(P1s)){
    return(Pts*ratio1_2)
  }
  else{
    return(P1s)
  }
}

##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.25##set to a moderately high value
kcat=10^2.25#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value 
PayOff_set<-c(10,1,0.1)
Tox_set=c(10^4)
Path_size_set=c(40)

N_reso=200#Resolution used for strategies (in terms of concentration)
Etot_back<-5.5*10^-3#other protein concentration
E_basal<-3*10^-3#scaling factor

#Adaptive dynamics constants
Phi_eq=10^-4
r_c=1e-5#dm
S_c=4*pi*r_c^2
V_c=4/3*pi*r_c^3#dm^3
beta=10^-3
alpha=10^-3
dl=1e-4#dm
V_env=dl^3

#Parameters
eta_set=10^seq(-4,-1.5,0.5)

#Variables
Etot_var<-10^seq(-6,-4,length=N_reso)

#Outcome
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
Pop_DA<-c()
tab_P2_eq <- c()
tab_S_eq <- c()
tab_N_eq <- list()
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40_noperm<-list()

#Parameter set 1
#PayOff<-PayOff_set[1]
Path_size=Path_size_set[1]
Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)
Tox<-Tox_set[1]

Enz=1
for(po in 1:length(PayOff_set)){
  PayOff<-PayOff_set[po]
  for (p in 1:length(eta_set)){
    eta<-eta_set[p]
    P_p2=P_p2_set[1]
    eta_p=eta+P_p2*S_c/V_c
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
      PO2<-1
    }
    else if (PayOff<1){
      PO1<-1
      PO2<-1/PayOff
    }
    ##Resident strategy and ecological equilibrium
    for (i in 1:N_reso){
      Etot0=10^-4.5##Set concentration approximately to its optimal level
      print(paste(po,p,i))
      Etot1r=Etot_var[i]#Molar
      Etot2r=Etot_var[i]#Molar
      Etot_conc=Etot_var[i]*(Path_size-1)+Etot0
      #Intermediate constants
      Vm0=kcat*Etot0
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kfact=kf*10^(-(Etot_conc+Etot_back)/(E_basal))
      KM=(kr+kcat)/kfact
      N_eq=2
      Phi_fit=0
      P2_out=1
      count=0
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-8) && N_eq>1){
        transport<-function(x){
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-Vm0*x[2]/(x[2]+KM)
          c(t1,t2)
        }
        pathway1<-function(y,Vm1c,Vm2c){
          a_1=eta*(KM+y)
          b_1=Vm2c*(y+KM)+eta*KM*(y+KM)-Vm1c*y
          c_1=-Vm1c*KM*y
          delta_1<-b_1^2-4*a_1*c_1
          P=(-b_1+delta_1^(1/2))/(2*a_1)
          return(P)
        }
        S<-nleqslv(c(1e-6,1e-6),transport,jac=NULL)$x
        Sin<-S[2]
        P1<-c()
        P1[1]<-pathway1(Sin,Vm0,Vm1r)
        for (p1 in 2:(Path1_size)){
          P1[p1]<-pathway1(P1[p1-1],Vm1r,Vm1r)
        }
        P2<-c()
        
        pathway2<-function(z){
          p1 <- Vm1r*P1[Path1_size-1]/(P1[Path1_size-1]+KM)-(Vm2r*z[1]/(z[1]+KM)+eta*z[1]+P_p2*S_c/V_c*(z[1]-z[2]))
          p2 <- N_eq*P_p2*S_c/V_env*(z[1]-z[2])-beta*z[2]
          c(p1,p2)
        }
        P2_sol<-nleqslv(c(1e-6,1e-6),pathway2,jac=NULL)$x
        P2[1]<-P2_sol[1]
        for (p2 in 2:(Path_size-Path1_size)){
          P2[p2]<-pathway1(P2[p2-1],Vm2r,Vm2r)
        }
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm0*Sin/(Sin+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(P1[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2r=Phi2r+(Vm2r*P2[p2]/(P2[p2]+KM))/(Path_size/2)
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc+Etot_back))*Tox/(Tox+sum(P1)+sum(P2))
        ##Algorithm to simulate the population size at demographic equilibrium for the resident strategy
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq/2)
          }
        }
        else if(count<200){
          N_eq=N_eq*Phi_fit/Phi_eq
        }
        else{
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq+1)
          }
          else{
            N_eq=(N_eq-1)
          }
          if(N_eq_p2==N_eq){
            break
          }
        }
        if(count%%2==0){
          N_eq_p2=N_eq
        }
        count=count+1
      }
      Se<-S[1]
      P2_out<-P2_sol[2]
      tab_S_eq[i]<-Se
      tab_P2_eq[i]<-P2_out
      Pop_DA[i]<-N_eq
      ##Recalculating resident fitness (for more accurracy)
      a_res<- (VTm*kfact+kfact*Etot0*kcat*(1+Se/KT))
      b_res<- (VTm*(kr+kcat-Se*kfact)+kfact*kcat*Etot0*(KT+Se))
      c_res<- -VTm*Se*(kr+kcat)
      delta_res<-b_res^2-4*a_res*c_res
      Sin_res=(-b_res+delta_res^(1/2))/(2*a_res)
      a1_res=eta*(KM+Sin_res)
      b1_res=Vm1r*(Sin_res+KM)+eta*KM*(Sin_res+KM)-Vm0*Sin_res
      c1_res=-Vm0*KM*Sin_res
      delta1_res<-b1_res^2-4*a1_res*c1_res
      P1_res<-c()
      P1_res[1]=(-b1_res+delta1_res^(1/2))/(2*a1_res)
      for (p1 in 2:(Path1_size-1)){
        a11_res=eta*(KM+P1_res[p1-1])
        b11_res=Vm1r*(P1_res[p1-1]+KM)+eta*KM*(P1_res[p1-1]+KM)-Vm1r*P1_res[p1-1]
        c11_res=-Vm1r*KM*P1_res[p1-1]
        delta11_res<-b11_res^2-4*a11_res*c11_res
        P1_res[p1]=(-b11_res+delta11_res^(1/2))/(2*a11_res)
      }
      a2_res=eta_p*(KM+P1_res[Path1_size-1])
      b2_res=Vm2r*(P1_res[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_res[Path1_size-1]+KM)-Vm1r*P1_res[Path1_size-1]
      c2_res=-Vm1r*KM*P1_res[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_res[Path1_size-1]+KM)*KM
      delta2_res<-b2_res^2-4*a2_res*c2_res
      P2_res<-c()
      P2_res[1]=(-b2_res+delta2_res^(1/2))/(2*a2_res)
      for (p2 in 2:(Path_size-Path1_size)){
        a22_res=eta*(KM+P2_res[p2-1])
        b22_res=Vm2r*(P2_res[p2-1]+KM)+eta*KM*(P2_res[p2-1]+KM)-Vm2r*P2_res[p2-1]
        c22_res=-Vm2r*KM*P2_res[p2-1]
        delta22_res<-b22_res^2-4*a22_res*c22_res
        P2_res[p2]=(-b22_res+delta22_res^(1/2))/(2*a22_res)
      }
      Phi1r=0
      Phi2r=0
      Phi1r=Phi1r+(Vm0*Sin_res/(Sin_res+KM))/(Path_size/2)
      for (p1 in 1:(Path1_size-1)){
        Phi1r=Phi1r+(Vm1r*P1_res[p1]/(P1_res[p1]+KM))/(Path_size/2)
      }
      for (p2 in 1:(Path_size-Path1_size-1)){
        Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
      }
      Phi_res<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))
      for (j in 1:N_reso){
        E_tot_conc0=10^-4.5
        E_tot_conc1<-Etot_var[j]
        E_tot_conc2<-Etot_var[j]
        Etot_conc<-Etot_var[j]*(Path_size-1)+E_tot_conc0
        kf_act<-kf*10^(-(Etot_conc+(Etot_back))/(E_basal))
        a<- (VTm*kf_act+kf_act*E_tot_conc0*kcat*(1+Se/KT))
        b<- (VTm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc0*(KT+Se))
        c<- -VTm*Se*(kr+kcat)
        delta<-b^2-4*a*c
        #First cellular concentration (substrate)
        if(a==0){
          S_conc_eq=0
        }else{
          S_conc_eq=(-b+delta^(1/2))/(2*a)
        }
        KM=(kcat+kr)/kf_act
        Vm0=kcat*E_tot_conc0
        Vm1=kcat*E_tot_conc1
        Vm2=kcat*E_tot_conc2
        a_1=eta*(KM+S_conc_eq)
        b_1=Vm1*(S_conc_eq+KM)+eta*KM*(S_conc_eq+KM)-Vm0*S_conc_eq
        c_1=-Vm0*KM*S_conc_eq
        delta_1<-b_1^2-4*a_1*c_1
        P1_conc_eq<-c()
        P1_conc_eq[1]=(-b_1+delta_1^(1/2))/(2*a_1)
        
        for (p1 in 2:(Path1_size-1)){
          a_11=eta*(KM+P1_conc_eq[p1-1])
          b_11=Vm1*(P1_conc_eq[p1-1]+KM)+eta*KM*(P1_conc_eq[p1-1]+KM)-Vm1*P1_conc_eq[p1-1]
          c_11=-Vm1*KM*P1_conc_eq[p1-1]
          delta_11<-b_11^2-4*a_11*c_11
          P1_conc_eq[p1]=(-b_11+delta_11^(1/2))/(2*a_11)
        }
        
        a_2=eta_p*(KM+P1_conc_eq[Path1_size-1])
        b_2=Vm2*(P1_conc_eq[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_conc_eq[Path1_size-1]+KM)-Vm1*P1_conc_eq[Path1_size-1]
        c_2=-Vm1*KM*P1_conc_eq[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_conc_eq[Path1_size-1]+KM)*KM
        delta_2<-b_2^2-4*a_2*c_2
        P2_conc_eq<-c()
        P2_conc_eq[1]=(-b_2+delta_2^(1/2))/(2*a_2)
        
        for (p2 in 2:(Path_size-Path1_size)){
          a_22=eta*(KM+P2_conc_eq[p2-1])
          b_22=Vm2*(P2_conc_eq[p2-1]+KM)+eta*KM*(P2_conc_eq[p2-1]+KM)-Vm2*P2_conc_eq[p2-1]
          c_22=-Vm2*KM*P2_conc_eq[p2-1]
          delta_22<-b_22^2-4*a_22*c_22
          P2_conc_eq[p2]=(-b_22+delta_22^(1/2))/(2*a_22)
        }
        Phi1=0
        Phi2=0
        Phi1=Phi1+(Vm0*S_conc_eq/(S_conc_eq+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1=Phi1+(Vm1*P1_conc_eq[p1]/(P1_conc_eq[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
        }
        if(as.numeric((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc+(Etot_back)))>0){
          P_Phi_eq=as.numeric(((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc+(Etot_back)))*Tox/(Tox+sum(P1_conc_eq)+sum(P2_conc_eq)))
        }
        else{
          P_Phi_eq=0
        }
        if(N_eq>1){
          fit=P_Phi_eq/Phi_res-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[i,j]=P_Phi_eq
        fit_inv[i,j]=fit
      }
    }
    list_fit_inv_S40_noperm[[p+length(eta_set)*(po-1)]]<-fit_inv
    tab_N_eq[[p+length(eta_set)*(po-1)]]<-max(Pop_DA)
  }
}
addtxt<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
for (i in 1:N_reso){
  print(fit_inv[i,i])
}
ncol=800
jet.colors <- colorRampPalette(c(rep("red",400),rep("green",400)))
palet<-jet.colors(ncol)
pal<-list(palet)
f_var_print<-c(0,1,0.2)
log10Etot1_var<-c(-5,-4,0.5)
sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
             expression(paste(eta,"=",10^-2,"/s")))
multiplePlot("","","","",c("","",""),c("","","","","","")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm,
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

Conc1_opt_DA_E0E1comb<-list()
Pox<-c("10","1","0.1")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for(po in 1:length(PayOff_set)){
  Pox_i<-as.character(PayOff_set[po])
  Conc1_opt_DA_E0E1comb[[Pox_i]]<-c()
  for (e in 1:length(eta_set)){
      i=1
      while(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]>0){
            Conc1_opt_DA_E0E1comb[[Tox_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-0
            break
          }
        }
        else{
          i=i+1
        }
      }
    }
  }

for(po in 1:length(PayOff_set)){
  Pox_i<-as.character(PayOff_set[po])
  Conc1_opt_DA_E0E1comb[[Pox_i]]<-c()
  for (e in 1:length(eta_set)){
    i=1
    while(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    while (i<N_reso){
      if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-1]<0){
        print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
        #print(Etot_var[S_opt])
        break
      }
      else if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+2]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-2]<0 && i>2 && i<N_reso-2){
        print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
        print(Etot_var[S_opt])
        break
      }
      else if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+3]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-3]<0 && i>3 && i<N_reso-3){
        print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
        print(Etot_var[S_opt])
        break
      }
      else if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+5]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-5]<0 && i>5 && i<N_reso-5){
        print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
        print(Etot_var[S_opt])
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]>0){
          Conc1_opt_DA_E0E1comb[[Pox_i]][e]=10^-4
          i=i+1
          break
        }
        else{
          Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-0
          break
        }
      }
      else{
        i=i+1
      }
    }
  }
}


#2.Adaptive dynamics with 2 sub-pathways

##Constant
ratio_set<-c(0.5) ##Size of each subpathway: 0.5 yields same size for both of them

##Redefining space of adaptive parameters to coincide with an allocation
N_reso=200
f_var<-seq(0,1,length=N_reso)

##New outcome list with new resolution
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40ratvar_Y1_eta1<-list()

for (po in 1:length(PayOff_set)){
  PayOff<-PayOff_set[po]
  ratio_12<-ratio_set[1]
  Path1_size=P1s_func(Pts=Path_size,ratio1_2=ratio_12)
  for (p in 1:length(eta_set)){
    Etot_conc=Conc1_opt_DA_E0E1comb[[po]][p]*2#Double so that investing half in each sub-pathway coincides with the optimum when considered as a single integrated pathway
    eta=eta_set[p]
    eta_p=eta+P_p2*S_c/V_c
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
      PO2<-1
    }
    else if (PayOff<1){
      PO1<-1
      PO2<-1/PayOff
    }
    ##Resident strategy and ecological equilibrium
    for (i in 1:N_reso){
      Etot0=10^-4.5 #First concentration set to avoid excess of upstreamn investment due to transport
      print(paste(po,p,i))
      Etot1r=Etot_conc*f_var[i]#Molar
      Etot2r=Etot_conc*(1-f_var[i])#Molar
      Etot_conc_r=Etot1r*(Path1_size-1)+Etot2r*(Path_size-Path1_size)+Etot0
      #Intermediate constants
      Vm0=kcat*Etot0
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kfact=kf*10^(-(Etot_conc_r+Etot_back)/(E_basal))
      KM=(kr+kcat)/kfact
      N_eq=10
      Phi_fit=0
      P2_out=1
      count=0
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-8) && N_eq>1){
        transport<-function(x){
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-Vm0*x[2]/(x[2]+KM)
          c(t1,t2)
        }
        pathway1<-function(y,Vm1c,Vm2c){
          a_1=eta*(KM+y)
          b_1=Vm2c*(y+KM)+eta*KM*(y+KM)-Vm1c*y
          c_1=-Vm1c*KM*y
          delta_1<-b_1^2-4*a_1*c_1
          P=(-b_1+delta_1^(1/2))/(2*a_1)
          return(P)
        }
        S<-nleqslv(c(1e-6,1e-6),transport,jac=NULL)$x
        Sin<-S[2]
        P1<-c()
        P1[1]<-pathway1(Sin,Vm0,Vm1r)
        for (p1 in 2:(Path1_size-1)){
          P1[p1]<-pathway1(P1[p1-1],Vm1r,Vm1r)
        }
        P2<-c()
        pathway2<-function(z){
          p1 <- Vm1r*P1[Path1_size-1]/(P1[Path1_size-1]+KM)-(Vm2r*z[1]/(z[1]+KM)+eta*z[1]+P_p2*S_c/V_c*(z[1]-z[2]))
          p2 <- N_eq*P_p2*S_c/V_env*(z[1]-z[2])-beta*z[2]
          c(p1,p2)
        }
        P2_sol<-nleqslv(c(1e-6,1e-6),pathway2,jac=NULL)$x
        P2[1]<-P2_sol[1]
        for (p2 in 2:(Path_size-Path1_size)){
          P2[p2]<-pathway1(P2[p2-1],Vm2r,Vm2r)
        }
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm0*Sin/(Sin+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(P1[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2r=Phi2r+(Vm2r*P2[p2]/(P2[p2]+KM))/(Path_size/2)
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1)+sum(P2))
        ##Algorithm to simulate the population size at demographic equilibrium for the resident strategy
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq/2)
          }
        }
        else if(count<200){
          N_eq=N_eq*Phi_fit/Phi_eq
        }
        else{
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq+1)
          }
          else{
            N_eq=(N_eq-1)
          }
          if(N_eq_p2==N_eq){
            break
          }
        }
        if(count%%2==0){
          N_eq_p2=N_eq
        }
        count=count+1
      }
      Se<-S[1]
      P2_out<-P2_sol[2]
      tab_S_eq[i]<-Se
      tab_P2_eq[i]<-P2_out
      Pop_DA[i]<-N_eq
      ##Recalculating resident fitness (for more accurracy)
      a_res<- (VTm*kfact+kfact*Etot0*kcat*(1+Se/KT))
      b_res<- (VTm*(kr+kcat-Se*kfact)+kfact*kcat*Etot0*(KT+Se))
      c_res<- -VTm*Se*(kr+kcat)
      delta_res<-b_res^2-4*a_res*c_res
      Sin_res=(-b_res+delta_res^(1/2))/(2*a_res)
      a1_res=eta*(KM+Sin_res)
      b1_res=Vm1r*(Sin_res+KM)+eta*KM*(Sin_res+KM)-Vm0*Sin_res
      c1_res=-Vm0*KM*Sin_res
      delta1_res<-b1_res^2-4*a1_res*c1_res
      P1_res<-c()
      P1_res[1]=(-b1_res+delta1_res^(1/2))/(2*a1_res)
      for (p1 in 2:(Path1_size-1)){
        a11_res=eta*(KM+P1_res[p1-1])
        b11_res=Vm1r*(P1_res[p1-1]+KM)+eta*KM*(P1_res[p1-1]+KM)-Vm1r*P1_res[p1-1]
        c11_res=-Vm1r*KM*P1_res[p1-1]
        delta11_res<-b11_res^2-4*a11_res*c11_res
        P1_res[p1]=(-b11_res+delta11_res^(1/2))/(2*a11_res)
      }
      a2_res=eta_p*(KM+P1_res[Path1_size-1])
      b2_res=Vm2r*(P1_res[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_res[Path1_size-1]+KM)-Vm1r*P1_res[Path1_size-1]
      c2_res=-Vm1r*KM*P1_res[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_res[Path1_size-1]+KM)*KM
      delta2_res<-b2_res^2-4*a2_res*c2_res
      P2_res<-c()
      P2_res[1]=(-b2_res+delta2_res^(1/2))/(2*a2_res)
      for (p2 in 2:(Path_size-Path1_size)){
        a22_res=eta*(KM+P2_res[p2-1])
        b22_res=Vm2r*(P2_res[p2-1]+KM)+eta*KM*(P2_res[p2-1]+KM)-Vm2r*P2_res[p2-1]
        c22_res=-Vm2r*KM*P2_res[p2-1]
        delta22_res<-b22_res^2-4*a22_res*c22_res
        P2_res[p2]=(-b22_res+delta22_res^(1/2))/(2*a22_res)
      }
      Phi1r=0
      Phi2r=0
      Phi1r=Phi1r+(Vm0*Sin_res/(Sin_res+KM))/(Path_size/2)
      for (p1 in 1:(Path1_size-1)){
        Phi1r=Phi1r+(Vm1r*P1_res[p1]/(P1_res[p1]+KM))/(Path_size/2)
      }
      for (p2 in 1:(Path_size-Path1_size-1)){
        Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
      }
      Phi_res<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))
      for (j in 1:N_reso){
        E_tot_conc0=10^-4.5
        E_tot_conc1<-Etot_conc*f_var[j]
        E_tot_conc2<-Etot_conc*(1-f_var[j])
        Etot_conc_m=E_tot_conc1*(Path1_size-1)+E_tot_conc2*(Path_size-Path1_size)+E_tot_conc0
        kf_act<-kf*10^(-(Etot_conc_m+Etot_back)/(E_basal))
        a<- (VTm*kf_act+kf_act*E_tot_conc0*kcat*(1+Se/KT))
        b<- (VTm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc0*(KT+Se))
        c<- -VTm*Se*(kr+kcat)
        delta<-b^2-4*a*c
        #First cellular concentration (substrate)
        if(a==0){
          S_conc_eq=0
        }else{
          S_conc_eq=(-b+delta^(1/2))/(2*a)
        }
        KM=(kcat+kr)/kf_act
        Vm0=kcat*E_tot_conc0
        Vm1=kcat*E_tot_conc1
        Vm2=kcat*E_tot_conc2
        a_1=eta*(KM+S_conc_eq)
        b_1=Vm1*(S_conc_eq+KM)+eta*KM*(S_conc_eq+KM)-Vm0*S_conc_eq
        c_1=-Vm0*KM*S_conc_eq
        delta_1<-b_1^2-4*a_1*c_1
        P1_conc_eq<-c()
        P1_conc_eq[1]=(-b_1+delta_1^(1/2))/(2*a_1)
        for (p1 in 2:(Path1_size-1)){
          a_11=eta*(KM+P1_conc_eq[p1-1])
          b_11=Vm1*(P1_conc_eq[p1-1]+KM)+eta*KM*(P1_conc_eq[p1-1]+KM)-Vm1*P1_conc_eq[p1-1]
          c_11=-Vm1*KM*P1_conc_eq[p1-1]
          delta_11<-b_11^2-4*a_11*c_11
          P1_conc_eq[p1]=(-b_11+delta_11^(1/2))/(2*a_11)
        }
        
        a_2=eta_p*(KM+P1_conc_eq[Path1_size-1])
        b_2=Vm2*(P1_conc_eq[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_conc_eq[Path1_size-1]+KM)-Vm1*P1_conc_eq[Path1_size-1]
        c_2=-Vm1*KM*P1_conc_eq[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_conc_eq[Path1_size-1]+KM)*KM
        delta_2<-b_2^2-4*a_2*c_2
        P2_conc_eq<-c()
        P2_conc_eq[1]=(-b_2+delta_2^(1/2))/(2*a_2)
        
        for (p2 in 2:(Path_size-Path1_size)){
          a_22=eta*(KM+P2_conc_eq[p2-1])
          b_22=Vm2*(P2_conc_eq[p2-1]+KM)+eta*KM*(P2_conc_eq[p2-1]+KM)-Vm2*P2_conc_eq[p2-1]
          c_22=-Vm2*KM*P2_conc_eq[p2-1]
          delta_22<-b_22^2-4*a_22*c_22
          P2_conc_eq[p2]=(-b_22+delta_22^(1/2))/(2*a_22)
        }
        Phi1=0
        Phi2=0
        Phi1=Phi1+(Vm0*S_conc_eq/(S_conc_eq+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1=Phi1+(Vm1*P1_conc_eq[p1]/(P1_conc_eq[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
        }
        
        if(as.numeric((Phi1*PO1+Phi2*PO2)-prot_cost*((E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size))+(Etot_back)))>0){
          P_Phi_eq=as.numeric(((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc_m+Etot_back))*Tox/(Tox+sum(P1_conc_eq)+sum(P2_conc_eq)))
        }
        else{
          P_Phi_eq=0
        }
        if(N_eq>1){
          fit=P_Phi_eq/Phi_res-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[i,j]=P_Phi_eq
        fit_inv[i,j]=fit
      }
    }
    list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]]<-fit_inv
  }
}

ncol=800
jet.colors <- colorRampPalette(c(rep("white",400),rep("grey30",400)))
palet<-jet.colors(ncol)
pal<-list(palet)
Pop_DA1<-Pop_DA
multiplePlot("Main=","","Deg. rate=","",c("1st","Both","2nd"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5")
             ,ncol=800,f_var_print,f_var_print,list_fit_inv_S40ratvar_Y1_eta1,
             abs="res. (% 1st pathway)",ord="mut. (% 1st pathway)",scale=c(-80,80),lev=c(0),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


dev.print(device = jpeg, file = "res-ToxDeg-2P-PO-var-enz25.jpeg", , width = 2100*3,height=800*3,res=300,type="cairo")
dev.off()
multiplePlot("Main=","","Deg. rate=","",c(""),c("")
             ,ncol=800,f_var_print,f_var_print,list(list_fit_inv_S40ratvar_Y1_eta1[[3]][152:159,152:159]),
             abs="res. (% 1st pathway)",ord="mut. (% 1st pathway)",scale=c(-80,80),lev=c(0),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


Rep_opt_DA_Deg_S40_eta_sens<-list()
for (po in 1:length(Pox)){
  Rep_opt_DA_Deg_S40_eta_sens[[po]]<-c(0)
  for (p in 1:length(eta_set)){
    #print(p)
    i=1
    while(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      #print(i)
      if(i>1 && i<(N_reso-1)){
        if(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i+1,i]>0 && list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i-1,i]>0){
          print(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i])
          S_opt=i
          print(S_opt)
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]<-f_var[S_opt]
          break
        }
      }
      else if (i==(N_reso-1)){
        if(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i+1,i]<0 && list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i+1]>0){
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]=1
          i=i+1
          break
        }
        else{
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]<-0
          break
        }
      }
      i=i+1
    }
  }
}

col_set<-c("darkgrey","pink","palegreen")
par(mfrow=c(1,1))
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875),xpd=FALSE)
for(po in 1:length(Pox)){
  Pox_i<-as.character(Pox[po])
  if (po==1){
    plot(Rep_opt_DA_Deg_S40_eta_sens[[1]]~log10(eta_set),col=col_set[po],pch=16,ylim=c(0.25,1),xlab="log10(Degradation rate)",ylab="(% invetsment in 1st pathway)",cex=2.25)
  }
  else{
    points(Rep_opt_DA_Deg_S40_eta_sens[[po]]~log10(eta_set),col=col_set[po],pch=16,cex=2.25)
  }
}
for(po in 1:length(Pox)){
  Pox_i<-as.character(Pox[po])
  points(Rep_opt_DA_Deg_S40_eta_sens[[po]]~log10(eta_set),col=po,pch=16,cex=1)
}
text(-4,1,"A",srt=addtxt$srt,font=2,col="black",cex=1)
abline(h=0.5,col="gray 10",lty=4)
text(-1.6,0.54,"50%",srt=addtxt$srt,font=1,col="black",cex=1)
legend("bottomright",title="Pay-off",legend=c("10:1","1:1","1:10"),bty="y",pch=c(16),col=c(1:3),ncol=1,cex=0.85)
legend("bottomleft",title="Singular strat.",legend=c("Loc.stable","Glob.stab."),bty="y",pch=c(16,16),col="gray20",ncol=1,pt.cex=c(0.85,1.5))

setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-ToxDeg-2P-PO-gen25.jpeg", width = 560*3,height=475*3,res=300,type="cairo")
dev.off()

###2.Influence of transporters on the allocation between subpathways

##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.25##set to a moderately high value
kcat=10^2.25#/s##set to a moderately high value
kr=kcat#
prot_cost_set=10^-2.5#Mmet/Menz##set to an average value where 
PayOff_set<-c(10,1,0.1)
Tox_set=c(10^4)
Path_size_set=c(40)

N_reso=200#Resolution used for strategies (in terms of concentration)
Etot_back<-5.5*10^-3#other protein concentration
E_basal<-3*10^-3#scaling factor

#Adaptive dynamics constants
Phi_eq=10^-4
r_c=1e-5#dm
S_c=4*pi*r_c^2
V_c=4/3*pi*r_c^3#dm^3
beta=10^-3
alpha=10^-3
dl=1e-4#dm
V_env=dl^3

#Parameters
eta_set=10^seq(-4,-1.5,0.5)

#Variables
Etot_var<-10^seq(-6,-4,length=N_reso)

#Outcome
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
Pop_DA<-c()
tab_P2_eq <- c()
tab_S_eq <- c()
tab_N_eq <- list()
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40_noperm<-list()

#Parameter set 1
#PayOff<-PayOff_set[1]
Path_size=Path_size_set[1]
Path1_size=Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)
Tox<-Tox_set[1]

Enz=1
for(po in 1:length(PayOff_set)){
  PayOff<-PayOff_set[po]
  for (p in 1:length(eta_set)){
    eta<-eta_set[p]
    P_p2=P_p2_set[1]
    eta_p=eta+P_p2*S_c/V_c
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
      PO2<-1
    }
    else if (PayOff<1){
      PO1<-1
      PO2<-1/PayOff
    }
    ##Resident strategy and ecological equilibrium
    for (i in 1:N_reso){
      Etot0=Etot_var[i]##Set concentration approximately to its optimal level
      print(paste(po,p,i))
      Etot1r=Etot_var[i]#Molar
      Etot2r=Etot_var[i]#Molar
      Etot_conc=Etot_var[i]*(Path_size-1)+Etot0
      #Intermediate constants
      Vm0=kcat*Etot0
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kfact=kf*10^(-(Etot_conc+Etot_back)/(E_basal))
      KM=(kr+kcat)/kfact
      N_eq=2
      Phi_fit=0
      P2_out=1
      count=0
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-8) && N_eq>1){
        transport<-function(x){
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-Vm0*x[2]/(x[2]+KM)
          c(t1,t2)
        }
        pathway1<-function(y,Vm1c,Vm2c){
          a_1=eta*(KM+y)
          b_1=Vm2c*(y+KM)+eta*KM*(y+KM)-Vm1c*y
          c_1=-Vm1c*KM*y
          delta_1<-b_1^2-4*a_1*c_1
          P=(-b_1+delta_1^(1/2))/(2*a_1)
          return(P)
        }
        S<-nleqslv(c(1e-6,1e-6),transport,jac=NULL)$x
        Sin<-S[2]
        P1<-c()
        P1[1]<-pathway1(Sin,Vm0,Vm1r)
        for (p1 in 2:(Path1_size)){
          P1[p1]<-pathway1(P1[p1-1],Vm1r,Vm1r)
        }
        P2<-c()
        
        pathway2<-function(z){
          p1 <- Vm1r*P1[Path1_size-1]/(P1[Path1_size-1]+KM)-(Vm2r*z[1]/(z[1]+KM)+eta*z[1]+P_p2*S_c/V_c*(z[1]-z[2]))
          p2 <- N_eq*P_p2*S_c/V_env*(z[1]-z[2])-beta*z[2]
          c(p1,p2)
        }
        P2_sol<-nleqslv(c(1e-6,1e-6),pathway2,jac=NULL)$x
        P2[1]<-P2_sol[1]
        for (p2 in 2:(Path_size-Path1_size)){
          P2[p2]<-pathway1(P2[p2-1],Vm2r,Vm2r)
        }
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm0*Sin/(Sin+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(P1[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2r=Phi2r+(Vm2r*P2[p2]/(P2[p2]+KM))/(Path_size/2)
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc+Etot_back))*Tox/(Tox+sum(P1)+sum(P2))
        ##Algorithm to simulate the population size at demographic equilibrium for the resident strategy
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq/2)
          }
        }
        else if(count<200){
          N_eq=N_eq*Phi_fit/Phi_eq
        }
        else{
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq+1)
          }
          else{
            N_eq=(N_eq-1)
          }
          if(N_eq_p2==N_eq){
            break
          }
        }
        if(count%%2==0){
          N_eq_p2=N_eq
        }
        count=count+1
      }
      Se<-S[1]
      P2_out<-P2_sol[2]
      tab_S_eq[i]<-Se
      tab_P2_eq[i]<-P2_out
      Pop_DA[i]<-N_eq
      ##Recalculating resident fitness (for more accurracy)
      a_res<- (VTm*kfact+kfact*Etot0*kcat*(1+Se/KT))
      b_res<- (VTm*(kr+kcat-Se*kfact)+kfact*kcat*Etot0*(KT+Se))
      c_res<- -VTm*Se*(kr+kcat)
      delta_res<-b_res^2-4*a_res*c_res
      Sin_res=(-b_res+delta_res^(1/2))/(2*a_res)
      a1_res=eta*(KM+Sin_res)
      b1_res=Vm1r*(Sin_res+KM)+eta*KM*(Sin_res+KM)-Vm0*Sin_res
      c1_res=-Vm0*KM*Sin_res
      delta1_res<-b1_res^2-4*a1_res*c1_res
      P1_res<-c()
      P1_res[1]=(-b1_res+delta1_res^(1/2))/(2*a1_res)
      for (p1 in 2:(Path1_size-1)){
        a11_res=eta*(KM+P1_res[p1-1])
        b11_res=Vm1r*(P1_res[p1-1]+KM)+eta*KM*(P1_res[p1-1]+KM)-Vm1r*P1_res[p1-1]
        c11_res=-Vm1r*KM*P1_res[p1-1]
        delta11_res<-b11_res^2-4*a11_res*c11_res
        P1_res[p1]=(-b11_res+delta11_res^(1/2))/(2*a11_res)
      }
      a2_res=eta_p*(KM+P1_res[Path1_size-1])
      b2_res=Vm2r*(P1_res[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_res[Path1_size-1]+KM)-Vm1r*P1_res[Path1_size-1]
      c2_res=-Vm1r*KM*P1_res[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_res[Path1_size-1]+KM)*KM
      delta2_res<-b2_res^2-4*a2_res*c2_res
      P2_res<-c()
      P2_res[1]=(-b2_res+delta2_res^(1/2))/(2*a2_res)
      for (p2 in 2:(Path_size-Path1_size)){
        a22_res=eta*(KM+P2_res[p2-1])
        b22_res=Vm2r*(P2_res[p2-1]+KM)+eta*KM*(P2_res[p2-1]+KM)-Vm2r*P2_res[p2-1]
        c22_res=-Vm2r*KM*P2_res[p2-1]
        delta22_res<-b22_res^2-4*a22_res*c22_res
        P2_res[p2]=(-b22_res+delta22_res^(1/2))/(2*a22_res)
      }
      Phi1r=0
      Phi2r=0
      Phi1r=Phi1r+(Vm0*Sin_res/(Sin_res+KM))/(Path_size/2)
      for (p1 in 1:(Path1_size-1)){
        Phi1r=Phi1r+(Vm1r*P1_res[p1]/(P1_res[p1]+KM))/(Path_size/2)
      }
      for (p2 in 1:(Path_size-Path1_size-1)){
        Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
      }
      Phi_res<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))
      for (j in 1:N_reso){
        E_tot_conc0=Etot_var[j]
        E_tot_conc1<-Etot_var[j]
        E_tot_conc2<-Etot_var[j]
        Etot_conc<-Etot_var[j]*(Path_size-1)+E_tot_conc0
        kf_act<-kf*10^(-(Etot_conc+(Etot_back))/(E_basal))
        a<- (VTm*kf_act+kf_act*E_tot_conc0*kcat*(1+Se/KT))
        b<- (VTm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc0*(KT+Se))
        c<- -VTm*Se*(kr+kcat)
        delta<-b^2-4*a*c
        #First cellular concentration (substrate)
        if(a==0){
          S_conc_eq=0
        }else{
          S_conc_eq=(-b+delta^(1/2))/(2*a)
        }
        KM=(kcat+kr)/kf_act
        Vm0=kcat*E_tot_conc0
        Vm1=kcat*E_tot_conc1
        Vm2=kcat*E_tot_conc2
        a_1=eta*(KM+S_conc_eq)
        b_1=Vm1*(S_conc_eq+KM)+eta*KM*(S_conc_eq+KM)-Vm0*S_conc_eq
        c_1=-Vm0*KM*S_conc_eq
        delta_1<-b_1^2-4*a_1*c_1
        P1_conc_eq<-c()
        P1_conc_eq[1]=(-b_1+delta_1^(1/2))/(2*a_1)
        
        for (p1 in 2:(Path1_size-1)){
          a_11=eta*(KM+P1_conc_eq[p1-1])
          b_11=Vm1*(P1_conc_eq[p1-1]+KM)+eta*KM*(P1_conc_eq[p1-1]+KM)-Vm1*P1_conc_eq[p1-1]
          c_11=-Vm1*KM*P1_conc_eq[p1-1]
          delta_11<-b_11^2-4*a_11*c_11
          P1_conc_eq[p1]=(-b_11+delta_11^(1/2))/(2*a_11)
        }
        
        a_2=eta_p*(KM+P1_conc_eq[Path1_size-1])
        b_2=Vm2*(P1_conc_eq[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_conc_eq[Path1_size-1]+KM)-Vm1*P1_conc_eq[Path1_size-1]
        c_2=-Vm1*KM*P1_conc_eq[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_conc_eq[Path1_size-1]+KM)*KM
        delta_2<-b_2^2-4*a_2*c_2
        P2_conc_eq<-c()
        P2_conc_eq[1]=(-b_2+delta_2^(1/2))/(2*a_2)
        
        for (p2 in 2:(Path_size-Path1_size)){
          a_22=eta*(KM+P2_conc_eq[p2-1])
          b_22=Vm2*(P2_conc_eq[p2-1]+KM)+eta*KM*(P2_conc_eq[p2-1]+KM)-Vm2*P2_conc_eq[p2-1]
          c_22=-Vm2*KM*P2_conc_eq[p2-1]
          delta_22<-b_22^2-4*a_22*c_22
          P2_conc_eq[p2]=(-b_22+delta_22^(1/2))/(2*a_22)
        }
        Phi1=0
        Phi2=0
        Phi1=Phi1+(Vm0*S_conc_eq/(S_conc_eq+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1=Phi1+(Vm1*P1_conc_eq[p1]/(P1_conc_eq[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
        }
        if(as.numeric((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc+(Etot_back)))>0){
          P_Phi_eq=as.numeric(((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc+(Etot_back)))*Tox/(Tox+sum(P1_conc_eq)+sum(P2_conc_eq)))
        }
        else{
          P_Phi_eq=0
        }
        if(N_eq>1){
          fit=P_Phi_eq/Phi_res-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[i,j]=P_Phi_eq
        fit_inv[i,j]=fit
      }
    }
    list_fit_inv_S40_noperm[[p+length(eta_set)*(po-1)]]<-fit_inv
    tab_N_eq[[p+length(eta_set)*(po-1)]]<-max(Pop_DA)
  }
}
addtxt<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
for (i in 1:N_reso){
  print(fit_inv[i,i])
}
ncol=800
jet.colors <- colorRampPalette(c(rep("red",400),rep("green",400)))
palet<-jet.colors(ncol)
pal<-list(palet)
f_var_print<-c(0,1,0.2)
log10Etot1_var<-c(-5,-4,0.5)
sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),
             expression(paste(eta,"=",10^-2,"/s")))
multiplePlot("","","","",c("","",""),c("","","","","","")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm,
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

Conc1_opt_DA_E0E1comb<-list()
Pox<-c("10","1","0.1")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for(po in 1:length(PayOff_set)){
  Pox_i<-as.character(PayOff_set[po])
  Conc1_opt_DA_E0E1comb[[Pox_i]]<-c()
  for (e in 1:length(eta_set)){
    i=1
    while(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i]==-100 && i<N_reso){
      i=i+1
      if(i==N_reso){
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-10^-10
      }
    }
    S_opt=0
    while (i<N_reso-1){
      if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i-1]<0){
        print(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-Etot_var[S_opt]
        print(Etot_var[S_opt])
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[e+length(eta_set)*(po-1)]][i,i+1]>0){
          Conc1_opt_DA_E0E1comb[[Tox_i]][p]=10^-4
          i=i+1
          break
        }
        else{
          Conc1_opt_DA_E0E1comb[[Pox_i]][e]<-0
          break
        }
      }
      else{
        i=i+1
      }
    }
  }
}


#2.Adaptive dynamics with 2 sub-pathways

##Constant
ratio_set<-c(0.5) ##Size of each subpathway: 0.5 yields same size for both of them

##Redefining space of adaptive parameters to coincide with an allocation
N_reso=200
f_var<-seq(0,1,length=N_reso)

##New outcome list with new resolution
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40ratvar_Y1_eta1<-list()

for (po in 1:length(PayOff_set)){
  PayOff<-PayOff_set[po]
  ratio_12<-ratio_set[1]
  Path1_size=P1s_func(Pts=Path_size,ratio1_2=ratio_12)
  for (p in 1:length(eta_set)){
    Etot_conc=Conc1_opt_DA_E0E1comb[[po]][p]*2#Double so that investing half in each sub-pathway coincides with the optimum when considered as a single integrated pathway
    eta=eta_set[p]
    eta_p=eta+P_p2*S_c/V_c
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
      PO2<-1
    }
    else if (PayOff<1){
      PO1<-1
      PO2<-1/PayOff
    }
    ##Resident strategy and ecological equilibrium
    for (i in 1:N_reso){
      Etot0=Etot_conc*f_var[i] #First concentration set to avoid excess of upstreamn investment due to transport
      print(paste(po,p,i))
      Etot1r=Etot_conc*f_var[i]#Molar
      Etot2r=Etot_conc*(1-f_var[i])#Molar
      Etot_conc_r=Etot1r*(Path1_size-1)+Etot2r*(Path_size-Path1_size)+Etot0
      #Intermediate constants
      Vm0=kcat*Etot0
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kfact=kf*10^(-(Etot_conc_r+Etot_back)/(E_basal))
      KM=(kr+kcat)/kfact
      N_eq=10
      Phi_fit=0
      P2_out=1
      count=0
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-8) && N_eq>1){
        transport<-function(x){
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-Vm0*x[2]/(x[2]+KM)
          c(t1,t2)
        }
        pathway1<-function(y,Vm1c,Vm2c){
          a_1=eta*(KM+y)
          b_1=Vm2c*(y+KM)+eta*KM*(y+KM)-Vm1c*y
          c_1=-Vm1c*KM*y
          delta_1<-b_1^2-4*a_1*c_1
          P=(-b_1+delta_1^(1/2))/(2*a_1)
          return(P)
        }
        S<-nleqslv(c(1e-6,1e-6),transport,jac=NULL)$x
        Sin<-S[2]
        P1<-c()
        P1[1]<-pathway1(Sin,Vm0,Vm1r)
        for (p1 in 2:(Path1_size-1)){
          P1[p1]<-pathway1(P1[p1-1],Vm1r,Vm1r)
        }
        P2<-c()
        pathway2<-function(z){
          p1 <- Vm1r*P1[Path1_size-1]/(P1[Path1_size-1]+KM)-(Vm2r*z[1]/(z[1]+KM)+eta*z[1]+P_p2*S_c/V_c*(z[1]-z[2]))
          p2 <- N_eq*P_p2*S_c/V_env*(z[1]-z[2])-beta*z[2]
          c(p1,p2)
        }
        P2_sol<-nleqslv(c(1e-6,1e-6),pathway2,jac=NULL)$x
        P2[1]<-P2_sol[1]
        for (p2 in 2:(Path_size-Path1_size)){
          P2[p2]<-pathway1(P2[p2-1],Vm2r,Vm2r)
        }
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm0*Sin/(Sin+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(P1[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2r=Phi2r+(Vm2r*P2[p2]/(P2[p2]+KM))/(Path_size/2)
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1)+sum(P2))
        ##Algorithm to simulate the population size at demographic equilibrium for the resident strategy
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq/2)
          }
        }
        else if(count<200){
          N_eq=N_eq*Phi_fit/Phi_eq
        }
        else{
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq+1)
          }
          else{
            N_eq=(N_eq-1)
          }
          if(N_eq_p2==N_eq){
            break
          }
        }
        if(count%%2==0){
          N_eq_p2=N_eq
        }
        count=count+1
      }
      Se<-S[1]
      P2_out<-P2_sol[2]
      tab_S_eq[i]<-Se
      tab_P2_eq[i]<-P2_out
      Pop_DA[i]<-N_eq
      ##Recalculating resident fitness (for more accurracy)
      a_res<- (VTm*kfact+kfact*Etot0*kcat*(1+Se/KT))
      b_res<- (VTm*(kr+kcat-Se*kfact)+kfact*kcat*Etot0*(KT+Se))
      c_res<- -VTm*Se*(kr+kcat)
      delta_res<-b_res^2-4*a_res*c_res
      Sin_res=(-b_res+delta_res^(1/2))/(2*a_res)
      a1_res=eta*(KM+Sin_res)
      b1_res=Vm1r*(Sin_res+KM)+eta*KM*(Sin_res+KM)-Vm0*Sin_res
      c1_res=-Vm0*KM*Sin_res
      delta1_res<-b1_res^2-4*a1_res*c1_res
      P1_res<-c()
      P1_res[1]=(-b1_res+delta1_res^(1/2))/(2*a1_res)
      for (p1 in 2:(Path1_size-1)){
        a11_res=eta*(KM+P1_res[p1-1])
        b11_res=Vm1r*(P1_res[p1-1]+KM)+eta*KM*(P1_res[p1-1]+KM)-Vm1r*P1_res[p1-1]
        c11_res=-Vm1r*KM*P1_res[p1-1]
        delta11_res<-b11_res^2-4*a11_res*c11_res
        P1_res[p1]=(-b11_res+delta11_res^(1/2))/(2*a11_res)
      }
      a2_res=eta_p*(KM+P1_res[Path1_size-1])
      b2_res=Vm2r*(P1_res[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_res[Path1_size-1]+KM)-Vm1r*P1_res[Path1_size-1]
      c2_res=-Vm1r*KM*P1_res[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_res[Path1_size-1]+KM)*KM
      delta2_res<-b2_res^2-4*a2_res*c2_res
      P2_res<-c()
      P2_res[1]=(-b2_res+delta2_res^(1/2))/(2*a2_res)
      for (p2 in 2:(Path_size-Path1_size)){
        a22_res=eta*(KM+P2_res[p2-1])
        b22_res=Vm2r*(P2_res[p2-1]+KM)+eta*KM*(P2_res[p2-1]+KM)-Vm2r*P2_res[p2-1]
        c22_res=-Vm2r*KM*P2_res[p2-1]
        delta22_res<-b22_res^2-4*a22_res*c22_res
        P2_res[p2]=(-b22_res+delta22_res^(1/2))/(2*a22_res)
      }
      Phi1r=0
      Phi2r=0
      Phi1r=Phi1r+(Vm0*Sin_res/(Sin_res+KM))/(Path_size/2)
      for (p1 in 1:(Path1_size-1)){
        Phi1r=Phi1r+(Vm1r*P1_res[p1]/(P1_res[p1]+KM))/(Path_size/2)
      }
      for (p2 in 1:(Path_size-Path1_size-1)){
        Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
      }
      Phi_res<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))
      for (j in 1:N_reso){
        E_tot_conc0=Etot_conc*f_var[j]
        E_tot_conc1<-Etot_conc*f_var[j]
        E_tot_conc2<-Etot_conc*(1-f_var[j])
        Etot_conc_m=E_tot_conc1*(Path1_size-1)+E_tot_conc2*(Path_size-Path1_size)+E_tot_conc0
        kf_act<-kf*10^(-(Etot_conc_m+Etot_back)/(E_basal))
        a<- (VTm*kf_act+kf_act*E_tot_conc0*kcat*(1+Se/KT))
        b<- (VTm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_tot_conc0*(KT+Se))
        c<- -VTm*Se*(kr+kcat)
        delta<-b^2-4*a*c
        #First cellular concentration (substrate)
        if(a==0){
          S_conc_eq=0
        }else{
          S_conc_eq=(-b+delta^(1/2))/(2*a)
        }
        KM=(kcat+kr)/kf_act
        Vm0=kcat*E_tot_conc0
        Vm1=kcat*E_tot_conc1
        Vm2=kcat*E_tot_conc2
        a_1=eta*(KM+S_conc_eq)
        b_1=Vm1*(S_conc_eq+KM)+eta*KM*(S_conc_eq+KM)-Vm0*S_conc_eq
        c_1=-Vm0*KM*S_conc_eq
        delta_1<-b_1^2-4*a_1*c_1
        P1_conc_eq<-c()
        P1_conc_eq[1]=(-b_1+delta_1^(1/2))/(2*a_1)
        for (p1 in 2:(Path1_size-1)){
          a_11=eta*(KM+P1_conc_eq[p1-1])
          b_11=Vm1*(P1_conc_eq[p1-1]+KM)+eta*KM*(P1_conc_eq[p1-1]+KM)-Vm1*P1_conc_eq[p1-1]
          c_11=-Vm1*KM*P1_conc_eq[p1-1]
          delta_11<-b_11^2-4*a_11*c_11
          P1_conc_eq[p1]=(-b_11+delta_11^(1/2))/(2*a_11)
        }
        
        a_2=eta_p*(KM+P1_conc_eq[Path1_size-1])
        b_2=Vm2*(P1_conc_eq[Path1_size-1]+KM)+(eta_p*KM-P_p2*S_c/V_c*P2_out)*(P1_conc_eq[Path1_size-1]+KM)-Vm1*P1_conc_eq[Path1_size-1]
        c_2=-Vm1*KM*P1_conc_eq[Path1_size-1]-P_p2*S_c/V_c*P2_out*(P1_conc_eq[Path1_size-1]+KM)*KM
        delta_2<-b_2^2-4*a_2*c_2
        P2_conc_eq<-c()
        P2_conc_eq[1]=(-b_2+delta_2^(1/2))/(2*a_2)
        
        for (p2 in 2:(Path_size-Path1_size)){
          a_22=eta*(KM+P2_conc_eq[p2-1])
          b_22=Vm2*(P2_conc_eq[p2-1]+KM)+eta*KM*(P2_conc_eq[p2-1]+KM)-Vm2*P2_conc_eq[p2-1]
          c_22=-Vm2*KM*P2_conc_eq[p2-1]
          delta_22<-b_22^2-4*a_22*c_22
          P2_conc_eq[p2]=(-b_22+delta_22^(1/2))/(2*a_22)
        }
        Phi1=0
        Phi2=0
        Phi1=Phi1+(Vm0*S_conc_eq/(S_conc_eq+KM))/(Path_size/2)
        for (p1 in 1:(Path1_size-1)){
          Phi1=Phi1+(Vm1*P1_conc_eq[p1]/(P1_conc_eq[p1]+KM))/(Path_size/2)
        }
        for (p2 in 1:(Path_size-Path1_size-1)){
          Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
        }
        
        if(as.numeric((Phi1*PO1+Phi2*PO2)-prot_cost*((E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size))+(Etot_back)))>0){
          P_Phi_eq=as.numeric(((Phi1*PO1+Phi2*PO2)-prot_cost*(Etot_conc_m+Etot_back))*Tox/(Tox+sum(P1_conc_eq)+sum(P2_conc_eq)))
        }
        else{
          P_Phi_eq=0
        }
        if(N_eq>1){
          fit=P_Phi_eq/Phi_res-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[i,j]=P_Phi_eq
        fit_inv[i,j]=fit
      }
    }
    list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]]<-fit_inv
  }
}

ncol=800
jet.colors <- colorRampPalette(c(rep("white",400),rep("grey30",400)))
palet<-jet.colors(ncol)
pal<-list(palet)
Pop_DA1<-Pop_DA
multiplePlot("Main=","","Deg. rate=","",c("1st","Both","2nd"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5")
             ,ncol=800,f_var_print,f_var_print,list_fit_inv_S40ratvar_Y1_eta1,
             abs="res. (% 1st pathway)",ord="mut. (% 1st pathway)",scale=c(-80,80),lev=c(0),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


dev.print(device = jpeg, file = "res-ToxDegTrans-2P-E0var25.jpeg", , width = 2100*3,height=800*3,res=300,type="cairo")
dev.off()

Rep_opt_DA_Deg_S40_eta_sens<-list()
for (po in 1:length(Pox)){
  Rep_opt_DA_Deg_S40_eta_sens[[po]]<-c(0)
  for (p in 1:length(eta_set)){
    #print(p)
    i=1
    while(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      #print(i)
      if(i>1 && i<(N_reso-1)){
        if(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i+1,i]>0 && list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i-1,i]>0){
          print(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i])
          S_opt=i
          print(S_opt)
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]<-f_var[S_opt]
          break
        }
      }
      else if (i==(N_reso-1)){
        if(list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i+1,i]<0 && list_fit_inv_S40ratvar_Y1_eta1[[p+length(eta_set)*(po-1)]][i,i+1]>0){
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]=1
          i=i+1
          break
        }
        else{
          Rep_opt_DA_Deg_S40_eta_sens[[po]][p]<-0
          break
        }
      }
      i=i+1
    }
  }
}

col_set<-c("darkgrey","pink","palegreen")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(po in 1:length(Pox)){
  Pox_i<-as.character(Pox[po])
  if (po==1){
    plot(Rep_opt_DA_Deg_S40_eta_sens[[1]]~log10(eta_set),col=col_set[po],pch=16,ylim=c(0.25,1),xlab="log10(Degradation rate)",ylab="(% invetsment in 1st pathway)",cex=2.25)
  }
  else{
    points(Rep_opt_DA_Deg_S40_eta_sens[[po]]~log10(eta_set),col=col_set[po],pch=16,cex=2.25)
  }
}
for(po in 1:length(Pox)){
  Pox_i<-as.character(Pox[po])
  points(Rep_opt_DA_Deg_S40_eta_sens[[po]]~log10(eta_set),col=po,pch=16,cex=1)
}
par(xpd=F)
text(-4,1.01,"B",srt=addtxt$srt,font=2,col="black",cex=1)
abline(h=0.5,col="gray 10",lty=4)
text(-1.6,0.54,"50%",srt=addtxt$srt,font=1,col="black",cex=1)
legend("bottomright",title="Pay-off",legend=c("10:1","1:1","1:10"),bty="y",pch=c(16),col=c(1:3),ncol=1,cex=0.85)
legend("bottomleft",title="Singular strat.",legend=c("Loc.stable","Glob.stab."),bty="y",pch=c(16,16),col="gray20",ncol=1,pt.cex=c(0.85,1.5))

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-ToxDegTrans-2P-EOvar-gen25.jpeg", width = 560*3,height=475*3,res=300,type="cairo")
dev.off()


delta_1<-function(Ks){
  del_opt<-2*(1+1/Ks)-sqrt((2*Ks+3)*(Ks+1))/Ks
  if(del_opt>1){
    del_opt<-1
  }
  return(del_opt)
}
delta_1(10^2)
delta_res<-c()
log10_KS_set<-seq(-3,3,0.01)
for (i in 1:length(log10_KS_set)){
  delta_res[i]=delta_1(10^log10_KS_set[i])
}
delta_res

par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875),xpd=FALSE)
plot(log10_KS_set,delta_res,type="l",col="red",lwd=2,ylim=c(0.5,1),xlab=expression(paste("log10 (",theta[s],")")),ylab=expression(delta),cex.lab=1.25)
abline(h=2-sqrt(2),lty=2)
text(-0,0.555,"Minimum overexpression
     in first subpathway",srt=addtxt$srt,font=1,col="black",cex=1)
par(xpd=TRUE)
text(-3,0.39,"Unsaturated
pathway",srt=addtxt$srt,font=1,col="grey20",cex=1)
text(3,0.39,"Saturated
pathway",srt=addtxt$srt,font=1,col="grey20",cex=1)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-Sat-ToyModel.jpeg", width = 560*3,height=475*3,res=300,type="cairo")
dev.off()

