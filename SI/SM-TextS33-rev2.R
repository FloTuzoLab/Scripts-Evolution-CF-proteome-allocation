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

##Preamble
options(digits=20)#To avoid problems with rounding numbers
#Defining features of plots
ncol=256
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato"))
palet<-jet.colors(ncol)
pal<-list(palet)

#3a.Adaptive dynamics with 2 sub-pathways, first pay-off value, high toxicity

Tox_set=c(10^-0.5,10^-1.5)
VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.25
kcat=10^2.25#/s
kr_set=c(kcat/9,kcat/3,kcat)#/s
kinh_set=c(kf,kf/3,kf/9)
kinh_spe_set=c(kf,kf/3,kf/9)
prot_cost=10^-2.5#Mmet/Menz
eta_set=c(10^-3)
PayOff_set<-c(1,0.1)
Path_size_set=c(40)
N_reso=50
Etot_back<-5.5*10^-3#other protein concentration
E_basal<-3*10^-3#scaling factor

P_p2_set<-c(0)

#Adaptive dynamics constants
Phi_eq=10^-4
r_c=1e-5#dm
S_c=4*pi*r_c^2
V_c=4/3*pi*r_c^3#dm^3
beta=10^-3
alpha=10^-3
dl=1e-4#dm
V_env=dl^3

#Variables
Etot_var<-10^seq(-5,-4,length=N_reso)
#Outcome
tab_N_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_S_eq <- c()
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_P2<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)

Neq<-list()
Phi_c<-list()
S_eq<-list()
Sin_eq<-list()
P1_eq<-list()
P2_eq<-list()
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40_rev_noperm<-list()

#Parameter set 1
Tox=Tox_set[2]
Path_size=Path_size_set[1]
Path1_size=Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)
eta<-eta_set[1]

init<-c()
for (i0 in 1:(Path_size+3)){
  init[i0]<-1e-5
}
initm<-c()
for (i0 in 1:(Path_size+1)){
  initm[i0]<-1e-5
}
Enz=1
for(re in 1:length(kr_set)){
  kr=kr_set[re]
  kinh=kinh_set[re]
  kinh_spe=kinh_spe_set[re]
  for (po in 1:length(PayOff_set)){
    #P_p2=P_p2_set[p]
    PayOff<-PayOff_set[po]
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    for (k in 1:N_reso){
      
      print(paste(re,po,k))
      Etot1r=Etot_var[k]#Molar
      Etot2r=Etot_var[k]#Molar
      #Intermediate constant
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kf_act=kf*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act=kinh*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act_spe=kinh_spe*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      KI=kinh_act/kf_act
      KI_spe=kinh_act_spe/kf_act
      KM1_f=(kr+kcat)/kf_act
      KM1_b=(kr+kcat)/kinh_act
      KM1_b_spe=(kr+kcat)/kinh_act_spe
      
      count=0
      N_eq_p=0
      N_eq=1
      Phi_fit=0
      #print(KM1_f)
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-6) && N_eq>0.99){
        metaboChain<-function(x){
          e1<-c()
          e2<-c()
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-(eta*x[2]+Vm1r*x[2]/(KM1_f+x[2]+KI*x[3])-Vm1r*(x[3])/(KM1_b+x[3]+x[2]/KI))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1r*x[i+1]/(KM1_f+x[i+1]+KI*x[i+2])-Vm1r*(x[i+2])/(KM1_b+x[i+2]+x[i+1]/KI))-(eta*x[i+2]+Vm1r*x[i+2]/(KM1_f+x[i+2]+KI*x[i+3])-Vm1r*(x[i+3])/(KM1_b+x[i+3]+x[i+2]/KI))
          }
          e1[Path1_size] <- (Vm1r*x[Path1_size+1]/(KM1_f+x[Path1_size+1]+KI*x[Path1_size+2])-Vm1r*(x[Path1_size+2])/(KM1_b+x[Path1_size+2]+x[Path1_size+1]/KI))-(eta*x[Path1_size+2]+P_p2*S_c/V_c*(x[Path1_size+2]-x[Path_size+3])+Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-(Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe)))
          e2[1] <- (Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe))-(eta*x[Path1_size+3]+Vm2r*x[Path1_size+3]/(KM1_f+x[Path1_size+3]+KI*x[Path1_size+4])-Vm2r*(x[Path1_size+4])/(KM1_b+x[Path1_size+4]+x[Path1_size+3]/KI))
          t3<-N_eq*P_p2*S_c/V_env*(x[Path1_size+2]-x[Path_size+3])-beta*x[Path_size+3]
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2r*x[Path1_size+2+j]/(KM1_f+x[Path1_size+2+j]+KI*x[Path1_size+3+j])-Vm2r*(x[Path1_size+3+j])/(KM1_b+x[Path1_size+3+j]+x[Path1_size+2+j]/KI))-(eta*x[Path1_size+3+j]+Vm2r*x[Path1_size+3+j]/(KM1_f+x[Path1_size+3+j]+KI*x[Path1_size+4+j])-Vm2r*(x[Path1_size+4+j])/(KM1_b+x[Path1_size+4+j]+x[Path1_size+3+j]/KI))
          }
          e2[Path_size-Path1_size]<-(Vm2r*x[Path_size+1]/(KM1_f+x[Path_size+1]+KI*x[Path_size+2])-Vm2r*(x[Path_size+2])/(KM1_b+x[Path_size+2]+x[Path_size+1]/KI))-(eta*x[Path_size+2]+Vm2r*x[Path_size+2]/(KM1_f+x[Path_size+2]))
          c(t1,t2,e1,t3,e2)
        }
        Reso<-nleqslv(init,metaboChain,jac=NULL)$x
        Sout<-Reso[1]
        Sin<-Reso[2]
        P1<-Reso[3:(Path_size+2)]
        P2_out<-Reso[Path_size+3]
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm1r*Sin/(KM1_f+Sin+KI*P1[1])-Vm1r*(P1[1])/(KM1_b+P1[1]+Sin/KI))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(KM1_f+P1[p1]+KI*P1[p1+1])-Vm1r*(P1[p1+1])/(KM1_b+P1[p1+1]+P1[p1]/KI))/((Path_size/2))
        }
        Phi2r=Phi2r+(Vm2r*P1[Path1_size]/(KM1_f+P1[Path1_size]+KI_spe*P1[Path1_size+1])-Vm2r*(P1[Path1_size+1])/(KM1_b_spe+P1[Path1_size+1]+P1[Path1_size]/KI_spe))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2r=Phi2r+(Vm2r*P1[p2]/(KM1_f+P1[p2]+KI*P1[p2+1])-Vm2r*(P1[p2+1])/(KM1_b+P1[p2+1]+P1[p2]/KI))/((Path_size/2))
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1))
        N_eq_p=N_eq
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq-2)
          }
        }
        else if(count<100){
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
      Pop_DA[k]<-N_eq
      P1_eq[[k]]<-P1
      P2_eq[[k]]<-P2_out
      for (l in 1:N_reso){
        Etot1m=Etot_var[l]#Molar
        Etot2m=Etot_var[l]#Molar
        Vm1m=kcat*Etot1m
        Vm2m=kcat*Etot2m
        kf_act_m=kf*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_m=kinh*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_spe_m=kinh_spe*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        KI_m=kinh_act_m/kf_act_m
        KI_spe_m=kinh_act_spe_m/kf_act_m##note: previously used to test the effect of a different reversibility level in the pathway
        KM1_f_m=(kr+kcat)/kf_act_m
        KM1_b_m=(kr+kcat)/kinh_act_m
        KM1_b_spe_m=(kr+kcat)/kinh_act_spe_m
        metaboChainm<-function(y){
          e1<-c()
          e2<-c()
          t2 <- VTm*(Sout-y[1])/(KT+(Sout+y[1])+Sout*y[1]/KT)-(eta*y[1]+Vm1m*y[1]/(KM1_f_m+y[1]+KI_m*y[2])-Vm1m*(y[2])/(KM1_b_m+y[2]+y[1]/KI_m))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1m*y[i]/(KM1_f_m+y[i]+KI_m*y[i+1])-Vm1m*(y[i+1])/(KM1_b_m+y[i+1]+y[i]/KI_m))-(eta*y[i+1]+Vm1m*y[i+1]/(KM1_f_m+y[i+1]+KI_m*y[i+2])-Vm1m*(y[i+2])/(KM1_b_m+y[i+2]+y[i+1]/KI_m))
          }
          e1[Path1_size] <- (Vm1m*y[Path1_size]/(KM1_f_m+y[Path1_size]+KI_m*y[Path1_size+1])-Vm1m*(y[Path1_size+1])/(KM1_b_m+y[Path1_size+1]+y[Path1_size]/KI_m))-(eta*y[Path1_size+1]+P_p2*S_c/V_c*(y[Path1_size+1]-P2_out)+Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-(Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m)))
          e2[1] <- (Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m))-(eta*y[Path1_size+2]+Vm2m*y[Path1_size+2]/(KM1_f_m+y[Path1_size+2]+KI_m*y[Path1_size+3])-Vm2m*(y[Path1_size+3])/(KM1_b_m+y[Path1_size+3]+y[Path1_size+2]/KI_m))
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2m*y[Path1_size+1+j]/(KM1_f_m+y[Path1_size+1+j]+KI_m*y[Path1_size+2+j])-Vm2m*(y[Path1_size+2+j])/(KM1_b_m+y[Path1_size+2+j]+y[Path1_size+1+j]/KI_m))-(eta*y[Path1_size+2+j]+Vm2m*y[Path1_size+2+j]/(KM1_f_m+y[Path1_size+2+j]+KI_m*y[Path1_size+3+j])-Vm2m*(y[Path1_size+3+j])/(KM1_b_m+y[Path1_size+3+j]+y[Path1_size+2+j]/KI_m))
          }
          e2[Path_size-Path1_size]<-(Vm2m*y[Path_size]/(KM1_f_m+y[Path_size]+KI_m*y[Path_size+1])-Vm2m*(y[Path_size+1])/(KM1_b_m+y[Path_size+1]+y[Path_size]/KI_m))-(eta*y[Path_size+1]+Vm2m*y[Path_size+1]/(KM1_f_m+y[Path_size+1]))
          c(t2,e1,e2)
        }
        Reso_m<-nleqslv(initm,metaboChainm,jac=NULL)$x
        Sin_m<-Reso_m[1]
        P1_m<-Reso_m[2:(Path_size+1)]
        Phi1m=0
        Phi2m=0
        Phi1m=Phi1m+(Vm1m*Sin_m/(KM1_f_m+Sin_m+KI_m*P1_m[1])-Vm1m*(P1_m[1])/(KM1_b_m+P1_m[1]+Sin_m/KI_m))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1m=Phi1m+(Vm1m*P1_m[p1]/(KM1_f_m+P1_m[p1]+KI_m*P1_m[p1+1])-Vm1m*(P1_m[p1+1])/(KM1_b_m+P1_m[p1+1]+P1_m[p1]/KI_m))/((Path_size/2))
        }
        Phi2m=Phi2m+(Vm2m*P1_m[Path1_size]/(KM1_f_m+P1_m[Path1_size]+KI_spe_m*P1_m[Path1_size+1])-Vm2m*(P1_m[Path1_size+1])/(KM1_b_spe_m+P1_m[Path1_size+1]+P1_m[Path1_size]/KI_spe_m))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2m=Phi2m+(Vm2m*P1_m[p2]/(KM1_f_m+P1_m[p2]+KI_m*P1_m[p2+1])-Vm2m*(P1_m[p2+1])/(KM1_b_m+P1_m[p2+1]+P1_m[p2]/KI_m))/((Path_size/2))
        }
        Phi_fit_m<-((Phi1m*PO1+Phi2m*PO2)-prot_cost*(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1_m))
        if(Phi_fit_m>0){}else{
          Phi_fit_m=0
        }
        if(N_eq>1){
          fit=Phi_fit_m/Phi_fit-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[k,l]=Phi_fit_m
        fit_inv[k,l]=fit
      }
    }
    list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]]<-fit_inv
  }
}
for (i in 1:N_reso){
  print(fit_inv[i,i])
}
ncol=800
jet.colors <- colorRampPalette(c(rep("white",400),rep("grey30",400)))
palet<-jet.colors(ncol)
pal<-list(palet)

addtxt<-list(l=0.95,h=0.05,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
log10Etot1_var<-c(-6,-4,0.5)
multiplePlot("","","","",c("","",""),c("","")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_rev_noperm,
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


Conc_opt_DA_Degperm_S40_rev_noperm<-list()
Reversibility<-c("kr","both","kinh")
for(re in 1:length(kr_set)){
  Rev_i<-as.character(Reversibility[re])
  print(Rev_i)
  Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]]<-c()
  for (po in 1:length(PayOff_set)){
    #print(p)
    i=1
    while(list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      if(list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i,i+2]<0 && list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i,i-2]<0){
        print(list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][po]<-Etot_var[S_opt]
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i+1,i]<0 && list_fit_inv_S40_rev_noperm[[po+length(PayOff_set)*(re-1)]][i,i+1]>0){
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][po]=1
          i=i+1
          break
        }
        else{
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][po]<-0
          break
        }
      }
      else{
        i=i+1
      }
    }
  }
}

#Constant
#Tox<-Tox_set[2]
PayOff<-PayOff_set[1]

#Parameter of interest
P_p2_set<-10^seq(-9,-4,1)

##Redefining space of adaptive parameters to coincide with an allocation
N_reso=25
f_var<-seq(0,1,length=N_reso)

##New outcome list with new resolution
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40_Yvar_Rev<-list()

for(re in 1:length(kr_set)){
  list_fit_inv_S40_Yvar_Rev[[re]]<-list()
  Rev_i<-as.character(Reversibility[re])
  kr=kr_set[re]
  kinh=kinh_set[re]
  kinh_spe=kinh_spe_set[re]
  for (pe in 1:length(P_p2_set)){
    P_p2=P_p2_set[pe]
    #PayOff<-PayOff_set[po]
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    for (k in 1:N_reso){
      Etot_conc=Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][1]*2
      print(paste(re,pe,k))
      Etot1r=Etot_conc*f_var[k]#Molar
      Etot2r=Etot_conc*(1-f_var[k])#Molar
      #Intermediate constant
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kf_act=kf*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act=kinh*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act_spe=kinh_spe*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      KI=kinh_act/kf_act
      KI_spe=kinh_act_spe/kf_act
      KM1_f=(kr+kcat)/kf_act
      KM1_b=(kr+kcat)/kinh_act
      KM1_b_spe=(kr+kcat)/kinh_act_spe
      count=0
      N_eq_p=0
      N_eq=1
      Phi_fit=0
      #print(KM1_f)
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-6) && N_eq>0.99){
        metaboChain<-function(x){
          e1<-c()
          e2<-c()
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-(eta*x[2]+Vm1r*x[2]/(KM1_f+x[2]+KI*x[3])-Vm1r*(x[3])/(KM1_b+x[3]+x[2]/KI))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1r*x[i+1]/(KM1_f+x[i+1]+KI*x[i+2])-Vm1r*(x[i+2])/(KM1_b+x[i+2]+x[i+1]/KI))-(eta*x[i+2]+Vm1r*x[i+2]/(KM1_f+x[i+2]+KI*x[i+3])-Vm1r*(x[i+3])/(KM1_b+x[i+3]+x[i+2]/KI))
          }
          e1[Path1_size] <- (Vm1r*x[Path1_size+1]/(KM1_f+x[Path1_size+1]+KI*x[Path1_size+2])-Vm1r*(x[Path1_size+2])/(KM1_b+x[Path1_size+2]+x[Path1_size+1]/KI))-(eta*x[Path1_size+2]+P_p2*S_c/V_c*(x[Path1_size+2]-x[Path_size+3])+Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-(Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe)))
          e2[1] <- (Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe))-(eta*x[Path1_size+3]+Vm2r*x[Path1_size+3]/(KM1_f+x[Path1_size+3]+KI*x[Path1_size+4])-Vm2r*(x[Path1_size+4])/(KM1_b+x[Path1_size+4]+x[Path1_size+3]/KI))
          t3<-N_eq*P_p2*S_c/V_env*(x[Path1_size+2]-x[Path_size+3])-beta*x[Path_size+3]
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2r*x[Path1_size+2+j]/(KM1_f+x[Path1_size+2+j]+KI*x[Path1_size+3+j])-Vm2r*(x[Path1_size+3+j])/(KM1_b+x[Path1_size+3+j]+x[Path1_size+2+j]/KI))-(eta*x[Path1_size+3+j]+Vm2r*x[Path1_size+3+j]/(KM1_f+x[Path1_size+3+j]+KI*x[Path1_size+4+j])-Vm2r*(x[Path1_size+4+j])/(KM1_b+x[Path1_size+4+j]+x[Path1_size+3+j]/KI))
          }
          e2[Path_size-Path1_size]<-(Vm2r*x[Path_size+1]/(KM1_f+x[Path_size+1]+KI*x[Path_size+2])-Vm2r*(x[Path_size+2])/(KM1_b+x[Path_size+2]+x[Path_size+1]/KI))-(eta*x[Path_size+2]+Vm2r*x[Path_size+2]/(KM1_f+x[Path_size+2]))
          c(t1,t2,e1,t3,e2)
        }
        Reso<-nleqslv(init,metaboChain,jac=NULL)$x
        Sout<-Reso[1]
        Sin<-Reso[2]
        P1<-Reso[3:(Path_size+2)]
        P2_out<-Reso[Path_size+3]
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm1r*Sin/(KM1_f+Sin+KI*P1[1])-Vm1r*(P1[1])/(KM1_b+P1[1]+Sin/KI))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(KM1_f+P1[p1]+KI*P1[p1+1])-Vm1r*(P1[p1+1])/(KM1_b+P1[p1+1]+P1[p1]/KI))/((Path_size/2))
        }
        Phi2r=Phi2r+(Vm2r*P1[Path1_size]/(KM1_f+P1[Path1_size]+KI_spe*P1[Path1_size+1])-Vm2r*(P1[Path1_size+1])/(KM1_b_spe+P1[Path1_size+1]+P1[Path1_size]/KI_spe))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2r=Phi2r+(Vm2r*P1[p2]/(KM1_f+P1[p2]+KI*P1[p2+1])-Vm2r*(P1[p2+1])/(KM1_b+P1[p2+1]+P1[p2]/KI))/((Path_size/2))
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1)+Sin)
        N_eq_p=N_eq
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq-2)
          }
        }
        else if(count<100){
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
      Pop_DA[k]<-N_eq
      P1_eq[[k]]<-P1
      P2_eq[[k]]<-P2_out
      for (l in 1:N_reso){
        Etot1m=Etot_conc*f_var[l]#Molar
        Etot2m=Etot_conc*(1-f_var[l])#Molar
        Vm1m=kcat*Etot1m
        Vm2m=kcat*Etot2m
        kf_act_m=kf*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_m=kinh*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_spe_m=kinh_spe*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        KI_m=kinh_act_m/kf_act_m
        KI_spe_m=kinh_act_spe_m/kf_act_m##note: previously used to test the effect of a different reversibility level in the pathway
        KM1_f_m=(kr+kcat)/kf_act_m
        KM1_b_m=(kr+kcat)/kinh_act_m
        KM1_b_spe_m=(kr+kcat)/kinh_act_spe_m
        metaboChainm<-function(y){
          e1<-c()
          e2<-c()
          t2 <- VTm*(Sout-y[1])/(KT+(Sout+y[1])+Sout*y[1]/KT)-(eta*y[1]+Vm1m*y[1]/(KM1_f_m+y[1]+KI_m*y[2])-Vm1m*(y[2])/(KM1_b_m+y[2]+y[1]/KI_m))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1m*y[i]/(KM1_f_m+y[i]+KI_m*y[i+1])-Vm1m*(y[i+1])/(KM1_b_m+y[i+1]+y[i]/KI_m))-(eta*y[i+1]+Vm1m*y[i+1]/(KM1_f_m+y[i+1]+KI_m*y[i+2])-Vm1m*(y[i+2])/(KM1_b_m+y[i+2]+y[i+1]/KI_m))
          }
          e1[Path1_size] <- (Vm1m*y[Path1_size]/(KM1_f_m+y[Path1_size]+KI_m*y[Path1_size+1])-Vm1m*(y[Path1_size+1])/(KM1_b_m+y[Path1_size+1]+y[Path1_size]/KI_m))-(eta*y[Path1_size+1]+P_p2*S_c/V_c*(y[Path1_size+1]-P2_out)+Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-(Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m)))
          e2[1] <- (Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m))-(eta*y[Path1_size+2]+Vm2m*y[Path1_size+2]/(KM1_f_m+y[Path1_size+2]+KI_m*y[Path1_size+3])-Vm2m*(y[Path1_size+3])/(KM1_b_m+y[Path1_size+3]+y[Path1_size+2]/KI_m))
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2m*y[Path1_size+1+j]/(KM1_f_m+y[Path1_size+1+j]+KI_m*y[Path1_size+2+j])-Vm2m*(y[Path1_size+2+j])/(KM1_b_m+y[Path1_size+2+j]+y[Path1_size+1+j]/KI_m))-(eta*y[Path1_size+2+j]+Vm2m*y[Path1_size+2+j]/(KM1_f_m+y[Path1_size+2+j]+KI_m*y[Path1_size+3+j])-Vm2m*(y[Path1_size+3+j])/(KM1_b_m+y[Path1_size+3+j]+y[Path1_size+2+j]/KI_m))
          }
          e2[Path_size-Path1_size]<-(Vm2m*y[Path_size]/(KM1_f_m+y[Path_size]+KI_m*y[Path_size+1])-Vm2m*(y[Path_size+1])/(KM1_b_m+y[Path_size+1]+y[Path_size]/KI_m))-(eta*y[Path_size+1]+Vm2m*y[Path_size+1]/(KM1_f_m+y[Path_size+1]))
          c(t2,e1,e2)
        }
        Reso_m<-nleqslv(initm,metaboChainm,jac=NULL)$x
        Sin_m<-Reso_m[1]
        P1_m<-Reso_m[2:(Path_size+1)]
        Phi1m=0
        Phi2m=0
        Phi1m=Phi1m+(Vm1m*Sin_m/(KM1_f_m+Sin_m+KI_m*P1_m[1])-Vm1m*(P1_m[1])/(KM1_b_m+P1_m[1]+Sin_m/KI_m))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1m=Phi1m+(Vm1m*P1_m[p1]/(KM1_f_m+P1_m[p1]+KI_m*P1_m[p1+1])-Vm1m*(P1_m[p1+1])/(KM1_b_m+P1_m[p1+1]+P1_m[p1]/KI_m))/((Path_size/2))
        }
        Phi2m=Phi2m+(Vm2m*P1_m[Path1_size]/(KM1_f_m+P1_m[Path1_size]+KI_spe_m*P1_m[Path1_size+1])-Vm2m*(P1_m[Path1_size+1])/(KM1_b_spe_m+P1_m[Path1_size+1]+P1_m[Path1_size]/KI_spe_m))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2m=Phi2m+(Vm2m*P1_m[p2]/(KM1_f_m+P1_m[p2]+KI_m*P1_m[p2+1])-Vm2m*(P1_m[p2+1])/(KM1_b_m+P1_m[p2+1]+P1_m[p2]/KI_m))/((Path_size/2))
        }
        Phi_fit_m<-((Phi1m*PO1+Phi2m*PO2)-prot_cost*(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1_m)+Sin_m)
        if(Phi_fit_m>0){}else{
          Phi_fit_m=0
        }
        if(N_eq>1){
          fit=Phi_fit_m/Phi_fit-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[k,l]=Phi_fit_m
        fit_inv[k,l]=fit
      }
    }
    list_fit_inv_S40_Yvar_Rev[[re]][[pe]]<-fit_inv
  }
}


multiplePlot("Main=","","Deg. rate=","",c("1:1","1:2"),c("1e-4","1e-3","&")
             ,ncol=800,f_var_print,f_var_print,list_fit_inv_S40_Yvar_Rev[[1]],
             abs="res. (% 1st pathway)",ord="mut. (% 1st pathway)",scale=c(-99,99),lev=c(0),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")
Rep_opt_DA_RevTox_S40_perm<-list()
for (re in 1:length(Reversibility)){
  Rev_i<-as.character(Reversibility[re])
  Rep_opt_DA_RevTox_S40_perm[[re]]<-c(0)
  for (pe in 1:length(P_p2_set)){
    i=1
    while(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      #print(i)
      if(i>1 && i<(N_reso-1)){
        if(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i+1,i]>0 && list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i-1,i]>0){
          print(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i])
          S_opt=i
          print(S_opt)
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]<-f_var[S_opt]
          break
        }
      }
      else if (i==(N_reso-1)){
        if(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i+1,i]<0 && list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i+1]>0){
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]=1
          i=i+1
          break
        }
        else{
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]<-0
          break
        }
      }
      i=i+1
    }
  }
}

Rep_opt_DA_RevTox_S40_perm[[1]][1]<-Rep_opt_DA_RevTox_S40_perm[[1]][2]#Manual identification

##Main plot about the adaptive dynamics outcomes with permeability
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1,0.5,1))
plot(Rep_opt_DA_RevTox_S40_perm[[1]]~log10(P_p2_set),cex=0.25,col=1,pch=15,ylim=c(-0,1),xlab="log10(Permeability rate)",ylab=expression(paste(delta," (% 1st pathway)")))
points(Rep_opt_DA_RevTox_S40_perm[[1]][1:3]~log10(P_p2_set)[1:3],col="light green",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[1]][1:3]~log10(P_p2_set)[1:3],col="dark green",pch=16,cex=1)
points(Rep_opt_DA_RevTox_S40_perm[[1]][4:6]~log10(P_p2_set)[4:6],col="light green",pch=1,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[1]][4:6]~log10(P_p2_set)[4:6],col="dark green",pch=16,cex=1)
points(Rep_opt_DA_RevTox_S40_perm[[2]]~log10(P_p2_set),col="light blue",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[2]]~log10(P_p2_set),col="blue",pch=16)
points(Rep_opt_DA_RevTox_S40_perm[[3]]~log10(P_p2_set),col="lightsalmon",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[3]]~log10(P_p2_set),col="brown",pch=16)
text(-9,1,"C",srt=addtxt$srt,font=2,col="black",cex=1)
#text(-3.8,0.18,"*",srt=addtxt$srt,font=2,col="black",cex=1)
#text(-3.8,0.1,"o",srt=addtxt$srt,font=2,col=" light blue",cex=1)
par(xpd=TRUE)
legend("bottomleft",title="Singular strat.",legend=c("Loc. stable","Loc. instable","Glob. stable","Glob. instable"),bty="y",pch=c(20,1,19,1),col=c(1,"dark grey",1,1),pt.cex=c(1,1,1.75,1.75),ncol=1,cex=0.9)
legend("bottomright",title="Reversibility spread",legend=c("mostly kinh","equal","mostly kr"),bty="y",pch=c(15),col=c("green","blue","red"),ncol=1,cex=0.9)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "Cross-feeding-PO1-toxhighrev.jpeg", width = 575*6,height=500*6,res=600,type="cairo")
dev.off()

#2b.Adaptive dynamics with 2 sub-pathways, first pay-off value, low toxicity

#Constant
PayOff<-PayOff_set[2]

#Parameter of interest
P_p2_set<-10^seq(-9,-4,1)

##Redefining space of adaptive parameters to coincide with an allocation
N_reso=25
f_var<-seq(0,1,length=N_reso)

##New outcome list with new resolution
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_S40_Yvar_Rev<-list()

for(re in 1:length(kr_set)){
  list_fit_inv_S40_Yvar_Rev[[re]]<-list()
  Rev_i<-as.character(Reversibility[re])
  kr=kr_set[re]
  kinh=kinh_set[re]
  kinh_spe=kinh_spe_set[re]
  for (pe in 1:length(P_p2_set)){
    P_p2=P_p2_set[pe]
    #PayOff<-PayOff_set[po]
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    for (k in 1:N_reso){
      Etot_conc=Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][2]*2
      print(paste(re,pe,k))
      Etot1r=Etot_conc*f_var[k]#Molar
      Etot2r=Etot_conc*(1-f_var[k])#Molar
      #Intermediate constant
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kf_act=kf*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act=kinh*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act_spe=kinh_spe*10^(-(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      KI=kinh_act/kf_act
      KI_spe=kinh_act_spe/kf_act
      KM1_f=(kr+kcat)/kf_act
      KM1_b=(kr+kcat)/kinh_act
      KM1_b_spe=(kr+kcat)/kinh_act_spe
      count=0
      N_eq_p=0
      N_eq=1
      Phi_fit=0
      #print(KM1_f)
      while((abs(Phi_fit-Phi_eq)/Phi_eq>10^-6) && N_eq>0.99){
        metaboChain<-function(x){
          e1<-c()
          e2<-c()
          t1 <- alpha-beta*x[1]-V_c/V_env*N_eq*VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-(eta*x[2]+Vm1r*x[2]/(KM1_f+x[2]+KI*x[3])-Vm1r*(x[3])/(KM1_b+x[3]+x[2]/KI))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1r*x[i+1]/(KM1_f+x[i+1]+KI*x[i+2])-Vm1r*(x[i+2])/(KM1_b+x[i+2]+x[i+1]/KI))-(eta*x[i+2]+Vm1r*x[i+2]/(KM1_f+x[i+2]+KI*x[i+3])-Vm1r*(x[i+3])/(KM1_b+x[i+3]+x[i+2]/KI))
          }
          e1[Path1_size] <- (Vm1r*x[Path1_size+1]/(KM1_f+x[Path1_size+1]+KI*x[Path1_size+2])-Vm1r*(x[Path1_size+2])/(KM1_b+x[Path1_size+2]+x[Path1_size+1]/KI))-(eta*x[Path1_size+2]+P_p2*S_c/V_c*(x[Path1_size+2]-x[Path_size+3])+Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-(Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe)))
          e2[1] <- (Vm2r*x[Path1_size+2]/(KM1_f+x[Path1_size+2]+KI_spe*x[Path1_size+3])-Vm2r*(x[Path1_size+3])/(KM1_b_spe+x[Path1_size+3]+x[Path1_size+2]/KI_spe))-(eta*x[Path1_size+3]+Vm2r*x[Path1_size+3]/(KM1_f+x[Path1_size+3]+KI*x[Path1_size+4])-Vm2r*(x[Path1_size+4])/(KM1_b+x[Path1_size+4]+x[Path1_size+3]/KI))
          t3<-N_eq*P_p2*S_c/V_env*(x[Path1_size+2]-x[Path_size+3])-beta*x[Path_size+3]
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2r*x[Path1_size+2+j]/(KM1_f+x[Path1_size+2+j]+KI*x[Path1_size+3+j])-Vm2r*(x[Path1_size+3+j])/(KM1_b+x[Path1_size+3+j]+x[Path1_size+2+j]/KI))-(eta*x[Path1_size+3+j]+Vm2r*x[Path1_size+3+j]/(KM1_f+x[Path1_size+3+j]+KI*x[Path1_size+4+j])-Vm2r*(x[Path1_size+4+j])/(KM1_b+x[Path1_size+4+j]+x[Path1_size+3+j]/KI))
          }
          e2[Path_size-Path1_size]<-(Vm2r*x[Path_size+1]/(KM1_f+x[Path_size+1]+KI*x[Path_size+2])-Vm2r*(x[Path_size+2])/(KM1_b+x[Path_size+2]+x[Path_size+1]/KI))-(eta*x[Path_size+2]+Vm2r*x[Path_size+2]/(KM1_f+x[Path_size+2]))
          c(t1,t2,e1,t3,e2)
        }
        Reso<-nleqslv(init,metaboChain,jac=NULL)$x
        Sout<-Reso[1]
        Sin<-Reso[2]
        P1<-Reso[3:(Path_size+2)]
        P2_out<-Reso[Path_size+3]
        Phi1r=0
        Phi2r=0
        Phi1r=Phi1r+(Vm1r*Sin/(KM1_f+Sin+KI*P1[1])-Vm1r*(P1[1])/(KM1_b+P1[1]+Sin/KI))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1r=Phi1r+(Vm1r*P1[p1]/(KM1_f+P1[p1]+KI*P1[p1+1])-Vm1r*(P1[p1+1])/(KM1_b+P1[p1+1]+P1[p1]/KI))/((Path_size/2))
        }
        Phi2r=Phi2r+(Vm2r*P1[Path1_size]/(KM1_f+P1[Path1_size]+KI_spe*P1[Path1_size+1])-Vm2r*(P1[Path1_size+1])/(KM1_b_spe+P1[Path1_size+1]+P1[Path1_size]/KI_spe))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2r=Phi2r+(Vm2r*P1[p2]/(KM1_f+P1[p2]+KI*P1[p2+1])-Vm2r*(P1[p2+1])/(KM1_b+P1[p2+1]+P1[p2]/KI))/((Path_size/2))
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot1r*Path1_size+Etot2r*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1)+Sin)
        N_eq_p=N_eq
        if(count<10){
          if(Phi_fit>Phi_eq){
            N_eq=(N_eq*2)
          }
          else{
            N_eq=(N_eq-2)
          }
        }
        else if(count<100){
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
      Pop_DA[k]<-N_eq
      P1_eq[[k]]<-P1
      P2_eq[[k]]<-P2_out
      for (l in 1:N_reso){
        Etot1m=Etot_conc*f_var[l]#Molar
        Etot2m=Etot_conc*(1-f_var[l])#Molar
        Vm1m=kcat*Etot1m
        Vm2m=kcat*Etot2m
        kf_act_m=kf*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_m=kinh*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_spe_m=kinh_spe*10^(-(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        KI_m=kinh_act_m/kf_act_m
        KI_spe_m=kinh_act_spe_m/kf_act_m##note: previously used to test the effect of a different reversibility level in the pathway
        KM1_f_m=(kr+kcat)/kf_act_m
        KM1_b_m=(kr+kcat)/kinh_act_m
        KM1_b_spe_m=(kr+kcat)/kinh_act_spe_m
        metaboChainm<-function(y){
          e1<-c()
          e2<-c()
          t2 <- VTm*(Sout-y[1])/(KT+(Sout+y[1])+Sout*y[1]/KT)-(eta*y[1]+Vm1m*y[1]/(KM1_f_m+y[1]+KI_m*y[2])-Vm1m*(y[2])/(KM1_b_m+y[2]+y[1]/KI_m))
          for (i in 1:(Path1_size-1)){
            e1[i] <- (Vm1m*y[i]/(KM1_f_m+y[i]+KI_m*y[i+1])-Vm1m*(y[i+1])/(KM1_b_m+y[i+1]+y[i]/KI_m))-(eta*y[i+1]+Vm1m*y[i+1]/(KM1_f_m+y[i+1]+KI_m*y[i+2])-Vm1m*(y[i+2])/(KM1_b_m+y[i+2]+y[i+1]/KI_m))
          }
          e1[Path1_size] <- (Vm1m*y[Path1_size]/(KM1_f_m+y[Path1_size]+KI_m*y[Path1_size+1])-Vm1m*(y[Path1_size+1])/(KM1_b_m+y[Path1_size+1]+y[Path1_size]/KI_m))-(eta*y[Path1_size+1]+P_p2*S_c/V_c*(y[Path1_size+1]-P2_out)+Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-(Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m)))
          e2[1] <- (Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m))-(eta*y[Path1_size+2]+Vm2m*y[Path1_size+2]/(KM1_f_m+y[Path1_size+2]+KI_m*y[Path1_size+3])-Vm2m*(y[Path1_size+3])/(KM1_b_m+y[Path1_size+3]+y[Path1_size+2]/KI_m))
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2m*y[Path1_size+1+j]/(KM1_f_m+y[Path1_size+1+j]+KI_m*y[Path1_size+2+j])-Vm2m*(y[Path1_size+2+j])/(KM1_b_m+y[Path1_size+2+j]+y[Path1_size+1+j]/KI_m))-(eta*y[Path1_size+2+j]+Vm2m*y[Path1_size+2+j]/(KM1_f_m+y[Path1_size+2+j]+KI_m*y[Path1_size+3+j])-Vm2m*(y[Path1_size+3+j])/(KM1_b_m+y[Path1_size+3+j]+y[Path1_size+2+j]/KI_m))
          }
          e2[Path_size-Path1_size]<-(Vm2m*y[Path_size]/(KM1_f_m+y[Path_size]+KI_m*y[Path_size+1])-Vm2m*(y[Path_size+1])/(KM1_b_m+y[Path_size+1]+y[Path_size]/KI_m))-(eta*y[Path_size+1]+Vm2m*y[Path_size+1]/(KM1_f_m+y[Path_size+1]))
          c(t2,e1,e2)
        }
        Reso_m<-nleqslv(initm,metaboChainm,jac=NULL)$x
        Sin_m<-Reso_m[1]
        P1_m<-Reso_m[2:(Path_size+1)]
        Phi1m=0
        Phi2m=0
        Phi1m=Phi1m+(Vm1m*Sin_m/(KM1_f_m+Sin_m+KI_m*P1_m[1])-Vm1m*(P1_m[1])/(KM1_b_m+P1_m[1]+Sin_m/KI_m))/((Path_size/2))
        for (p1 in 1:(Path1_size-1)){
          Phi1m=Phi1m+(Vm1m*P1_m[p1]/(KM1_f_m+P1_m[p1]+KI_m*P1_m[p1+1])-Vm1m*(P1_m[p1+1])/(KM1_b_m+P1_m[p1+1]+P1_m[p1]/KI_m))/((Path_size/2))
        }
        Phi2m=Phi2m+(Vm2m*P1_m[Path1_size]/(KM1_f_m+P1_m[Path1_size]+KI_spe_m*P1_m[Path1_size+1])-Vm2m*(P1_m[Path1_size+1])/(KM1_b_spe_m+P1_m[Path1_size+1]+P1_m[Path1_size]/KI_spe_m))/((Path_size/2))
        for (p2 in (Path1_size+1):(Path_size-1)){
          Phi2m=Phi2m+(Vm2m*P1_m[p2]/(KM1_f_m+P1_m[p2]+KI_m*P1_m[p2+1])-Vm2m*(P1_m[p2+1])/(KM1_b_m+P1_m[p2+1]+P1_m[p2]/KI_m))/((Path_size/2))
        }
        Phi_fit_m<-((Phi1m*PO1+Phi2m*PO2)-prot_cost*(Etot1m*Path1_size+Etot2m*(Path_size-Path1_size)+Etot_back))*Tox/(Tox+sum(P1_m)+Sin_m)
        if(Phi_fit_m>0){}else{
          Phi_fit_m=0
        }
        if(N_eq>1){
          fit=Phi_fit_m/Phi_fit-1
        }
        else{
          fit=-100
        }
        tab_Phi_inv[k,l]=Phi_fit_m
        fit_inv[k,l]=fit
      }
    }
    list_fit_inv_S40_Yvar_Rev[[re]][[pe]]<-fit_inv
  }
}


multiplePlot("Main=","","Deg. rate=","",c("1:1","1:2"),c("1e-4","1e-3","&")
             ,ncol=800,f_var_print,f_var_print,list_fit_inv_S40_Yvar_Rev[[3]],
             abs="res. (% 1st pathway)",ord="mut. (% 1st pathway)",scale=c(-99,99),lev=c(0),palette=pal,cextext=1.5,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

Rep_opt_DA_RevTox_S40_perm<-list()
for (re in 1:length(Reversibility)){
  Rev_i<-as.character(Reversibility[re])
  Rep_opt_DA_RevTox_S40_perm[[re]]<-c(0)
  for (pe in 1:length(P_p2_set)){
    i=1
    while(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      #print(i)
      if(i>1 && i<(N_reso-1)){
        if(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i+1,i]>0 && list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i-1,i]>0){
          print(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i])
          S_opt=i
          print(S_opt)
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]<-f_var[S_opt]
          break
        }
      }
      else if (i==(N_reso-1)){
        if(list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i+1,i]<0 && list_fit_inv_S40_Yvar_Rev[[re]][[pe]][i,i+1]>0){
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]=1
          i=i+1
          break
        }
        else{
          Rep_opt_DA_RevTox_S40_perm[[re]][pe]<-0
          break
        }
      }
      i=i+1
    }
  }
}

##Main plot about the adaptive dynamics outcomes with permeability
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1,0.5,1))
plot(Rep_opt_DA_RevTox_S40_perm[[1]]~log10(P_p2_set),cex=0.02,col=1,pch=15,ylim=c(-0,1),xlab="log10(Permeability rate)",ylab=expression(paste(delta," (% 1st pathway)")))
points(Rep_opt_DA_RevTox_S40_perm[[1]][1:3]~log10(P_p2_set)[1:3],col="light green",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[1]][1:3]~log10(P_p2_set)[1:3],col="dark green",pch=16,cex=1)
points(Rep_opt_DA_RevTox_S40_perm[[1]][4]~log10(P_p2_set)[4],col="light green",pch=1,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[1]][4]~log10(P_p2_set)[4],col="dark green",pch=16,cex=1)
points(Rep_opt_DA_RevTox_S40_perm[[1]][5:6]~log10(P_p2_set)[5:6],col="light green",pch=1,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[1]][5:6]~log10(P_p2_set)[5:6],col="dark green",pch=1,cex=1)
points(Rep_opt_DA_RevTox_S40_perm[[2]]~log10(P_p2_set),col="light blue",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[2]]~log10(P_p2_set),col="blue",pch=16)
points(Rep_opt_DA_RevTox_S40_perm[[3]]~log10(P_p2_set),col="lightsalmon",pch=19,cex=2.5)
points(Rep_opt_DA_RevTox_S40_perm[[3]]~log10(P_p2_set),col="brown",pch=16)
text(-9,1,"D",srt=addtxt$srt,font=2,col="black",cex=1)
#text(-3.8,0.18,"*",srt=addtxt$srt,font=2,col="black",cex=1)
#text(-3.8,0.1,"o",srt=addtxt$srt,font=2,col=" light blue",cex=1)
par(xpd=TRUE)
legend("bottomleft",title="Singular strat.",legend=c("Loc. stable","Loc. instable","Glob. stable","Glob. instable"),bty="y",pch=c(20,1,19,1),col=c(1,"dark grey",1,1),pt.cex=c(1,1,1.75,1.75),ncol=1,cex=0.9)
legend("topright",title="Reversibility spread",legend=c("mostly kinh","equal","mostly kr"),bty="y",pch=c(15),col=c("green","blue","red"),ncol=1,cex=0.9)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "Cross-feeding-PO1_10-toxhighrev.jpeg", width = 575*6,height=500*6,res=600,type="cairo")
dev.off()
