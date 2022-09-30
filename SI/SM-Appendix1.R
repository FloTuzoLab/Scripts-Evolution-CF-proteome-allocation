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

options(digits=10)#To avoid problems with rounding numbers
#Defining features of plots
ncol=256
ncontour=5
#palet<-colorRampPalette(c("white","black"))(ncol)
jet.colors <- colorRampPalette(c("steelblue1", "yellow", "tomato"))
palet<-jet.colors(ncol)
pal<-list(palet)

#1.Trait Evolution Plot TEP for PO=1:10, eta=1e-1.5 and P>=10^-5 (case studied in the main paper)
##a.Determining the total concentration

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
kf=10^6.5##set to a moderately high value
kcat=10^2.5#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value 
PayOff_set<-c(0.1)
P_p2_set<-c(0)
Tox_set=c(10^10)
Path_size_set=c(40)

N_reso=100#Resolution used for strategies (in terms of concentration)
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
eta_set=c(10^-1.5)

#Variables
Etot_var<-10^seq(-5,-4,length=N_reso)

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
P_p2=P_p2_set[1]
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

Conc1_opt_DA_E0E1comb<-list()
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

#b.Adaptive dynamics with 2 sub-pathways to determine areas of coexistence

##Constant
ratio_set<-c(0.5) ##Size of each subpathway: 0.5 yields same size for both of them
#PayOff_set<-c(0.1)
#VTm<-10^-3

#Parameter of interest
VTm_max<-2*10^-3
Carr_cost<-10^-1.5
P_p2_set<-c(1e-5,1e-4)

##Redefining space of adaptive parameters to coincide with an allocation
N_reso=25
f_var<-seq(0,1,length=N_reso)

##New outcome list with new resolution
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
fit_inv_rec<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)
list_fit_inv_coex_5<-list()
list_fit_inv_coex_5_rec<-list()

for (p in 1:length(eta_set)){
  list_fit_inv_coex_5[[p]]<-list()
  list_fit_inv_coex_5_rec[[p]]<-list()
  eta<-eta_set[p]
  PayOff<-PayOff_set[1]
  ratio_12<-ratio_set[1]
  Path1_size=P1s_func(Pts=Path_size,ratio1_2=ratio_12)
  for (pe in 1:length(P_p2_set)){
    Etot_conc=Conc1_opt_DA_E0E1comb[[as.character(PayOff)]][p]*2# Double so that investing half in each sub-pathway coincides with the optimum when considered as a single integrated pathway
    P_p2=P_p2_set[pe]
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
      VTm=VTm_max*f_var[i]
      Etot0=Etot_conc*f_var[i] #First concentration set to avoid excess of upstreamn investment due to transport
      print(paste(p,pe,i))
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
        for (p2 in 1:(Path_size-Path1_size)){
          Phi2r=Phi2r+(Vm2r*P2[p2]/(P2[p2]+KM))/(Path_size/2)
        }
        Phi_fit<-((Phi1r*PO1+Phi2r*PO2)-VTm*Carr_cost-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1)+sum(P2))#-VTm*10^-1.5
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
      for (p2 in 1:(Path_size-Path1_size)){
        Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
      }
      Phi_res<-((Phi1r*PO1+Phi2r*PO2)-VTm*Carr_cost-prot_cost*(Etot_conc_r+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))#-VTm*10^-1.5
      for (j in 1:N_reso){
        VTm=VTm_max*f_var[j]
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
        for (p2 in 1:(Path_size-Path1_size)){
          Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
        }
        
        if(as.numeric((Phi1*PO1+Phi2*PO2)-prot_cost*((E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size))+(Etot_back)))>0){
          P_Phi_eq=as.numeric(((Phi1*PO1+Phi2*PO2)-VTm*Carr_cost-prot_cost*(Etot_conc_m+Etot_back))*Tox/(Tox+sum(P1_conc_eq)+sum(P2_conc_eq)))#
        }
        else{
          P_Phi_eq=0
        }
        if(N_eq>1){
          fit=P_Phi_eq/Phi_res-1
        }
        else{
          fit=-1000
        }
        tab_Phi_inv[i,j]=P_Phi_eq
        fit_inv[i,j]=fit
        fit_inv_rec[j,i]=fit
      }
    }
    list_fit_inv_coex_5[[p]][[pe]]<-fit_inv
    list_fit_inv_coex_5_rec[[p]][[pe]]<-fit_inv_rec
  }
}

jet.colors <- colorRampPalette(c(rep("red",1),rep("green",1)))
palet<-jet.colors(2)
pal<-list(palet)
multiplePlot("","","","",c(""),c("P=1e-5","P=1e-4")
             ,ncol=2,f_var_print,f_var_print,list_fit_inv_coex_5[[1]],
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-300,300),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="FALSE")

list_fit_inv_coex_5_area<-list()
fit_inv_coex<-matrix(nrow=N_reso , ncol=N_reso)

for (p in 1:length(eta_set)){
  list_fit_inv_coex_5_area[[p]]<-list()
  for (pe in 1:length(P_p2_set)){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        print(paste(p,pe,i))
        if(abs(list_fit_inv_coex_5[[p]][[pe]][i,j])<10^-5 && abs(list_fit_inv_coex_5_rec[[p]][[pe]][i,j])<10^-5){
          fit_inv_coex[i,j]=0
        }
        else if(list_fit_inv_coex_5[[p]][[pe]][i,j]==-1000 && list_fit_inv_coex_5_rec[[p]][[pe]][i,j]==-1000){
          fit_inv_coex[i,j]=-100
        }
        else if(list_fit_inv_coex_5_rec[[p]][[pe]][i,j]==-1000 && list_fit_inv_coex_5[[p]][[pe]][i,j]>0){
          fit_inv_coex[i,j]=1
        }
        else if(list_fit_inv_coex_5[[p]][[pe]][i,j]==-1000 && list_fit_inv_coex_5_rec[[p]][[pe]][i,j]>0){
          fit_inv_coex[i,j]=1
        }
        else if(list_fit_inv_coex_5[[p]][[pe]][i,j]>0 && list_fit_inv_coex_5_rec[[p]][[pe]][i,j]>0){
          fit_inv_coex[i,j]=1
        }
        else{
          fit_inv_coex[i,j]=-1
        }
      }
    }
    list_fit_inv_coex_5_area[[p]][[pe]]<-fit_inv_coex
  }
}

jet.colors <- colorRampPalette(c("tomato","steelblue1"))
palet<-jet.colors(2)
pal<-list(palet)
addtxt<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
multiplePlot("","","","",c(""),c("P=1e-5","P=1e-4")
             ,ncol=2,f_var_print,f_var_print,list_fit_inv_coex_5_area[[1]],
             abs=expression(paste(delta," resident")),ord=expression(paste(delta," mutant")),scale=c(-99,99),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="FALSE")


##c.Demographic equilibrium and evolutionary trajectories or coexistence areas
f_var<-seq(0,1,length=N_reso)
Neq1_tab<-matrix(nrow=N_reso,ncol=N_reso)
Neq2_tab<-matrix(nrow=N_reso,ncol=N_reso)
Neq1_list<-list()
Neq2_list<-list()
list_fit_gradinv_P5<-list()
fit_gradinv_P5<-list()
for (i in 1:N_reso){
  fit_gradinv_P5[[i]]<-list()
  for (j in 1:N_reso){
    fit_gradinv_P5[[i]][[j]]<-c(0)
  }
}

rm_coex<-c(1,-1)
for (p in 1:length(eta_set)){
  Etot_conc=Conc1_opt_DA_E0E1comb[[as.character(PayOff)]][p]*2
  Path1_size=P1s_func(Pts=Path_size,ratio1_2=ratio_12)
  for (pe in 1:length(P_p2_set)){
    P_p2=P_p2_set[pe]
    eta_p=eta+P_p2*S_c/V_c
    PO1<-1
    PO2<-1
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    r1=0
    while(r1<(N_reso)){
      r1=r1+1
      r2=0
      while(r2<(N_reso)){
        r2=r2+1
        print(paste(r1,r2))
        if(list_fit_inv_coex_5_area[[p]][[pe]][r1,r2]==1){
          #Etot_conc=2*Conc_opt_DA_Degperm_S40_noperm_eta1[1]
          if(r1==1 || r2==1){
            VTm<-VTm_max*c(f_var[r1],f_var[r2])+VTm_max/10^4
          }
          else{
            VTm<-VTm_max*c(f_var[r1],f_var[r2])
          }
          Etot0r1=Etot_conc*f_var[r1]
          Etot1r1=Etot_conc*f_var[r1]#Molar
          Etot2r1=Etot_conc*(1-f_var[r1])#Molar
          Etot0r2=Etot_conc*f_var[r2]
          Etot1r2=Etot_conc*f_var[r2]#Molar
          Etot2r2=Etot_conc*(1-f_var[r2])#Molar
          #Intermediate constants
          Vm0r=kcat*c(Etot0r1,Etot0r2)
          Vm1r=kcat*c(Etot1r1,Etot1r2)
          Vm2r=kcat*c(Etot2r1,Etot2r2)
          kfactr=kf*10^(-(c(Etot1r1,Etot1r2)*Path1_size+c(Etot2r1,Etot2r2)*(Path_size-Path1_size)+Etot_back)/(E_basal))
          KMr=(kr+kcat)/kfactr
          count=0
          N_eq1_p=0
          N_eq2_p=0
          N_eq1=2
          N_eq2=2
          Phi_fit<-c(0,0)
          P2_out=1
          while((abs(Phi_fit[1]-Phi_eq)/Phi_eq>10^-6) && (N_eq1>1 || N_eq2>1) && (abs(Phi_fit[2]-Phi_eq)/Phi_eq>10^-6)){
            transport<-function(x){
              t1 <- alpha-beta*x[3]-V_c/V_env*(N_eq1*VTm[1]*(x[3]-x[1])/(KT+(x[1]+x[3])+x[1]*x[3]/KT)+N_eq2*VTm[2]*(x[3]-x[2])/(KT+(x[2]+x[3])+x[2]*x[3]/KT))
              type=1
              t2<-c()
              while(type<3){
                t2[type] <- VTm[type]*(x[3]-x[type])/(KT+(x[type]+x[3])+x[3]*x[type]/KT)-Vm0r[type]*x[type]/(x[type]+KMr[type])
                type=type+1
              }
              c(t1,t2)
            }
            S<-nleqslv(c(1e-4,1e-6,1e-6),transport,jac=NULL)$x
            pathway1<-function(y,Vm1c,Vm2c,KM){
              a_1=eta*(KM+y)
              b_1=Vm2c*(y+KM)+eta*KM*(y+KM)-Vm1c*y
              c_1=-Vm1c*KM*y
              delta_1<-b_1^2-4*a_1*c_1
              P=(-b_1+delta_1^(1/2))/(2*a_1)
              return(P)
            }
            Sin<-c(S[1],S[2])
            P1_1<-c()
            P1_2<-c()
            P1_1[1]<-pathway1(Sin[1],Vm0r[1],Vm1r[1],KMr[1])
            P1_2[1]<-pathway1(Sin[2],Vm0r[2],Vm1r[2],KMr[2])
            for (p1 in 2:(Path1_size-1)){
              P1_1[p1]<-pathway1(P1_1[p1-1],Vm1r[1],Vm1r[1],KMr[1])
              P1_2[p1]<-pathway1(P1_2[p1-1],Vm1r[2],Vm1r[2],KMr[2])
            }
            P2_1<-c()
            P2_2<-c()
            pathway2<-function(z){
              type=1
              p1<-c()
              p1[type] <- Vm1r[type]*P1_1[Path1_size-1]/(P1_1[Path1_size-1]+KMr[type])-(Vm2r[type]*z[type]/(z[type]+KMr[type])+eta*z[type]+P_p2*S_c/V_c*(z[type]-z[3]))
              type=type+1
              p1[type] <- Vm1r[type]*P1_2[Path1_size-1]/(P1_2[Path1_size-1]+KMr[type])-(Vm2r[type]*z[type]/(z[type]+KMr[type])+eta*z[type]+P_p2*S_c/V_c*(z[type]-z[3]))
              p2 <- N_eq1*P_p2*S_c/V_env*(z[1]-z[3])+N_eq2*P_p2*S_c/V_env*(z[2]-z[3])-beta*z[3]
              c(p1,p2)
            }
            P2_sol<-nleqslv(c(1e-4,1e-6,1e-6),pathway2,jac=NULL)$x
            P2_1<-c()
            P2_2<-c()
            P2_1[1]<-P2_sol[1]
            P2_2[1]<-P2_sol[2]
            for (p2 in 2:(Path_size-Path1_size)){
              P2_1[p2]<-pathway1(P2_1[p2-1],Vm2r[1],Vm2r[1],KMr[1])
              P2_2[p2]<-pathway1(P2_2[p2-1],Vm2r[2],Vm2r[2],KMr[2])
            }
            Phi1r<-c()
            Phi2r-c()
            type=1
            while(type<3){
              Phi1r[type]=0
              Phi2r[type]=0
              Phi1r[type]=Phi1r[type]+(Vm0r[type]*Sin[type]/(Sin[type]+KMr[type]))/(Path_size/2)
              type=type+1
            }
            for (p1 in 1:(Path1_size-1)){
              Phi1r[1]=Phi1r[1]+(Vm1r[1]*P1_1[p1]/(P1_1[p1]+KMr[1]))/(Path_size/2)
              Phi1r[2]=Phi1r[2]+(Vm1r[2]*P1_2[p1]/(P1_2[p1]+KMr[2]))/(Path_size/2)
            }
            for (p2 in 1:(Path_size-Path1_size)){
              Phi2r[1]=Phi2r[1]+(Vm2r[1]*P2_1[p2]/(P2_1[p2]+KMr[1]))/(Path_size/2)
              Phi2r[2]=Phi2r[2]+(Vm2r[2]*P2_2[p2]/(P2_1[p2]+KMr[2]))/(Path_size/2)
            }
            
            Phi_fit[1]<-(Phi1r[1]*PO1+Phi2r[1]*PO2)-VTm*Carr_cost-prot_cost*(Etot1r1*Path1_size+Etot2r1*(Path_size-Path1_size)+Etot_back)
            Phi_fit[2]<-(Phi1r[2]*PO1+Phi2r[2]*PO2)-VTm*Carr_cost-prot_cost*(Etot1r2*Path1_size+Etot2r2*(Path_size-Path1_size)+Etot_back)
            if(count<10){
              if(Phi_fit[1]>Phi_eq){
                N_eq1=N_eq1*2
              }
              else{
                N_eq1=N_eq1/2
              }
              if(Phi_fit[2]>Phi_eq){
                N_eq2=N_eq2*2
              }
              else{
                N_eq2=N_eq2/2
              }
            }
            else if(count<30){
              N_eq1=N_eq1*Phi_fit[1]/Phi_eq
              N_eq2=N_eq2*Phi_fit[2]/Phi_eq
            }
            else{
              if(count==200){
                N_eq1=N_eq1*Phi_fit[1]/Phi_eq
                N_eq2=N_eq2*Phi_fit[2]/Phi_eq
              }
              print(paste(N_eq1,N_eq2))
              if(N_eq1<0){N_eq1=0}
              else if(Phi_fit[1]>Phi_fit[2]){
                if(Phi_fit[1]>Phi_eq){N_eq1=(N_eq1+1)}
                else{
                  N_eq1=(N_eq1-1)
                }
              }
              if(N_eq2<0){N_eq2=0}
              else if(Phi_fit[2]>Phi_fit[1]){
                if(Phi_fit[2]>Phi_eq){N_eq2=(N_eq2+1)}
                else{
                  N_eq2=(N_eq2-1)
                }
              }
              if(N_eq2_p==N_eq2 && N_eq1_p==N_eq1){
                break
              }
              if(count%%10==0){
                N_eq2_p=N_eq2
                N_eq1_p=N_eq1
              }
            }
            count=count+1
          }
          if(N_eq1<1){
            Neq1_tab[r1,r2]<-0.99
          }
          else{
            Neq1_tab[r1,r2]<-N_eq1
          }
          if(N_eq2<1){
            Neq2_tab[r1,r2]<-0.99
          }
          else{
            Neq2_tab[r1,r2]<-N_eq2
          }
          print(paste(r1,r2,N_eq1,N_eq2))
          Se<-S[3]
          P2_out<-P2_sol[3]
          ar<- (VTm*kfactr+kfactr*c(Etot0r1,Etot0r2)*kcat*(1+Se/KT))
          br<- (VTm*(kr+kcat-Se*kfactr)+kfactr*kcat*c(Etot0r1,Etot0r2)*(KT+Se))
          cr<- -VTm*Se*(kr+kcat)
          deltar<-br^2-4*ar*cr
          #First cellular concentration (substrate)
          S_conc_eq<-c()
          for (r in 1:2){
            if(ar[r]==0){
              S_conc_eq[r]=0
            }else{
              S_conc_eq[r]=(-br[r]+deltar[r]^(1/2))/(2*ar[r])
            }
          }
          
          a_1r=eta*(KMr+S_conc_eq)
          b_1r=Vm1r*(S_conc_eq+KMr)+eta*KMr*(S_conc_eq+KMr)-Vm0r*S_conc_eq
          c_1r=-Vm0r*KMr*S_conc_eq
          delta_1r<-b_1r^2-4*a_1r*c_1r
          P1_conc_eq<-list()
          P1_conc_eq[[1]]<-c()
          P1_conc_eq[[1]]=(-b_1r+delta_1r^(1/2))/(2*a_1r)
          
          for (p1 in 2:(Path1_size-1)){
            P1_conc_eq[[p1]]<-c()
            a_11r=eta*(KMr+P1_conc_eq[[p1-1]])
            b_11r=Vm1r*(P1_conc_eq[[p1-1]]+KMr)+eta*KMr*(P1_conc_eq[[p1-1]]+KMr)-Vm1r*P1_conc_eq[[p1-1]]
            c_11r=-Vm1r*KMr*P1_conc_eq[[p1-1]]
            delta_11r<-b_11r^2-4*a_11r*c_11r
            P1_conc_eq[[p1]]=(-b_11r+delta_11r^(1/2))/(2*a_11r)
          }
          
          a_2r=eta_p*(KMr+P1_conc_eq[[Path1_size-1]])
          b_2r=Vm2r*(P1_conc_eq[[Path1_size-1]]+KMr)+(eta_p*KMr-P_p2*S_c/V_c*P2_out)*(P1_conc_eq[[Path1_size-1]]+KMr)-Vm1r*P1_conc_eq[[Path1_size-1]]
          c_2r=-Vm1r*KMr*P1_conc_eq[[Path1_size-1]]-P_p2*S_c/V_c*P2_out*(P1_conc_eq[[Path1_size-1]]+KMr)*KMr
          delta_2r<-b_2r^2-4*a_2r*c_2r
          P2_conc_eq<-list()
          P2_conc_eq[[1]]<-c()
          P2_conc_eq[[1]]=(-b_2r+delta_2r^(1/2))/(2*a_2r)
          
          for (p2 in 2:(Path_size-Path1_size)){
            P2_conc_eq[[p2]]<-c()
            a_22r=eta*(KMr+P2_conc_eq[[p2-1]])
            b_22r=Vm2r*(P2_conc_eq[[p2-1]]+KMr)+eta*KMr*(P2_conc_eq[[p2-1]]+KMr)-Vm2r*P2_conc_eq[[p2-1]]
            c_22r=-Vm2r*KMr*P2_conc_eq[[p2-1]]
            delta_22r<-b_22r^2-4*a_22r*c_22r
            P2_conc_eq[[p2]]=(-b_22r+delta_22r^(1/2))/(2*a_22r)
          }
          Phi1rr=0
          Phi2rr=0
          Phi1rr=Phi1rr+(Vm0r*S_conc_eq/(S_conc_eq+KMr))/(Path_size/2)
          for (p1 in 1:(Path1_size-1)){
            Phi1rr=Phi1rr+(Vm1r*P1_conc_eq[[p1]]/(P1_conc_eq[[p1]]+KMr))/(Path_size/2)
          }
          for (p2 in 1:(Path_size-Path1_size)){
            Phi2rr=Phi2rr+(Vm2r*P2_conc_eq[[p2]]/(P2_conc_eq[[p2]]+KMr))/(Path_size/2)
          }
          Phirr<-(Phi1rr*PO1+Phi2rr*PO2)-VTm*Carr_cost-prot_cost*(c(Etot1r1,Etot1r2)*Path1_size+c(Etot2r1,Etot2r2)*(Path_size-Path1_size)+(Etot_back))
          
          max_rm<-4
          for (rm in 1:max_rm){
            if(rm<3){
              if(r1==1){
                co_m<-r1+1
              }
              else if(r1==N_reso){
                co_m=r1-1
              }
              else{
                co_m<-r1+rm_coex[rm%%2+1]
              }
            }else{
              if(r2==1){
                co_m<-r2+1
              }
              else if(r2==N_reso){
                co_m=r2-1
              }
              else{
                co_m<-r2+rm_coex[rm%%2+1]
              }
            }
            print(paste(r1,r2,co_m))
            if(r1==1 || r2==1){
              VTm=VTm_max*f_var[co_m]+VTm_max/10^4
            }
            else{
              VTm=VTm_max*f_var[co_m]
            }
            E_tot_conc0=Etot_conc*f_var[co_m]
            E_tot_conc1<-Etot_conc*f_var[co_m]
            E_tot_conc2<-Etot_conc*(1-f_var[co_m])
            kf_act<-kf*10^(-(E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size)+(Etot_back))/(E_basal))#+Etot_back
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
            for (p2 in 1:(Path_size-Path1_size)){
              Phi2=Phi2+(Vm2*P2_conc_eq[p2]/(P2_conc_eq[p2]+KM))/(Path_size/2)
            }
            
            if(as.numeric((Phi1*PO1+Phi2*PO2)-VTm*Carr_cost-prot_cost*(E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size)+(Etot_back)))>0){
              P_Phi_eq=as.numeric((Phi1*PO1+Phi2*PO2)-VTm*Carr_cost-prot_cost*(E_tot_conc1*Path1_size+E_tot_conc2*(Path_size-Path1_size)+(Etot_back)))
            }
            else{
              P_Phi_eq=0
            }
            if(rm<3){
              fit=(P_Phi_eq-Phirr[1])/Phi_fit
            }else{
              fit=(P_Phi_eq-Phirr[2])/Phi_fit
            }
            
            if(fit>0){
              fit_gradinv_P5[[r1]][[r2]][rm]<-1
            }
            else{
              fit_gradinv_P5[[r1]][[r2]][rm]<--1
            }
            if(r1==1){
              fit_gradinv_P5[[r1]][[r2]][1]<--1
            }
            if(r1==N_reso){
              fit_gradinv_P5[[r1]][[r2]][2]<--1
            }
            if(r2==1){
              fit_gradinv_P5[[r1]][[r2]][3]<--1
            }
            if(r2==N_reso){
              fit_gradinv_P5[[r1]][[r2]][4]<--1
            }
          }
        }
      }
    }
    Neq1_list[[pe+length(P_p2_set)*(p-1)]]<- Neq1_tab
    Neq2_list[[pe+length(P_p2_set)*(p-1)]]<- Neq2_tab
    list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]]<-fit_gradinv_P5
  }
}


##d.Plotting coalition dynamics (TEP)

trueVec<-function(vector1,vector2){
  for (v in 1:length(vector1)){
    if (vector1[v]!=vector2[v]){
      return(FALSE)
      break
    }
  }
  return(TRUE)
}
list_fit_coAlinv_P5<-list()
fit_coAlinv_P5<-matrix(nrow=N_reso,ncol=N_reso)

for (p in 1:length(eta_set)){
  for (pe in 1:length(P_p2_set)){
    for (i in  1:N_reso){
      for (j in 1:N_reso){
        #print()
        print(paste(i,j))
        if(length(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]])==1){
          fit_coAlinv_P5[i,j]=-1
        }
        else if(list_fit_inv_coex_5_area[[p]][[pe]][i,j]==-1){
          fit_coAlinv_P5[i,j]=-1
        }
        else if(length(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]])==4){
          if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,1,1,1))){
            fit_coAlinv_P5[i,j]=-1#repell
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,-1,-1,-1))){
            fit_coAlinv_P5[i,j]=0#css
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,-1,-1,-1))){
            fit_coAlinv_P5[i,j]=1#l
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,1,-1,-1))){
            fit_coAlinv_P5[i,j]=2#r
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,-1,1,-1))){
            fit_coAlinv_P5[i,j]=3#d
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,-1,-1,1))){
            fit_coAlinv_P5[i,j]=4#t
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,1,-1,-1))){
            fit_coAlinv_P5[i,j]=5#lr
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,-1,1,-1))){
            fit_coAlinv_P5[i,j]=6#ld
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,-1,-1,1))){
            fit_coAlinv_P5[i,j]=7#lt
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,1,1,-1))){
            fit_coAlinv_P5[i,j]=8#rd
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,1,-1,1))){
            fit_coAlinv_P5[i,j]=9#rt
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,-1,1,1))){
            fit_coAlinv_P5[i,j]=10#td
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,1,1,-1))){
            fit_coAlinv_P5[i,j]=11#lrd
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,-1,1,1))){
            fit_coAlinv_P5[i,j]=12#ldt
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(1,1,-1,1))){
            fit_coAlinv_P5[i,j]=13#lrt
          }
          else if(trueVec(list_fit_gradinv_P5[[pe+length(P_p2_set)*(p-1)]][[i]][[j]],c(-1,1,1,1))){
            fit_coAlinv_P5[i,j]=14#rtd
          }
        }
      }
    }
    list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]]<-fit_coAlinv_P5
  }
}


dev.off()
par(mar=c(2.5, 2.5, 2, 2), mfrow=c(1,1), las=1,xpd=TRUE)
addtxt<-list(l=0.05,h=0.95,txt=c("A","B","C"),srt = 0,font=2,col="black")
size=0.02
delta_arr=4
for (p in 1:length(eta_set)){
  for (pe in 1:1){
    multiplePlot("","","","",c(""),c("")
                 ,ncol=2,f_var_print,f_var_print,list(list_fit_inv_coex_5_area[[p]][[2]]),
                 abs=expression(paste(delta," resident 1 (% 1st pathway)")),ord=expression(paste(delta," resident 2 (% 1st pathway)")),scale=c(-80,80),palette=pal,cextext=2,TEXT_to_Add=addtxt,
                 image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=1,globcex=0.65,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="FALSE")
    i=1
    while(i<N_reso+1){
      j=1
      while(j<N_reso+1){
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==0){
          points(f_var[i], f_var[j],col="green",pch=16,cex=3)
        }
        j=j+1
      }
      i=i+1
    }
    i=1
    while(i<N_reso+1){
      j=1
      while(j<N_reso+1){
        print(paste(i,j))
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==-1){
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==1){
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==2){
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==3){
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==4){
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==5){#lr
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==6){#ld
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==7){#lt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==8){#rd
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==9){#rt
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==10){#td
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==11){#lrd
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==12){#ldt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==13){#lrt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==14){#rtd
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        j=j+delta_arr
      }
      i=i+delta_arr
    }
    par(xpd=FALSE)
    abline(a=0,b=1)
  }
}

par(xpd=TRUE)
points(x=0,y=1,col="black",pch=19,cex=0.5)
points(x=1,y=0,col="black",pch=19,cex=0.5)
points(x=0.66,y=0.66,,col="black",pch=13,cex=3.5)
legend("bottomleft",title="Singular Strategy",legend=c("Loc. stable","Glob. stable","Glob. instable"),bty="y",pch=c(20,19,13,19),col=c(1,"green",1,"blue"),ncol=1,cex=1.5)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "Coex-perm-5.jpeg", width = 550*6,height=500*6,res=600,type="cairo")
dev.off()   

par(mar=c(2.5, 2.5, 2, 2), mfrow=c(1,1), las=1,xpd=TRUE)
addtxt<-list(l=0.05,h=0.95,txt=c("B","C"),srt = 0,font=2,col="black")
size=0.02
delta_arr=4
for (p in 1:length(eta_set)){
  for (pe in 1:1){
    multiplePlot("","","","",c(""),c("")
                 ,ncol=2,f_var_print,f_var_print,list(list_fit_inv_coex_5_area[[p]][[2]]),
                 abs=expression(paste(delta," resident 1 (% 1st pathway)")),ord=expression(paste(delta," resident 2 (% 1st pathway)")),scale=c(-80,80),palette=pal,cextext=2,TEXT_to_Add=addtxt,
                 image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=1,globcex=0.65,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="FALSE")
    i=1
    while(i<N_reso+1){
      j=1
      while(j<N_reso+1){
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==0){
          points(f_var[i], f_var[j],col="green",pch=16,cex=3)
        }
        j=j+1
      }
      i=i+1
    }
    i=1
    while(i<N_reso+1){
      j=1
      while(j<N_reso+1){
        print(paste(i,j))
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==-1){
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==1){
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==2){
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==3){
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==4){
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==5){#lr
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==6){#ld
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==7){#lt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==8){#rd
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==9){#rt
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==10){#td
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==11){#lrd
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==12){#ldt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==13){#lrt
          arrows(f_var[i], f_var[j], f_var[i]-size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        if(list_fit_coAlinv_P5[[pe+length(P_p2_set)*(p-1)]][i,j]==14){#rtd
          arrows(f_var[i], f_var[j], f_var[i]+size, f_var[j], length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]-size, length=0.025, col="black", lwd=1)
          arrows(f_var[i], f_var[j], f_var[i], f_var[j]+size, length=0.025, col="black", lwd=1)
        }
        j=j+delta_arr
      }
      i=i+delta_arr
    }
    par(xpd=FALSE)
    abline(a=0,b=1)
  }
}

par(xpd=TRUE)
points(x=0,y=1,col="black",pch=19,cex=0.5)
points(x=1,y=0,col="black",pch=19,cex=0.5)
points(x=0.66,y=0.66,,col="black",pch=13,cex=3.5)
legend("bottomleft",title="Singular Strategy",legend=c("Loc. stable","Glob. stable","Glob. instable"),bty="y",pch=c(20,19,13,19),col=c(1,"green",1,"blue"),ncol=1,cex=1.5)
setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "Coex-perm-4.jpeg", width = 550*6,height=500*6,res=600,type="cairo")
dev.off()  