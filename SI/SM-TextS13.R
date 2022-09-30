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

#1.Adaptive dynamics showing concentration ceilings for different toxicity levels##SM Figure S4-A
##Costs sensitivity

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
prot_cost_set=10^-2.5#Mmet/Menz##set to an average value where 
PayOff_set<-c(1)
Path_size_set=c(40)

N_reso=150#Resolution used for strategies (in terms of concentration)
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

##Parameters

Etot0_set<-c(0)
Tox_set=c(10^-2,10^-1,1)#not used here
eta_set=10^seq(-4,-1,0.5)
P_p2_set<-c(0)#not used here

#Pathway parameter set
PayOff<-PayOff_set[1]
Path_size=Path_size_set[1]
Path1_size=Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)

##Variables

#Focal
Etot_var<-10^seq(-5.5,-4,length=N_reso)

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

for (etot in 1:length(Etot0_set)){
  list_fit_inv_S40_noperm[[etot]]<-list()
  Enz=1
  for(vtox in 1:length(Tox_set)){
    Tox<-Tox_set[vtox]
    for (p in 1:length(eta_set)){
      print(paste(etot,vtox,p))
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
        if(etot<length(Etot0_set)){
          Etot0=Etot0_set[etot]
        }
        else{
          Etot0=Etot_var[i]
        }
        Etot1r=Etot_var[i]#Molar
        Etot2r=Etot_var[i]#Molar
        Etot_conc=Etot_var[i]*(Path_size-1)+Etot0
        #Intermediate constants
        Vm0=kcat*Etot0
        Vm1r=kcat*Etot1r
        Vm2r=kcat*Etot2r
        kfact=kf*10^(-(Etot_conc+Etot_back)/(E_basal))
        KM=(kr+kcat)/kfact
        #KM=(kr+kcat)/kf
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
        for (p2 in 1:(Path_size-Path1_size)){
          Phi2r=Phi2r+(Vm2r*P2_res[p2]/(P2_res[p2]+KM))/(Path_size/2)
        }
        Phi_res<-((Phi1r*PO1+Phi2r*PO2)-prot_cost*(Etot_conc+Etot_back))*Tox/(Tox+sum(P1_res)+sum(P2_res))
        ##Determining invasion fitness
        for (j in 1:N_reso){
          if(etot<length(Etot0_set)){
            E_tot_conc0=Etot0_set[etot]
          }
          else{
            E_tot_conc0=Etot_var[j]
          }
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
          for (p2 in 1:(Path_size-Path1_size)){
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(vtox-1)]]<-max(Pop_DA)
    }
  }
}


addtxt<-list(l=0.075,h=1.00,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"),srt = 0,font=2,col="gray")
for (i in 1:N_reso){
  print(fit_inv[i,i])
}
ncol=800
jet.colors <- colorRampPalette(c(rep("white",400),rep("grey20",400)))
#jet.colors <- colorRampPalette(c(rep("red",49),rep("green",49)))
palet<-jet.colors(ncol)
pal<-list(palet)
f_var_print<-c(0,1,0.2)
log10Etot1_var<-c(-5.5,-4,0.5)
sublegend<-c(expression(paste(eta,"=",10^-4,"/s")),           
             expression(paste(eta,"=",10^-2,"/s")))

##Concentration of 1e-4.75 for E0 (enzyme after transporter)
#multiplePlot("Cost=","","Deg. rate=","",c("1e-3","1e-2.5","1e-2"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
#             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
#             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
#             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-475.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()

##Concentration of E0 (enzyme after transporter) similar to that downstream 
multiplePlot("Cost=","","Deg. rate=","",c("1e-3","1e-2.5","1e-2"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-var.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()


Conc1_opt_DA_eff_E045E1comb<-list()
Toxicity<-c("High","Moderate","Low")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
  for(vtox in 1:length(Tox_set)){
    Tox_i<-as.character(Toxicity[vtox])
    print(Tox_i)
    Conc1_opt_DA_eff_E045E1comb[[etot]][[Tox_i]]<-c()
    for (p in 1:length(eta_set)){
      i=2
      while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Tox_i]][p]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Tox_i]][p]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vtox-1)]][i,i+1]>0){
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Tox_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Tox_i]][p]<-0
            break
          }
        }
        else{
          i=i+1
        }
      }
    }
  }
}
Conc1_opt_DA_eff_E045E1comb
#Conc1_opt_DA_eff_E045E1comb[["Moderately low"]][5]<-0
Etot0_set_graph<-c("var")
txt<-c("A","B","C","D")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(vtox in 1:length(Tox_set)){
    Tox_i<-as.character(Toxicity[vtox])
    if (vtox==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Tox_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Tox_i)]])~log10(eta_set),col=vtox,pch=14+vtox)
    }
  }
  percent=c(0.05,0.1,0.15,0.2)
  for (p in 1:length(percent)){
    abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
  }
  legend("bottomright",title="Toxicity",legend=c("High","Moderate","Low"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.75)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  text(-3.25,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-4,-4.4,txt[etot],srt=addtxt$srt,font=2,col="black",cex=1)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-toxicity-firstenzconc.jpeg", , width = 500*3,height=415*3,res=300,type="cairo")
dev.off()


#2.Adaptive dynamics showing concentration ceilings for different reversibility levels##SM Figure S4-B

Tox=1
VTm=1e-3#M/s
KT=1e-2#M
#Se=KT*2
kf=10^6.5
kcat=10^2.5#/s
kr_set=c(kcat/9,kcat/3,kcat)#/s
kinh_set=c(kf,kf/3,kf/9)
kinh_spe_set=c(kf,kf/3,kf/9)
prot_cost=10^-2.5#Mmet/Menz
eta_set=10^seq(-4,-1,0.5)
PayOff_set<-c(1)
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
Etot_var<-10^seq(-5.5,-4,length=N_reso)
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
PayOff<-PayOff_set[1]
Path_size=Path_size_set[1]
Path1_size=Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)

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
  for (p in 1:length(eta_set)){
    print(paste(re,p))
    #P_p2=P_p2_set[p]
    eta<-eta_set[p]
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    for (k in 1:N_reso){
      #k=18
      #print(paste(k))
      Etot1r=Etot_var[k]#Molar
      Etot2r=Etot_var[k]#Molar
      #Intermediate constants
      #Vm0=kcat*Etot0
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
      #print(KI)
      #KM=(kr+kcat)/kf
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
    list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]]<-fit_inv
  }
}
for (i in 1:N_reso){
  print(fit_inv[i,i])
}


addtxt<-list(l=0.95,h=0.05,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
log10Etot1_var<-c(-6,-4,0.5)
multiplePlot("","","","",c("","",""),c("","","","","","")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_rev_noperm,
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


Conc_opt_DA_Degperm_S40_rev_noperm<-list()
Reversibility<-c("kr","both","kinh")
for(re in 1:length(kr_set)){
  Rev_i<-as.character(Reversibility[re])
  print(Rev_i)
  Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]]<-c()
  for (p in 1:length(eta_set)){
    #print(p)
    i=1
    while(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      if(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i+2]<0 && list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i-2]<0){
        print(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]<-Etot_var[S_opt]
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i+1,i]<0 && list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i+1]>0){
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]=1
          i=i+1
          break
        }
        else{
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]<-0
          break
        }
      }
      else{
        i=i+1
      }
    }
  }
}


#dev.off()
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(re in 1:length(kr_set)){
  Rev_i<-as.character(Reversibility[re])
  if (re==2){
    plot(log10(Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][1:7])~log10(eta_set[1:7]),col=2,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=","var"))
  }
  else{
    points(log10(Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][1:7])~log10(eta_set[1:7]),col=re,pch=14+re)
  }
}
percent=c(0.05,0.1,0.15,0.2)
for (p in 1:length(percent)){
  abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
}
legend("bottomright",title="Reversibility spread",legend=c("equal","mostly kr"),bty="y",pch=c(16:18),col=c(2:4),ncol=1,cex=1)
#text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
#text(-2.5,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-4,-4.4,"B",srt=addtxt$srt,font=2,col="black",cex=1)

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-toxrev-e0var.jpeg", , width = 500*3,height=415*3,res=300,type="cairo")
dev.off()

#3.Adaptive dynamics showing concentration ceilings for different reversibility levels - reversibility with first enzyme set

Tox=1
VTm=1e-3#M/s
KT=1e-2#M
#Se=KT*2
kf=10^6.5
kcat=10^2.5#/s
kr_set=c(kcat/9,kcat/3,kcat)#/s
kinh_set=c(kf,kf/3,kf/9)
kinh_spe_set=c(kf,kf/3,kf/9)
prot_cost=10^-2.5#Mmet/Menz
eta_set=10^seq(-3.5,-1.5,0.5)
PayOff_set<-c(1)
Path_size_set=c(40)
N_reso=50
Etot_back<-5.5*10^-3#other protein concentration
E_basal<-3*10^-3#scaling factor
#Parameters
#eta=10^-0.5

P_p2_set<-c(0)
#eta=eta_set[2]
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
Etot_var<-10^seq(-5.5,-4,length=N_reso)
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

#Outcome
tab_N_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
tab_S_eq <- c()
tab_Sin_eq<-matrix(nrow=N_reso , ncol=N_reso)
tab_P1<-matrix(nrow=N_reso , ncol=N_reso)
tab_Phi_inv<-matrix(nrow=N_reso , ncol=N_reso)

Neq<-list()
Phi_c<-list()
S_eq<-list()
Sin_eq<-list()
P1_eq<-list()
fit_inv<-matrix(nrow=N_reso , ncol=N_reso)
#list_fit_inv_S20rat14_Y1<-list()
Pop_DA<-c()
#Parameter set 1
PayOff<-PayOff_set[1]
Path_size=Path_size_set[1]
Path1_size=Path1_size=P1s_func(Pts=Path_size,ratio1_2=1/2)

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
  for (p in 1:length(eta_set)){
    print(paste(re,p))
    #P_p2=P_p2_set[p]
    eta<-eta_set[p]
    PO1<-5
    PO2<-5
    if(PayOff>1){
      PO1<-PayOff
    }
    else if (PayOff<1){
      PO2<-1/PayOff
    }
    for (k in 1:N_reso){
      #k=18
      #print(paste(k))
      Etot0=10^-4.5
      Etot1r=Etot_var[k]#Molar
      Etot2r=Etot_var[k]#Molar
      #Intermediate constants
      Vm0=kcat*Etot0
      Vm1r=kcat*Etot1r
      Vm2r=kcat*Etot2r
      kf_act=kf*10^(-(Etot0+Etot1r*(Path1_size-1)+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act=kinh*10^(-(Etot0+Etot1r*(Path1_size-1)+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      kinh_act_spe=kinh_spe*10^(-(Etot0+Etot1r*(Path1_size-1)+Etot2r*(Path_size-Path1_size)+Etot_back)/(E_basal))
      KI=kinh_act/kf_act
      KI_spe=kinh_act_spe/kf_act
      KM1_f=(kr+kcat)/kf_act
      KM1_b=(kr+kcat)/kinh_act
      KM1_b_spe=(kr+kcat)/kinh_act_spe
      #print(KI)
      #KM=(kr+kcat)/kf
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
          t2 <- VTm*(x[1]-x[2])/(KT+(x[1]+x[2])+x[1]*x[2]/KT)-(eta*x[2]+Vm0*x[2]/(KM1_f+x[2]+KI*x[3])-Vm0*(x[3])/(KM1_b+x[3]+x[2]/KI))
          e1[1] <- (Vm1r*x[2]/(KM1_f+x[2]+KI*x[3])-Vm1r*(x[3])/(KM1_b+x[3]+x[2]/KI))-(eta*x[3]+Vm0*x[3]/(KM1_f+x[3]+KI*x[4])-Vm0*(x[4])/(KM1_b+x[4]+x[3]/KI))
          for (i in 2:(Path1_size-1)){
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
        Phi1r=Phi1r+(Vm0*Sin/(KM1_f+Sin+KI*P1[1])-Vm0*(P1[1])/(KM1_b+P1[1]+Sin/KI))/((Path_size/2))
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
          ##print(paste(p,N_eq,N_eq_p))
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
        #print(paste(Phi_fit,N_eq))
      }
      Pop_DA[k]<-N_eq
      P1_eq[[k]]<-P1
      P2_eq[[k]]<-P2_out
      for (l in 1:N_reso){
        #l=25
        Etot1m=Etot_var[l]#Molar
        Etot2m=Etot_var[l]#Molar
        #Etot_conc=Etot0*Path_size
        #Intermediate constants
        #Vm0=kcat*Etot0
        Vm1m=kcat*Etot1m
        Vm2m=kcat*Etot2m
        kf_act_m=kf*10^(-(Etot0+Etot1m*(Path1_size-1)+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_m=kinh*10^(-(Etot0+Etot1m*(Path1_size-1)+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        kinh_act_spe_m=kinh_spe*10^(-(Etot0+Etot1m*(Path1_size-1)+Etot2m*(Path_size-Path1_size)+Etot_back)/(E_basal))
        KI_m=kinh_act_m/kf_act_m
        KI_spe_m=kinh_act_spe_m/kf_act_m##note: previously used to test the effect of a different reversibility level in the pathway
        KM1_f_m=(kr+kcat)/kf_act_m
        KM1_b_m=(kr+kcat)/kinh_act_m
        KM1_b_spe_m=(kr+kcat)/kinh_act_spe_m
        #KM=(kr+kcat)/kf
        metaboChainm<-function(y){
          e1<-c()
          e2<-c()
          #t1 <- alpha-beta*Sout-V_c/V_env*N_eq*VTm*(Sout-y[1])/(KT+(Sout+y[1])+Sout*y[1]/KT)
          t2 <- VTm*(Sout-y[1])/(KT+(Sout+y[1])+Sout*y[1]/KT)-(eta*y[1]+Vm0*y[1]/(KM1_f_m+y[1]+KI_m*y[2])-Vm0*(y[2])/(KM1_b_m+y[2]+y[1]/KI_m))
          e1[1] <- (Vm1m*y[1]/(KM1_f_m+y[1]+KI_m*y[2])-Vm1m*(y[2])/(KM1_b_m+y[2]+y[1]/KI_m))-(eta*y[2]+Vm0*y[2]/(KM1_f_m+y[2]+KI_m*y[3])-Vm0*(y[3])/(KM1_b_m+y[3]+y[2]/KI_m))
          for (i in 2:(Path1_size-1)){
            e1[i] <- (Vm1m*y[i]/(KM1_f_m+y[i]+KI_m*y[i+1])-Vm1m*(y[i+1])/(KM1_b_m+y[i+1]+y[i]/KI_m))-(eta*y[i+1]+Vm1m*y[i+1]/(KM1_f_m+y[i+1]+KI_m*y[i+2])-Vm1m*(y[i+2])/(KM1_b_m+y[i+2]+y[i+1]/KI_m))
          }
          e1[Path1_size] <- (Vm1m*y[Path1_size]/(KM1_f_m+y[Path1_size]+KI_m*y[Path1_size+1])-Vm1m*(y[Path1_size+1])/(KM1_b_m+y[Path1_size+1]+y[Path1_size]/KI_m))-(eta*y[Path1_size+1]+P_p2*S_c/V_c*(y[Path1_size+1]-P2_out)+Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-(Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m)))
          e2[1] <- (Vm2m*y[Path1_size+1]/(KM1_f_m+y[Path1_size+1]+KI_spe_m*y[Path1_size+2])-Vm2m*(y[Path1_size+2])/(KM1_b_spe_m+y[Path1_size+2]+y[Path1_size+1]/KI_spe_m))-(eta*y[Path1_size+2]+Vm2m*y[Path1_size+2]/(KM1_f_m+y[Path1_size+2]+KI_m*y[Path1_size+3])-Vm2m*(y[Path1_size+3])/(KM1_b_m+y[Path1_size+3]+y[Path1_size+2]/KI_m))
          #t3<-N_eq*P_p2*S_c/V_env*(x[Path1_size+2]-x[Path_size+3])-beta*x[Path_size+3]
          for (j in 1:(Path_size-(Path1_size+2))){
            e2[j+1] <- (Vm2m*y[Path1_size+1+j]/(KM1_f_m+y[Path1_size+1+j]+KI_m*y[Path1_size+2+j])-Vm2m*(y[Path1_size+2+j])/(KM1_b_m+y[Path1_size+2+j]+y[Path1_size+1+j]/KI_m))-(eta*y[Path1_size+2+j]+Vm2m*y[Path1_size+2+j]/(KM1_f_m+y[Path1_size+2+j]+KI_m*y[Path1_size+3+j])-Vm2m*(y[Path1_size+3+j])/(KM1_b_m+y[Path1_size+3+j]+y[Path1_size+2+j]/KI_m))
          }
          e2[Path_size-Path1_size]<-(Vm2m*y[Path_size]/(KM1_f_m+y[Path_size]+KI_m*y[Path_size+1])-Vm2m*(y[Path_size+1])/(KM1_b_m+y[Path_size+1]+y[Path_size]/KI_m))-(eta*y[Path_size+1]+Vm2m*y[Path_size+1]/(KM1_f_m+y[Path_size+1]))
          c(t2,e1,e2)
        }
        Reso_m<-nleqslv(initm,metaboChainm,jac=NULL)$x
        Sin_m<-Reso_m[1]
        P1_m<-Reso_m[2:(Path_size+1)]
        #print(paste(l,P1_m))
        #P2_out<-Reso[Path_size+3]
        Phi1m=0
        Phi2m=0
        Phi1m=Phi1m+(Vm0*Sin_m/(KM1_f_m+Sin_m+KI_m*P1_m[1])-Vm0*(P1_m[1])/(KM1_b_m+P1_m[1]+Sin_m/KI_m))/((Path_size/2))
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
    list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]]<-fit_inv
  }
}
for (i in 1:N_reso){
  print(fit_inv[i,i])
}


addtxt<-list(l=0.95,h=0.05,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
log10Etot1_var<-c(-6,-4,0.5)
multiplePlot("","","","",c("","",""),c("","","","","")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_rev_noperm,
             abs="resident (logEtot1)",ord="mutant (logEtot1)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,TEXT_to_Add=addtxt,
             image=TRUE,pcex=1,subcex=1,labcex=1.5,axcex=0.75,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


Conc_opt_DA_Degperm_S40_rev_noperm<-list()
Reversibility<-c("kr","both","kinh")
for(re in 1:length(kr_set)){
  Rev_i<-as.character(Reversibility[re])
  print(Rev_i)
  Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]]<-c()
  for (p in 1:length(eta_set)){
    #print(p)
    i=1
    while(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i]==-100 && i<N_reso){
      i=i+1
    }
    S_opt=0
    print(i)
    while (i<N_reso){
      if(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i+2]<0 && list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i-2]<0){
        print(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]<-Etot_var[S_opt]
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i+1,i]<0 && list_fit_inv_S40_rev_noperm[[p+length(eta_set)*(re-1)]][i,i+1]>0){
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]=1
          i=i+1
          break
        }
        else{
          Conc_opt_DA_Degperm_S40_rev_noperm[[Rev_i]][p]<-0
          break
        }
      }
      else{
        i=i+1
      }
    }
  }
}


dev.off()
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(re in 1:length(kr_set)){
  Rev_i<-as.character(Reversibility[re])
  if (re==1){
    plot(log10(Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][1:5])~log10(eta_set[1:5]),col=1,pch=15,ylim=c(-5.2,-4.0),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)",main=paste("[E0]=","var"))
  }
  else{
    points(log10(Conc_opt_DA_Degperm_S40_rev_noperm[[as.character(Rev_i)]][1:5])~log10(eta_set[1:5]),col=re,pch=14+re)
  }
}
percent=c(0.05,0.1,0.15,0.2)
for (p in 1:length(percent)){
  abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
}
legend("bottomright",title="Reversibility spread",legend=c("equal","mostly kr"),bty="y",pch=c(16:18),col=c(2:4),ncol=1,cex=1)
#text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
text(-2.5,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-1.55,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
text(-3.5,-4.4,"C",srt=addtxt$srt,font=2,col="black",cex=1)

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-toxrev-e0var.jpeg", , width = 560*3,height=475*3,res=300,type="cairo")
dev.off()
