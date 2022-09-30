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

#1.Adaptive dynamics showing concentration ceilings for different protein costs##SM Figure S3-A
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
prot_cost_set=c(10^-3,10^-2.5,10^-2)#Mmet/Menz##set to an average value where 
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
Tox=10^10#not used here
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
  for(vcost in 1:length(prot_cost_set)){
    prot_cost=prot_cost_set[vcost]
    for (p in 1:length(eta_set)){
      print(paste(etot,vcost,p))
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(vcost-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(vcost-1)]]<-max(Pop_DA)
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
 #            image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

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
Protein_cost<-c("Low","Moderate","High")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
for(cost in 1:length(prot_cost_set)){
  Cost_i<-as.character(Protein_cost[cost])
  print(Cost_i)
  Conc1_opt_DA_eff_E045E1comb[[etot]][[Cost_i]]<-c()
  for (p in 1:length(eta_set)){
    i=2
    while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i,i]==-100 && i<N_reso){
      i=i+1
      if(i==N_reso){
        Conc1_opt_DA_eff_E045E1comb[[etot]][[Cost_i]][p]<-10^-10
      }
    }
    S_opt=0
    while (i<N_reso-1){
      if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i,i-1]<0){
        print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i,i])
        S_opt=i
        print(S_opt)
        Conc1_opt_DA_eff_E045E1comb[[etot]][[Cost_i]][p]<-Etot_var[S_opt]
        print(Etot_var[S_opt])
        break
      }
      else if (i==N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(cost-1)]][i,i+1]>0){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Cost_i]][p]=10^-4
          i=i+1
          break
        }
        else{
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Cost_i]][p]<-0
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
Etot0_set_graph<-c("1e-4.5","var")
txt<-c("A","B","C","D")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(cost in 1:length(prot_cost_set)){
    Cost_i<-as.character(Protein_cost[cost])
    if (cost==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Cost_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Cost_i)]])~log10(eta_set),col=cost,pch=14+cost)
    }
  }
  percent=c(0.05,0.1,0.15,0.2)
  for (p in 1:length(percent)){
    abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
  }
  legend("bottomright",title="Protein cost",legend=c("Low","Moderate","High"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.75)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  text(-3.25,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-4,-4.4,txt[1],srt=addtxt$srt,font=2,col="black",cex=1)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-protcost-firstenzconc.jpeg", width = 500*3,height=415*3,res=300,type="cairo")
dev.off()


#2.Adaptive dynamics showing concentration ceilings for different protein costs##SM Figure S3-B
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

##2.Adaptive dynamics for total content with regard to background concentration (basal metabolism) - Figure S3-B

##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.5##set to a moderately high value
kcat=10^2.5#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value where 
PayOff_set<-c(1)
Path_size_set=c(40)

N_reso=150#Resolution used for strategies (in terms of concentration)
Etot_back_set<-c(5*10^-3,5.5*10^-3,6*10^-3)#other protein concentration
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
Tox=10^10#not used here
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
  for(back in 1:length(Etot_back_set)){
    Etot_back<-Etot_back_set[back]
    for (p in 1:length(eta_set)){
      print(paste(etot,back,p))
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(back-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(back-1)]]<-max(Pop_DA)
    }
  }
}


addtxt<-list(l=0.075,h=1.00,txt=c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"),srt = 0,font=2,col="gray")
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
#multiplePlot("Basal met.=","","Deg. rate=","",c("5e-3","5.5e-3","6e-3"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
#             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
#             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
#             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-475.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()

##Concentration of E0 (enzyme after transporter) similar to that downstream 
multiplePlot("Basal met.=","","Deg. rate=","",c("5e-3","5.5e-3","6e-3"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-var.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()


Conc1_opt_DA_eff_E045E1comb<-list()
Bas_metabo<-c("Low","Moderate","High")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
  for(basal in 1:length(Etot_back_set)){
    Basal_i<-as.character(Bas_metabo[basal])
    print(Basal_i)
    Conc1_opt_DA_eff_E045E1comb[[etot]][[Basal_i]]<-c()
    for (p in 1:length(eta_set)){
      i=1
      if(etot==2){
        i=50
      }
      while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Basal_i]][p]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Basal_i]][p]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(basal-1)]][i,i+1]>0){
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Basal_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Basal_i]][p]<-0
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
Etot0_set_graph<-c("1e-4.5","var")
txt<-c("B","C")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(basal in 1:length(Etot_back_set)){
    Basal_i<-as.character(Bas_metabo[basal])
    if (basal==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Basal_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Basal_i)]])~log10(eta_set),col=basal,pch=14+basal)
    }
  }
  Etot_back_set2=c(5e-3,5.5e-3,6e-3)
  for (p in 1:length(Etot_back_set2)){
    abline(h=log10(0.15/(1-0.15)*Etot_back_set2[p]/Path_size),col=p,lty=4)
    abline(h=log10(0.1/(1-0.1)*Etot_back_set2[p]/Path_size),col=p,lty=4)
  }
  legend("bottomright",title="Basal protein content",legend=c("Low","Moderate","High"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.85)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  text(-3,-4.595,"15%-Mb=6mM",srt=addtxt$srt,font=1,col=3,cex=0.85)
  text(-3,-4.635,"15%-Mb=5.5mM",srt=addtxt$srt,font=1,col=2,cex=0.85)
  text(-3,-4.675,"15%-Mb=5mM",srt=addtxt$srt,font=1,col=1,cex=0.85)
  text(-4,-4.4,txt[etot],srt=addtxt$srt,font=2,col="black",cex=1)
  text(-3,-4.795,"10%-Mb=6mM",srt=addtxt$srt,font=1,col=3,cex=0.85)
  text(-3,-4.835,"10%-Mb=5.5mM",srt=addtxt$srt,font=1,col=2,cex=0.85)
  text(-3,-4.875,"10%-Mb=5mM",srt=addtxt$srt,font=1,col=1,cex=0.85)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-basalmet-firstenzconc.jpeg", width = 500*3,height=415*3,res=300,type="cairo")
dev.off()

##3.Adaptive dynamics for total content with regard to path size - Figure S2-C
##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.5##set to a moderately high value
kcat=10^2.5#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value where 
PayOff_set<-c(1)
Path_size_set=c(20,30,40)

N_reso=150#Resolution used for strategies (in terms of concentration)
Etot_back_set<-5.5*10^-3#other protein concentration
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
Tox=10^10#not used here
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
  for(path in 1:length(Path_size_set)){
    Path_size<-Path_size_set[path]
    for (p in 1:length(eta_set)){
      print(paste(etot,path,p))
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(path-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(path-1)]]<-max(Pop_DA)
    }
  }
}


addtxt<-list(l=0.075,h=1.00,txt=c("E","F","G","H","I","J","K","L","M","N","O","P","Q","R"),srt = 0,font=2,col="gray")
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
#multiplePlot("Path size=","","Deg. rate=","",c("20","30","40"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
#             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
#             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
#             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-475.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()

##Concentration of E0 (enzyme after transporter) similar to that downstream 
multiplePlot("Path size=","","Deg. rate=","",c("20","30","40"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-var.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()


Conc1_opt_DA_eff_E045E1comb<-list()
PathwaySize<-c("20","30","40")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
  for(psize in 1:length(Path_size_set)){
    Psize_i<-as.character(PathwaySize[psize])
    print(Psize_i)
    Conc1_opt_DA_eff_E045E1comb[[etot]][[Psize_i]]<-c()
    for (p in 1:length(eta_set)){
      i=20
      while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Psize_i]][p]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Psize_i]][p]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(psize-1)]][i,i+1]>0){
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Psize_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Psize_i]][p]<-0
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
Etot0_set_graph<-c("1e-4.5","var")
txt<-c("C","F")
addtxt<-list(l=-1.6,h=-2.8,txt=c("E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(psize in 1:length(Path_size_set)){
    Psize_i<-as.character(PathwaySize[psize])
    if (psize==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Psize_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.075),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Psize_i)]])~log10(eta_set),col=psize,pch=14+psize)
    }
  }
  size=c(20,30,40)
  for (p in 1:length(size)){
    abline(h=log10(0.2/(1-0.2)*Etot_back/size[p]),col=p,lty=4)
    abline(h=log10(0.1/(1-0.1)*Etot_back/size[p]),col=p,lty=4)
  }
  legend("bottomright",title="Pathway size",legend=c("20","30","40"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.85)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  #text(-3.25,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-1.25,-4.4,"20%-40R",srt=addtxt$srt,font=1,col=3,cex=0.65)
  text(-1.25,-4.275,"20%-30R",srt=addtxt$srt,font=1,col=2,cex=0.65)
  text(-1.25,-4.1,"20%-20R",srt=addtxt$srt,font=1,col=1,cex=0.65)
  text(-4,-4.075,txt[etot],srt=addtxt$srt,font=2,col="black",cex=1)
  text(-3.25,-4.8,"10%-40R",srt=addtxt$srt,font=1,col=3,cex=0.65)
  text(-3.25,-4.675,"10%-30R",srt=addtxt$srt,font=1,col=2,cex=0.65)
  text(-3.25,-4.5,"10%-20R",srt=addtxt$srt,font=1,col=1,cex=0.65)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-pathsize-firstenzconc.jpeg", , width = 500*3,height=415*3,res=300,type="cairo")
dev.off()

##4.Adaptive dynamics for total content with regard to pathway yield

##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.5##set to a moderately high value
kcat=10^2.5#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value where 
PayOff_set<-c(2.5,5,10)
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
Tox=10^10#not used here
eta_set=10^seq(-4,-1,0.5)
P_p2_set<-c(0)#not used here

#Pathway parameter set
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
  for(po in 1:length(PayOff_set)){
    PO1<-PayOff_set[po]
    PO2<-PayOff_set[po]
    for (p in 1:length(eta_set)){
      print(paste(etot,po,p))
      eta<-eta_set[p]
      P_p2=P_p2_set[1]
      eta_p=eta+P_p2*S_c/V_c
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(po-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(po-1)]]<-max(Pop_DA)
    }
  }
}


addtxt<-list(l=0.075,h=1.00,txt=c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"),srt = 0,font=2,col="gray")
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
#multiplePlot("Path yield=","","Deg. rate=","",c("1:5","1:10","1:20"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
#             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
#             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
#             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-475.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()

##Concentration of E0 (enzyme after transporter) similar to that downstream 
multiplePlot("Path yield=","","Deg. rate=","",c("1:5","1:10","1:20"),c("1e-4","1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-var.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()


Conc1_opt_DA_eff_E045E1comb<-list()
PathYield<-c("1:5","1:10","1:20")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
  for(yield in 1:length(PayOff_set)){
    PayOff_i<-as.character(PathYield[yield])
    print(PayOff_i)
    Conc1_opt_DA_eff_E045E1comb[[etot]][[PayOff_i]]<-c()
    for (p in 1:length(eta_set)){
      i=1
      if(etot==2){
        i=50
      }
      while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[PayOff_i]][p]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_eff_E045E1comb[[etot]][[PayOff_i]][p]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(yield-1)]][i,i+1]>0){
            Conc1_opt_DA_eff_E045E1comb[[etot]][[PayOff_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_eff_E045E1comb[[etot]][[PayOff_i]][p]<-0
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
Etot0_set_graph<-c("1e-4.5","var")
txt<-c("D","H")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,1),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(yield in 1:length(PayOff_set)){
    PayOff_i<-as.character(PathYield[yield])
    if (yield==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(PayOff_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)")#,main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(PayOff_i)]])~log10(eta_set),col=yield,pch=14+yield)
    }
  }
  percent=c(0.05,0.1,0.15,0.2)
  for (p in 1:length(percent)){
    abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
  }
  legend("bottomright",title="Pathway yield",legend=c("1G:5E","1G:10E","1G:20E"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.75)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  text(-3.25,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-4,-4.4,txt[etot],srt=addtxt$srt,font=2,col="black",cex=1)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-pathyield-firstenzconc.jpeg", , width = 500*3,height=415*3,res=300,type="cairo")
dev.off()

##5.Adaptive dynamics for total content with regard to environment replenishment

##Constants

VTm=1e-3#M/s
KT=1e-2#M
kf=10^6.5##set to a moderately high value
kcat=10^2.5#/s##set to a moderately high value
kr=kcat#
prot_cost=10^-2.5#Mmet/Menz##set to an average value where
PayOff_set<-c(5)
Path_size_set=c(40)

N_reso=150#Resolution used for strategies (in terms of concentration)
Etot_back<-5.5*10^-3#other protein concentration
E_basal<-3*10^-3#scaling factor

#Adaptive dynamics constants

Phi_eq=10^-4
r_c=1e-5#dm
S_c=4*pi*r_c^2
V_c=4/3*pi*r_c^3#dm^3
beta_set=c(10^-3,10^-4,10^-5)
alpha_set=c(10^-3,10^-4,10^-5)
dl=1e-4#dm
V_env=dl^3

##Parameters

Etot0_set<-c(10^-4.5,0)
Tox=10^10#not used here
eta_set=10^seq(-3.5,-1,0.5)
P_p2_set<-c(0)#not used here

#Pathway parameter set
Path_size=Path_size_set[1]
PO1<-PayOff_set[1]
PO2<-PayOff_set[1]
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
  for(rep in 1:length(beta_set)){
    alpha<-alpha_set[rep]
    beta<-beta_set[rep]
    for (p in 1:length(eta_set)){
      print(paste(etot,rep,p))
      eta<-eta_set[p]
      P_p2=P_p2_set[1]
      eta_p=eta+P_p2*S_c/V_c
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
      list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep-1)]]<-fit_inv
      tab_N_eq[[p+length(eta_set)*(rep-1)]]<-max(Pop_DA)
    }
  }
}


addtxt<-list(l=0.075,h=1.00,txt=c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"),srt = 0,font=2,col="gray")
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
multiplePlot("Alpha=","","Deg. rate=","",c("1e-3","1e-4","1e-5"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[1]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")

#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-475.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()

##Concentration of E0 (enzyme after transporter) similar to that downstream 
multiplePlot("Alpha=","","Deg. rate=","",c("1e-3","1e-4","1e-5"),c("1e-3.5","1e-3","1e-2.5","1e-2","1e-1.5","1e-1")
             ,ncol=ncol,log10Etot1_var,log10Etot1_var,list_fit_inv_S40_noperm[[2]],
             abs="resident (logEtot)",ord="mutant (logEtot)",scale=c(-80,80),lev=c(0),palette=pal,cextext=2,
             image=TRUE,pcex=0.75,subcex=0.75,labcex=1.5,axcex=0.65,globcex=0.5,legcex=1,contourlab=TRUE,meth="edge",contcex=0.5,colorkey="COMMON")


#setwd(dir="~")
#setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
#dev.print(device = jpeg, file = "PIP-enzeff-firstenzconc-var.jpeg", , width = 1600*2,height=740*2,res=200,type="cairo")
#dev.off()


Conc1_opt_DA_eff_E045E1comb<-list()
Replenishment_rate<-c("1e-3","1e-4","1e-5")
#Conc1_opt_DA_E0E1comb[["2*10^-3"]]<-c()
for (etot in 1:length(Etot0_set)){
  Conc1_opt_DA_eff_E045E1comb[[etot]]<-list()
  for(rep_r in 1:length(beta_set)){
    Rep_i<-as.character(Replenishment_rate[rep_r])
    print(Rep_i)
    Conc1_opt_DA_eff_E045E1comb[[etot]][[Rep_i]]<-c()
    for (p in 1:length(eta_set)){
      i=1
      if(etot==2){
        i=50
      }
      while(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i,i]==-100 && i<N_reso){
        i=i+1
        if(i==N_reso){
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Rep_i]][p]<-10^-10
        }
      }
      S_opt=0
      while (i<N_reso-1){
        if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i,i+1]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i,i-1]<0){
          print(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i,i])
          S_opt=i
          print(S_opt)
          Conc1_opt_DA_eff_E045E1comb[[etot]][[Rep_i]][p]<-Etot_var[S_opt]
          print(Etot_var[S_opt])
          break
        }
        else if (i==N_reso-1){
          if(list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i+1,i]<0 && list_fit_inv_S40_noperm[[etot]][[p+length(eta_set)*(rep_r-1)]][i,i+1]>0){
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Rep_i]][p]=10^-4
            i=i+1
            break
          }
          else{
            Conc1_opt_DA_eff_E045E1comb[[etot]][[Rep_i]][p]<-0
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
Etot0_set_graph<-c("1e-4.5","var")
txt<-c("I","J")
addtxt<-list(l=-1.6,h=-2.8,txt=c("A","B","C","D","E","F","G","H","I"),srt = 0,font=2,col="black")
par(mfrow=c(1,2),mai=c(1,1.125,0.5,0.875))
for(etot in 1:length(Etot0_set)){
  for(rep_r in 1:length(beta_set)){
    Rep_i<-as.character(Replenishment_rate[rep_r])
    if (rep_r==1){
      plot(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Rep_i)]])~log10(eta_set),col=1,pch=15,ylim=c(-5.2,-4.4),xlab="log10(Degradation rate)",ylab="log10(Enzyme concentration)",main=paste("[E0]=",Etot0_set_graph[etot]))
    }
    else{
      points(log10(Conc1_opt_DA_eff_E045E1comb[[etot]][[as.character(Rep_i)]])~log10(eta_set),col=rep_r,pch=14+rep_r)
    }
  }
  percent=c(0.05,0.1,0.15,0.2)
  for (p in 1:length(percent)){
    abline(h=log10(percent[p]/(1-percent[p])*Etot_back/Path_size),col="dark grey",lty=4)
  }
  legend("bottomright",title="Pathway yield",legend=c("1e-3","1e-4","1e-5"),bty="y",pch=c(15:18),col=c(1:4),ncol=1,cex=0.75)
  #text(addtxt$l,addtxt$h,addtxt$txt[1],srt=addtxt$srt,font=addtxt$font,col=addtxt$col,cex=1)
  text(-3.25,-5.12,"5%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.79,"10%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.59,"15%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.25,-4.435,"20%",srt=addtxt$srt,font=1,col="black",cex=1)
  text(-3.5,-4.4,txt[etot],srt=addtxt$srt,font=2,col="black",cex=1)
}

setwd(dir="~")
setwd(dir="/Users/florianlabourel/Desktop/Ongoing-projects/cross-feeding/Proteome allocation/Draft/SM/Figures-February")
dev.print(device = jpeg, file = "res-pathyield-firstenzconc.jpeg", , width = 1075*3,height=460*3,res=300,type="cairo")
dev.off()



