##Function to plot multiple images
library(lattice)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(plot.matrix)
library(fields)
library(plot3D)

##Function to plot multiple images
multiplePlot<-function(nvar1,unitV1,nvar2,unitV2,v1,v2,ncol,xlim,ylim,xyZ,tit,
                       abs,ord,palette,scale,ncont,lev,xyData,rego,TEXT_to_Add,
                       col_data,image,ncat,pcex,labcex,legtext,legpos,globcex,legcex,legtitle,
                       sub,subcex,axcex,cextext,pch_dat,meth,colorkey,contourlab,addline,contcex){
  len1<-length(v1)
  len2<-length(v2)
  if(length(palette)==1){
    if(colorkey=="COMMON"){
      t<-c()
      for(i in 1:len1){
        if(i<2){
          t<-c(t,(i-1)*(len2)+seq(1,len2),len1*len2+1)
        }
        else{
          t<-c(t,(i-1)*len2+seq(1,len2),len1*len2+1)
        }
      }
      layout(matrix(t, len1, len2+1, byrow = TRUE), 
             widths=c(rep(1,len2),0.3), heights=c(rep(1,len1)))
    }
    else{
      par(mfrow=c(len1,len2))
    }
    pchlegend<-c()
    collegend<-c()
    for (l in 1:len1){
      for (s in 1:len2){
        if(colorkey==TRUE){
          if(l==1 && s==1){
            if(len1==1 && len2==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
                }
                else{
                  par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
                }
              }
              else{
                par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
              }
            }
            else{
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.75,1,0.5,1),cex=globcex)
                }
                else{
                  par(mai = c(0.75,1,0.5,1),cex=globcex)
                }
              }#mar=c(6,6,2,1.5),
              else{
                par(mai = c(0.75,1,1.5,1),cex=globcex)
              }
            }
          }
          else if(l==1){
            if(missing(tit)){
              par(mai = c(0.75, 0.75, 0.5, 1.25),cex=globcex)
            }
            else{
              par(mai = c(0.75, 0.75, 1.5, 1.25),cex=globcex)
            }
          }
          else if(s==1){
            if(missing(tit)){
              par(mai = c(0.75, 1, 0.5, 1),cex=globcex)
            }
            else{
              par(mai = c(0.75, 1, 1.5, 1),cex=globcex)
            }
          }
          else{
            if(missing(tit)){
              par(mai = c(0.75, 0.75, 0.5, 1.25),cex=globcex)
            }
            else{
              par(mai = c(0.75, 0.75, 1.5, 1.25),cex=globcex)
            }
          }
        }
        
        else if(colorkey==FALSE){
          if(l==1 && s==1){
            if(len1==1 && len2==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(legtitle)){
                    par(mai = c(0.75,1.25,0.75,0.75),cex=globcex)
                  }
                  else{
                    par(mai = c(0.75,0.75,0.75,1.5),cex=globcex)
                  }
                }
                else{
                  par(mai = c(0.75,1,0.5,0.75),cex=globcex)
                }
              }
              else{
                par(mai = c(0.75,1,0.5,0.75),cex=globcex)
              }
            }
            else{
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.75,1,0.5,0.75),cex=globcex)
                }
                else{
                  par(mai = c(0.75,1,0.5,0.75),cex=globcex)
                }
              }#mar=c(6,6,2,1.5),
              else{
                par(mai = c(0.75,1,1.5,0.75),cex=globcex)
              }
            }
          }
          else if(l==1){
            if(missing(tit)){
              par(mai = c(0.75, 1, 0.5, 0.75),cex=globcex)
            }
            else{
              par(mai = c(0.75, 1, 1.5, 0.75),cex=globcex)
            }
          }
          else if(s==1){
            if(missing(tit)){
              par(mai = c(0.75, 1, 0.5, 0.75),cex=globcex)
            }
            else{
              par(mai = c(0.75, 1, 1.5, 0.75),cex=globcex)
            }
          }
          else{
            if(missing(tit)){
              par(mai = c(0.75, 1, 0.5,0.75),cex=globcex)
            }
            else{
              par(mai = c(0.75, 1, 1.5, 0.75),cex=globcex)
            }
          }
        }
        
        else{
          if(l==1 && s==1){
            if(len1==1 && len2==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
                }
                else{
                  par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
                }
              }
              else{
                par(mai = c(0.75,0.75,0.75,0.75),cex=globcex)
              }
            }
            else if(len2==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.5,0.75,0.5,0.25),cex=globcex)
                }
                else{
                  par(mai = c(0.5,0.75,0.5,0.25),cex=globcex)
                }
              }
              else{
                par(mai = c(0.75,1,1.5,0.05),cex=globcex)
              }
            }
            else{
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.6,0.75,0.25,0.25),cex=globcex)
                }
                else{
                  par(mai = c(0.6,0.75,0.25,0.25),cex=globcex)
                }
              }#mar=c(6,6,2,1.5),
              else{
                par(mai = c(0.6,0.75,1.5,0.25),cex=globcex)
              }
            }
          }
          else if(l==1){
            if(missing(tit)){
              par(mai = c(0.6, 0.75, 0.25, 0.25),cex=globcex)
            }
            else{
              par(mai = c(0.6, 0.75, 1.5, 0.25),cex=globcex)
            }
          }
          else if(s==1){
            if(len2==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  par(mai = c(0.5,0.75,0.5,0.25),cex=globcex)
                }
                else{
                  par(mai = c(0.5,0.75,0.5,0.25),cex=globcex)
                }
              }
              else{
                par(mai = c(0.75,1,1.5,0.05),cex=globcex)
              }
            }
          else{
              if(missing(tit)){
                par(mai = c(0.6, 0.75, 0.25, 0.25),cex=globcex)
              }
              else{
                par(mai = c(0.6, 0.75, 1.5, 0.25),cex=globcex)
              }
          }
        }
        else{
            if(missing(tit)){
              par(mai = c(0.6, 0.75, 0.25, 0.25),cex=globcex)
            }
            else{
              par(mai = c(0.6, 0.75, 1.5, 0.25),cex=globcex)
            }
          }
        }
        if(image==TRUE){
          if(missing(scale)){
            Maximal<-max(xyZ[[s+(l-1)*len2]])
            Minimal<-min(xyZ[[s+(l-1)*len2]])
            if(Minimal<0){
              newz.na <- 0-(Maximal)/ncol
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
            #ncol=as.integer((Maximal-Minimal)/Minimal)
            #palette<-colorRampPalette(palette)(ncol-1)
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=1),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=1),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                  }
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),
                             legend.shrink=1,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),
                               legend.shrink=1,legend.cex=legcex)
                  }
                }
              }
            }
          }
          else if(length(scale)==1){
            if(scale=="AUTO"){
              Maximal<-max(xyZ[[s+(l-1)*len2]])
              Minimal<-min(xyZ[[s+(l-1)*len2]])
              if(Minimal<0){
                newz.na <- 0-(Maximal)/ncol
                xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<0)] <- newz.na
              #ncol=as.integer((Maximal-Minimal)/Minimal)
              #palette<-colorRampPalette(palette)(ncol-1)
                if(missing(sub)){
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    if(missing(colorkey) || colorkey==TRUE){
                      image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                    }
                  }
                  else{
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),
                               legend.shrink=1,legend.cex=legcex,
                               cex.main=subcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),
                                 legend.shrink=1,legend.cex=legcex,
                                 cex.main=subcex)
                    }
                  }
                }
                else{
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    }
                  }
                  else{
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2)
                               ,legend.cex=legcex,legend.shrink=1)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),nlevel=ncol,col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2)
                                 ,legend.cex=legcex,legend.shrink=1)
                    }
                  }
                }
              }
              else{
                if(missing(sub)){
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                    }
                  }
                  else{
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2))
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2))
                    }
                  }
                }
                else{
                  if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    }
                  }
                  else{
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol+1,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),
                               legend.cex=legcex,
                               legend.shrink=1)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],zlim=c(Minimal,Maximal),nlevel=ncol+1,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),
                                 legend.cex=legcex,
                                 legend.shrink=1)
                    }
                  }
                }
              }
            }
          }
          else{
            if(length(which(xyZ[[s+(l-1)*len2]][]>scale[2]))>0){
              newz.na <- scale[2]+(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]>scale[2])] <- newz.na
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",legend.cex=legcex)
                    
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],
                               zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                               #breaks=seq(min(scale[1]),max(scale[2])),
                               nlevel=ncol,
                               col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                    
                    
                  }
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],
                               zlim=c(scale[1],(scale[2]-scale[1])/ncol+scale[2]),
                               #breaks=seq(min(scale[1]),max(scale[2])),
                               nlevel=ncol,
                               col=c(palette[[1]],"light gray"),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                    
                  }
                }
              }
            ##mtext()
            }
            else if(length(which(xyZ[[s+(l-1)*len2]][]<scale[1]))>0){
              newz.na <- scale[1]-(scale[2]-scale[1])/ncol # new z for NA
              xyZ[[s+(l-1)*len2]][][which(xyZ[[s+(l-1)*len2]][]<scale[1])] <- newz.na
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,scale[2]),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,scale[2]),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                    
                  }
                }
                  else{
                    if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],
                               zlim=c(newz.na,scale[2]),
                               #breaks=seq(min(scale[1]),max(scale[2])),
                               nlevel=ncol,
                               col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                    }
                    else{
                      image(xyZ[[s+(l-1)*len2]],
                                 zlim=c(newz.na,scale[2]),
                                 #breaks=seq(min(scale[1]),max(scale[2])),
                                 nlevel=ncol,
                                 col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                      
                    }
                  }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,scale[2]),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,scale[2]),col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    
                  }
                    }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                    image.plot(xyZ[[s+(l-1)*len2]],
                             zlim=c(newz.na,scale[2]),
                             #breaks=seq(min(scale[1]),max(scale[2])),
                             nlevel=ncol,
                             col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],
                               zlim=c(newz.na,scale[2]),
                               #breaks=seq(min(scale[1]),max(scale[2])),
                               nlevel=ncol,
                               col=c("light gray",palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                    
                  }
                    }
              }
            }
            else{
              if(missing(sub)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",legend.cex=legcex)
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",axis.args=list(cex.axis=2),legend.shrink=1,cex.main=subcex,legend.cex=legcex)
                  }
                }
              }
              else{
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(newz.na,Maximal),col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",legend.cex=legcex)
                    
                  }
                }
                else{
                  if(missing(colorkey) || colorkey==TRUE){
                  image.plot(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,
                             axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                  }
                  else{
                    image(xyZ[[s+(l-1)*len2]],zlim=c(scale[1],scale[2]),nlevel=ncol,col=c(palette[[1]]),xaxt="n",yaxt="n",main=sub[s+(l-1)*len2],cex.main=subcex,col.main="black",font.main=1,
                               axis.args=list(cex.axis=2),legend.shrink=1,legend.cex=legcex)
                  }
                }
              }
            }
          }
          if(missing(ncont)){
          }
          else{
            if(missing(meth)){
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont, add = TRUE,labcex=contcex,method="edge",drawlabels = contourlab)
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],nlevels=ncont, add = TRUE,labcex=contcex,method=meth,drawlabels = contourlab)
            }
          }
          if(missing(lev)){
          }
          else{
            if(missing(meth)){
            contour(xyZ[[s+(l-1)*len2]],levels=lev, add = TRUE,labcex=contcex,method="edge",drawlabels = contourlab)
            }
            else{
              contour(xyZ[[s+(l-1)*len2]],levels=lev, add = TRUE,labcex=contcex,method=meth,drawlabels = contourlab)
            }
          }
        }
        else{
          par(mfrow=c(1,1))
          if(missing(ncont)){
          }
          else{
            if(s==1 & l==1){
              if(missing(meth)){
                contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,xaxt="n",yaxt="n",labcex=contcex,col=1,lty=4,method = "edge",drawlabels = contourlab)
                collegend[1]=2
              }
              else{
                contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,xaxt="n",yaxt="n",labcex=contcex,col=1,lty=4,method = meth,drawlabels = contourlab)
                collegend[1]=2
              }
            }
            else{
              if(missing(meth)){
                contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,add=TRUE,xaxt="n",yaxt="n",labcex=contcex,col=s+(l-1)*len2,lty=4,method = "edge",drawlabels = contourlab)
                collegend[s+(l-1)*len2]=1+s+(l-1)*len2
              }
              else{
                contour(xyZ[[s+(l-1)*len2]],nlevels=ncont,add=TRUE,xaxt="n",yaxt="n",labcex=contcex,col=s+(l-1)*len2,lty=4,method = meth,drawlabels = contourlab)
                collegend[s+(l-1)*len2]=1+s+(l-1)*len2
              }
            }
          }
          if(missing(lev)){
          }
          else{
            if(s==1 & l==1){
              if(missing(meth)){
                contour(xyZ[[s+(l-1)*len2]],levels=lev,xaxt="n",yaxt="n",labcex=contcex,col=2,lty=4,method = "simple",drawlabels = contourlab)
                collegend[1]=2
              }
              else{
                contour(xyZ[[s+(l-1)*len2]],levels=lev,xaxt="n",yaxt="n",labcex=contcex,col=2,lty=4,method = meth,drawlabels = contourlab)
                collegend[1]=2
              }
            }
            else{
              if(missing(meth)){
                contour(xyZ[[s+(l-1)*len2]],levels=lev,add=TRUE,xaxt="n",yaxt="n",labcex=contcex,col=1+s+(l-1)*len2,lty=4,method = "edge",drawlabels = contourlab)
                collegend[s+(l-1)*len2]=1+s+(l-1)*len2
              }
              else{
                contour(xyZ[[s+(l-1)*len2]],levels=lev,add=TRUE,xaxt="n",yaxt="n",labcex=contcex,col=1+s+(l-1)*len2,lty=4,method = meth,drawlabels = contourlab)
                collegend[s+(l-1)*len2]=1+s+(l-1)*len2
              }
            }
          }
        }
        
        if(missing(addline)){
        }
        else{
          abline(a=addline[1],b=addline[2],col=addline[3],lty=addline[4],lwd=addline[5])
        }
      
        if(missing(TEXT_to_Add)){
        }
        else{
          if(missing(cextext)){
            text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
          }
          else{
            text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=cextext)
          }
        }
        
        if(missing(rego)){
        }
        else{
          abline(coef=c(((rego$coefficients[1]-rego$coefficients[2]
                          *(0-xlim[1]))-ylim[1])/(ylim[2]-ylim[1]),
                        (rego$coefficients[2])),
                lwd=1,col="white")
        }
    
        if(missing(xyData)){
        }
        else{
          if(missing(col_data)){
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=i+1,pch=20)
                  pchlegend[1]=20
                  collegend[i]=1
                }
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=i+1,pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=1
                }
              }
            }
            else{
              #for (i in 1:(length(xyData)/ncat)){
                for (j in 1:length(xyData)){
                  #for(i in 1:length(xyData[[j]][[1]])){
                    if(missing(pcex)){
                      points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                          (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                          col=j,pch=j)
                      pchlegend[j]=j
                      collegend[j]=j
                     }
                    else{
                      points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                          (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                          col=j,pch=j,cex=pcex)
                      pchlegend[j]=j
                      collegend[j]=j
                    }
                  #}
              }
            }
          }
          else{
            if(missing(ncat)){
              for (i in 1:length(xyData)){
                if(missing(pcex)){
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                     (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                     col=col_data[i],pch=20)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
              
                else{
                  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                       (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                       col=col_data[i],pch=20,cex=pcex)
                  pchlegend[1]=20
                  collegend[i]=col_data[i]
                }
              }
            }
            else{
              #for (i in 1:(length(xyData)/ncat)){
                for (j in 1:length(xyData)){
                    if(missing(pcex)){
                      if(missing(pch_dat)){
                        points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                               (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                               col=col_data[j],pch=j)
                        pchlegend[j]=j
                        collegend[j]=col_data[j]
                      }
                     else{
                       points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                              (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                              col=col_data[j],pch=pch_dat[j])
                       pchlegend[j]=pch_dat[j]
                       collegend[j]=col_data[j]
                     }
                    }
                    else{
                      if(missing(pch_dat)){
                      points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                          (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                          col=col_data[j],pch=j,cex=pcex)
                      pchlegend[j]=j
                      collegend[j]=col_data[j]
                      }
                      else{
                        points((xyData[[j]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
                               (xyData[[j]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
                               col=col_data[j],pch=pch_dat[j],cex=pcex)
                        pchlegend[j]=pch_dat[j]
                        collegend[j]=col_data[j]
                      }
                  }
                #}
              }
            }
          }
        #print(pchlegend[1])
          if(missing(TEXT_to_Add)){
          
          }
          else{
            if(missing(cextext)){
              text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=1)
            }
            else{
              text(TEXT_to_Add$l,TEXT_to_Add$h,TEXT_to_Add$txt[s+(l-1)*len2],srt=TEXT_to_Add$srt,font=TEXT_to_Add$font,col=TEXT_to_Add$col,cex=cextext)
            }
          }
        }
        if(missing(legtext)){
        }
        else{
          if(image==TRUE){
            if(is.character(legpos)){
              if(missing(legtitle)){
                legend(legpos,legend=legtext,col=collegend,pch=pchlegend, cex=pcex*legcex)
              }
              else{
                legend(legpos,legend=legtext,col=collegend,pch=pchlegend,title=legtitle, cex=pcex*legcex)
              }
            }
            else{
              if(missing(legtitle)){
              legend(legpos[1],legpos[2],xpd=NA,legend=legtext,col=collegend,pch=pchlegend, cex=pcex*legcex)
              }
              else{
                legend(legpos[1],legpos[2],xpd=NA,legend=legtext,col=collegend,pch=pchlegend,title=legtitle, cex=pcex*legcex)
              }
            }
          }
          else{
            if(missing(legtitle)){
            legend(legpos,legend=legtext,col=collegend,lty=6, cex=pcex*legcex)
            }
            else{
              legend(legpos,legend=legtext,col=collegend,lty=6, cex=pcex*legcex,title=legtitle)
            }
          }
        }
        if(image==TRUE){
          if(missing(axcex)){
            if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
            axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,tck=-0.03)
            mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=2.25*2,side=1,cex=1.25)
            axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,tck=-0.03)
            mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=1.75*2,side=2,cex=1.25)
            }
            else{
              if(l==len1){
                axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,tck=-0.03,line=0)
                mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=4*axcex/3,side=1,cex=1.25)
              }
              if(s==1){
                axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,tck=-0.03,line=0)
                mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=2*axcex/3,side=2,cex=1.25)
              }
            }
          }
          else{
            if(colorkey==TRUE){
                  if(l==len1){
                  axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,tck=-0.03,line=0)
                  mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=2*axcex/3,side=1,cex=axcex)
                  }
                  if(s==1){
                  axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,tck=-0.03,line=0)
                  mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=2*axcex/3,side=2,cex=axcex)
                  }
            }
            else{
              if(l==len1){
                axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,tck=-0.03,line=0)
                if(missing(legtitle)){
                  mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=4*axcex/3,side=1,cex=axcex)
                }
                else{
                  mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=6*axcex/3,side=1,cex=axcex)
                }
              }
              if(s==1){
                axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,tck=-0.03,line=0)
                mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=4*axcex/3,side=2,cex=axcex)
              }
            }
          }
        }
        else{
          if(s==1 & l==1){
            if(missing(axcex)){
              if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                if(colorkey==TRUE){
              axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,cex.axis=1.25)
              mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=4*1.25/3,side=1,cex=1.25)
              axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,cex.axis=1.25)
              mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=2*1.25/3,side=2,cex=1.25)
                }
              }
              else{
                if(colorkey==TRUE){
                axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=1.25)
                axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=1.25)
                }
              }
            }
            else{
              if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                if(colorkey==TRUE){
              axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=FALSE,cex.axis=axcex,line=0)
              mtext(text=seq(xlim[1],xlim[2],xlim[3]),at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),line=4*axcex/3,side=1,cex=axcex)
              axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=FALSE,cex.axis=axcex,line=0)
              mtext(text=seq(ylim[1],ylim[2],ylim[3]),at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),line=3*axcex/3,side=2,cex=axcex)
                }
              }
              else{
                    if(l==len1){
                      axis(1, at=seq(0,1,xlim[3]/abs(xlim[2]-xlim[1])),labels=seq(xlim[1],xlim[2],xlim[3]),cex.axis=axcex)
                    }
                    if(s==len2){
                      axis(2, at=seq(0,1,ylim[3]/abs(ylim[2]-ylim[1])),labels=seq(ylim[1],ylim[2],ylim[3]),cex.axis=axcex)
                    }
                
              }
            }
          }
        }
        if(missing(labcex)){
          if(image==TRUE){
            if(missing(tit)){
              if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                title(xlab=abs,line=6,cex.lab=1.25)
                title(ylab=ord,line=5,cex.lab=1.25)
              }
              else{
                title(xlab=abs,ylab=ord,line=2.5)
              }
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=2)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                if(colorkey==TRUE){
                  title(xlab=abs,ylab=ord,line=2.5)
                }
                else{
                  title(xlab=abs,line=1.5,cex.lab=1.25)
                  title(ylab=ord,line=3.5,cex.lab=1.25)
                }
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord)
              }
            }
          }
        }
        else{
          if(image==TRUE){
            if(missing(tit)){
              if(colorkey==TRUE){
                if(s==1){
                title(ylab=ord,line=2.25,cex.lab=labcex)
                }
                if(l==len1){
                title(xlab=abs,line=2.25,cex.lab=labcex)
                }
              }
              else if(colorkey==FALSE){
                  if(s==1){
                    title(ylab=ord,line=4.5,cex.lab=labcex)
                  }
                  if(l==len1){
                    title(xlab=abs,line=4,cex.lab=labcex)
                  }
              }
              else{
                if(s==1){
                  if(missing(legtitle)){
                    title(ylab=ord,line=(labcex+axcex)*1.1,cex.lab=labcex)
                  }
                  else{
                    title(ylab=ord,line=(labcex+axcex)*1.1,cex.lab=labcex)
                  }
                }
                if(l==len1){
                  if(missing(legtitle)){
                    title(xlab=abs,line=(labcex+axcex)*1.1,cex.lab=labcex)
                  }
                else{
                  title(xlab=abs,line=(labcex+axcex)*1.1,cex.lab=labcex)
                }
                }
              }
            }
            else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
              title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
            }
            else{
              title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
            }
          }
          else{
            if(s==1 & l==1){
              if(missing(tit)){
                if(pal[[1]][1]==pal[[1]][length(pal[[1]])]){
                  title(ylab=ord,line=5,cex.lab=labcex)
                  title(xlab=abs,line=6,cex.lab=labcex)
                }
                else{
                  title(xlab=abs,ylab=ord,line=2.5,cex.lab=labcex)
                }
              }
              else if (grepl(pattern="\n",tit)){## the ideal would be to find the carriage return character
                title(main=tit,outer=T,line=-3,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
              }
              else{
                title(main=tit,outer=T,line=-1.5,cex.main=1.25,font=2,xlab=abs,ylab=ord,cex.lab=labcex)
              }
            }
          }
        }
      
        if (len2==1){
          if(image==TRUE){
            if(colorkey=="COMMON"){
              mtext(paste(nvar1,v1[l],unitV1),side=3,line=1,cex=1,font=2)
            }
            else{
              mtext(paste(nvar1,v1[l],unitV1),side=3,line=1,cex=1,font=2)
            }
          }
          else{
          }
        }
        else{
          if(l==len1){
            if(image==TRUE){
              if(colorkey=="COMMON"){
                mtext(paste(nvar2,v2[s],unitV2),side=1,line=4.5,cex=1,font=2)
              }
              else{
                mtext(paste(nvar2,v2[s],unitV2),side=1,line=4.5,cex=1,font=2)
              }
            }
            else{
            }
          }
          if(s==1){
            if(image==TRUE){
              if(colorkey=="COMMON"){
              mtext(paste(nvar1,v1[l],unitV1),side=2,line=6,cex=1,font=2)
              }
              else{
                mtext(paste(nvar1,v1[l],unitV1),side=2,line=6,cex=1,font=2)
              }
            }
            else{
            }
          }
        }
      }
    }
    if(colorkey=="COMMON"){
      colkey(clim=c(0,1),col=palette[[1]],side=2,length=1,width=5,cex.axis=1.75,add=FALSE)
    }
  }
}

#Defining colors
#ncol=64
#pal<-colorRampPalette(c("hot pink","red","green"))(ncol)

