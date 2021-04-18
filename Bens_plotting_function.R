
Item_Parameters<-readRDS("Item_Parameters.rds")
PF8c <- c("PFA9","PFA15","PFA16r1","PFA21","PFA55","PFA23","PFA1","PFB50")
ipar<-Item_Parameters$Physical_Function[PF8c,-1]

Data<-matrix(sample(c(1:5),size=8*500,replace=T),ncol=8)

#comment things - data author, credits to seung for functs, formats of parameters going in and out
###BS wants to share with other analysts who are working sf data or calibrating
##Also make a test case to go with the script

BFF<-function(ipar, Data_to_score=NA, maxCat=5)
{
        
        
        marginal.reliability<-function(ipar,model=1,maxCat=5,D=1,minTheta=-4.0,maxTheta=4.0,inc=0.1) {
          #computes marginal reliability based on IRT item parameter estimates based on GRM or GPCM
          #ipar = a data frame containing item parameters without a header in the order of 
          #		a, cb1-cb[maxCat], ncat, where maxCat is defined below
          #maxCat = maximum number of response categories across items (default: maxCat=5)
          #model = IRT model by which item parameters have been calibrated, 1 for GRM or 2 for GPCM (default: model=1)
          
          NCAT<-ipar[,"NCAT"];
          DISC<-ipar[,"a"];
          CB<-ipar[paste("cb",1:(maxCat-1),sep="")];
          
          theta<-seq(minTheta,maxTheta,inc)
          nq=length(theta)
          ni<-nrow(ipar)
          
          prep.prob.info<-function(){
            pp<-array(0,c(nq,ni,maxCat))
            matrix.info<-matrix(0,nq,ni)
            if (model==1) {
              for (i in 1:ni) {
                ps<-matrix(0,nq,NCAT[i]+1)
                ps[,1]<-1
                ps[,NCAT[i]+1]<-0
                for (k in 1:(NCAT[i]-1)) {
                  ps[,k+1]<-1/(1+exp(-D*DISC[i]*(theta-CB[i,k])))
                }
                pp[,i,1]<-1-ps[,1]
                pp[,i,NCAT[i]]<-ps[,NCAT[i]]
                for (k in 1:NCAT[i]) {
                  pp[,i,k]=ps[,k]-ps[,k+1]
                  matrix.info[,i]<-matrix.info[,i]+(D*DISC[i]*(ps[,k]*(1-ps[,k])-ps[,k+1]*(1-ps[,k+1])))^2/pp[,i,k]
                }
              }
            } else if (model==2) {
              for (i in 1:ni) {
                cb<-unlist(CB[i,])
                cb<-c(0,cb)
                zz<-matrix(0,nq,NCAT[i])
                sdsum<-0
                den<-rep(0,nq)
                
                for (k in 1:NCAT[i]) {
                  sdsum<-sdsum+cb[k]
                  zz[,k]<-exp(D*DISC[i]*(k*theta-sdsum))
                  den<-den+zz[,k]
                }
                AX<-rep(0,nq); BX<-rep(0,nq)
                for (k in 1:NCAT[i]) {
                  pp[,i,k]<-zz[,k]/den
                  AX<-AX+k^2*pp[,i,k]
                  BX<-BX+k*pp[,i,k]
                }
                matrix.info[,i]<-D^2*DISC[i]^2*(AX-BX^2)
              }
            }
            list(pp=pp,matrix.info=matrix.info)
          }
          matrix.prob.info<-prep.prob.info()
          matrix.info<-matrix.prob.info$matrix.info
          
          scale.info<-rowSums(matrix.info)+1
          scale.SE<-1/sqrt(scale.info)
          g.theta<-dnorm(theta)
          rel<-1-sum(scale.SE^2*g.theta)/sum(g.theta)
          return (list(ipar=ipar,theta=theta,item.info=matrix.info,scale.info=scale.info,marginal.reliability=rel))
          #return (list(marginal.reliability=rel))
        }
        
        
        thetaSE.eap<-function(ipar,resp.data,maxCat=5,model=1,minTheta=-4.0,maxTheta=4.0,inc=0.1,prior.dist=1,prior.mean=0.0,prior.sd=1.0,D=1.0) {
          ni<-nrow(ipar); #number of items
          nExaminees<<-nrow(resp.data); #number of examiness
          NCAT<-ipar$NCAT;
          
          theta<-seq(minTheta,maxTheta,inc);
          nq<-length(theta);
          
          if (prior.dist==1) {
            prior<-dnorm((theta-prior.mean)/prior.sd); #normal prior
          } else if (prior.dist==2) {
            prior<-exp((theta-prior.mean)/prior.sd)/(1+exp((theta-prior.mean)/prior.sd))^2; #logistic prior
          } else if (prior.dist==3) {
            prior<-rep(1,nq); #uniform prior
          }
          DISC<-ipar[["a"]];
          CB<-ipar[paste("cb",1:(maxCat-1),sep="")];
          
          prep.prob<-function(){
            pp<-array(0,c(nq,ni,maxCat));
            if (model==1) {
              for (i in 1:ni) {
                ps<-matrix(0,nq,NCAT[i]+1);
                ps[,1]<-1;
                ps[,NCAT[i]+1]<-0;
                for (k in 1:(NCAT[i]-1)) {
                  ps[,k+1]<-1/(1+exp(-D*DISC[i]*(theta-CB[i,k])));
                }
                #pp[,i,1]<-1-ps[,1];
                #pp[,i,NCAT[i]]<-ps[,NCAT[i]];
                for (k in 1:NCAT[i]) {
                  pp[,i,k]=ps[,k]-ps[,k+1];
                }
              }
            } else if (model==2) {
              for (i in 1:ni) {
                cb<-unlist(CB[i,]);
                cb<-c(0,cb);
                zz<-matrix(0,nq,NCAT[i]);
                sdsum<-0;
                den<-rep(0,nq);
                
                for (k in 1:NCAT[i]) {
                  sdsum<-sdsum+cb[k];
                  zz[,k]<-exp(D*DISC[i]*(k*theta-sdsum));
                  den<-den+zz[,k];
                }
                for (k in 1:NCAT[i]) {
                  pp[,i,k]<-zz[,k]/den;
                }
              }
            }
            
            return(pp);
          }
          
          pp<-prep.prob();
          
          calcEAP<-function() {
            posterior<-matrix(rep(prior,nExaminees),nExaminees,nq,byrow=T);
            
            for (i in 1:ni) {
              resp<-matrix(resp.data[,i],nExaminees,1);
              prob<-t(pp[,i,resp]);
              prob[is.na(prob)]<-1.0
              posterior<-posterior*prob;
            }
            EAP<-as.vector(posterior%*%theta/rowSums(posterior));
            SE<-as.vector(sqrt(rowSums(posterior*(matrix(theta,nExaminees,nq,byrow=T)-matrix(EAP,nExaminees,nq))^2)/rowSums(posterior)));
            return(list(theta=EAP,SE=SE))
          }
          return(data.frame(calcEAP()));
        }
        
  
  
  
  
  
  
  
  Info_Matrix= marginal.reliability(ipar, maxCat = maxCat)$scale.info
  
  Tscore_UL=50-((length(Info_Matrix))-1)/2
  Tscore_LL=((length(Info_Matrix))-1)/2+50
  names(Info_Matrix)<-c(Tscore_UL:Tscore_LL)
  
  SE_Matrix= 1/sqrt(Info_Matrix)
  Reliability_Matrix= 1-(SE_Matrix^2)

  if(all(is.na(Data_to_score))){
    par(mar=c(4,4,1,1))
    plot(names(Info_Matrix), Info_Matrix,  type = "l", xlim=c(Tscore_UL,Tscore_LL), yaxt="n", xlab="T Score", ylab="Information")
    abline(h=10)
    axis(side=1, at=10*(1:9))
    axis(side=2, at=5*(1:5))
  }  
  

  else{
    
    Scored_data<-thetaSE.eap(ipar, Data_to_score)$theta*10+50
    
    layout.show( layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1), heights=c(1,2)))
    par(mar=c(0,4,0,0));
    plot(names(Info_Matrix), Info_Matrix,  type = "l", xlim=c(Tscore_UL,Tscore_LL), xaxt="n", ylab="Information")
    abline(h=10)
    axis(side=2, at=5*(1:5))
    par(mar=c(5,4,0,0));
    h<-hist(Scored_data, xlim=c(Tscore_UL,Tscore_LL) ,ylim=c(110,0), main="",  xlab="T Score", ylab=" ", yaxt="n")
    text(h$mids,h$counts,labels=h$counts, adj=c(0.5, 1.25), cex=0.5)
    axis(side=1, at=10*(1:9))
    }
  
    
}