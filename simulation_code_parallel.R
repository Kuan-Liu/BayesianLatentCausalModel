#rm(list=ls())

### for SciNet one needs these parameters from outside the file:
OUT="/scratch/k/kuan/kuan/paper2"; WD="~/bayeslatent/paper2";

setwd(OUT)

ncores<-40

options(warn=-1, scipen = 999)

if( !file.exists("model_norm_indicators3.txt") ){
  cat( " model {
     #N = nobs
     for (i in 1:N) {

     #outcome model;
     y[i] ~ dnorm(mu[i], tau)
     mu[i] <- t0+t1*z1[i]+t2*z2[i]+t3*(u1[i]==2)+t4*(u1[i]==3)+t5*(u2[i]==2)+t6*(u2[i]==3)

     # treatment visit 2 model;
     z2[i] ~ dbern(p6[i])
     logit(p6[i]) <- e0 +  e1*z1[i] + e2*(u2[i]==2) + e3*(u2[i]==3)

     #indicators at visit 2;
     for (j in 1:nindic) {
     logit(p5[i,j]) <- b0[j]+b1[j]*(u2[i]==2)+b2[j]*(u2[i]==3)
     indic2[i,j] ~ dbern(p5[i,j])
     }


     # latent class visit 2 model;
     u2[i] ~ dcat(prob2[i,1:3])
     prob2[i,2]<- exp(d20+d21*c_a[i]+d22*c_s[i]+d23*(u1[i]==2)+d24*(u1[i]==3)+d25*z1[i])*prob2[i,1]
     prob2[i,3]<- exp(d30+d31*c_a[i]+d32*c_s[i]+d33*(u1[i]==2)+d34*(u1[i]==3)+d35*z1[i])*prob2[i,1]
     prob2[i,1]<-1/(1+exp(d20+d21*c_a[i]+d22*c_s[i]+d23*(u1[i]==2)+d24*(u1[i]==3)+d25*z1[i])+exp(d30+d31*c_a[i]+d32*c_s[i]+d33*(u1[i]==2)+d34*(u1[i]==3)+d35*z1[i]))

     # treatment visit 1 model;
     z1[i] ~ dbern(p3[i])
     logit(p3[i]) <- c0 + c1*(u1[i]==2)+c2*(u1[i]==3)

     #indicators at visit 1;
     for (k in 1:nindic) {
     logit(p2[i,k]) <- b0[k]+b1[k]*(u1[i]==2)+b2[k]*(u1[i]==3)
     indic1[i,k] ~ dbern(p2[i,k])
     }

     # latent class visit 1 model;
     u1[i] ~ dcat(prob1[i, 1:3])
     prob1[i,2]<- exp(a20+a21*c_a[i]+a22*c_s[i])*prob1[i,1]
     prob1[i,3] <- exp(a30+a31*c_a[i]+a32*c_s[i])*prob1[i,1]
     prob1[i,1]<-1/(1+exp(a20+a21*c_a[i]+a22*c_s[i])+exp(a30+a31*c_a[i]+a32*c_s[i]))

     }

     # Priors
     tau ~ dgamma(.01,.01)
     a20 ~ dnorm(0,0.1)
     a21 ~ dnorm(0,0.1)
     a22 ~ dnorm(0,0.1)
     a30 ~ dnorm(0,0.1)
     a31 ~ dnorm(0,0.1)
     a32 ~ dnorm(0,0.1)

    for (h in 1:nindic) {
     b0[h] ~ dnorm(0,0.1)
     b1[h] ~ dnorm(0,0.1)
     b2[h] ~ dnorm(0,0.1)
    }

     c0 ~ dnorm(0,0.1)
     c1 ~ dnorm(0,0.1)
     c2 ~ dnorm(0,0.1)
     d20 ~ dnorm(0,0.1)
     d21 ~ dnorm(0,0.1)
     d22 ~ dnorm(0,0.1)
     d23 ~ dnorm(0,0.1)
     d24 ~ dnorm(0,0.1)
     d25 ~ dnorm(0,0.1)
     d30 ~ dnorm(0,0.1)
     d31 ~ dnorm(0,0.1)
     d32 ~ dnorm(0,0.1)
     d33 ~ dnorm(0,0.1)
     d34 ~ dnorm(0,0.1)
     d35 ~ dnorm(0,0.1)
     e0 ~ dnorm(0,0.1)
     e1 ~ dnorm(0,0.1)
     e2 ~ dnorm(0,0.1)
     e3 ~ dnorm(0,0.1)
     t0 ~ dnorm(0,0.1)
     t1 ~ dnorm(0,0.1)
     t2 ~ dnorm(0,0.1)
     t3 ~ dnorm(0,0.1)
     t4 ~ dnorm(0,0.1)
     t5 ~ dnorm(0,0.1)
     t6 ~ dnorm(0,0.1)
     }",

       file = "model_norm_indicators3.txt")
}


mysim <- function(outfile, from=1, to=100, nobs=500, nindic=10, qindic=1, confounding="h", samplesize=10000) {

  # from = 1;to=1;nobs=250;nindic=10;qindic=1;confounding="h";samplesize = 10000;

  expit<-function(x){1/(1+exp(-x))}

  nrep<-1000

  library(doParallel)
  registerDoParallel(ncores)

  results<-foreach(i=from:to, .combine='rbind',.inorder=F, .noexport=NULL, .verbose=T) %dopar% {

    #.libPaths(c(.libPaths(),"/hpf/tools/centos6/R/3.6.0/lib64/R/library"))

    library(R2jags)
    library(geepack)
    library(gfoRmula)
    library(ltmle)
    set.seed(1989+i)

    results.it <- matrix(NA, 1, 41)

    #let's first simulate confounder, C, imaging this as sex and age 2 variables;
    c_a <- rnorm(nobs,mean=10, sd=3)
    c_s <- rbinom(nobs,1,prob=0.6)

    #first latent variable (3 classes);
    pu1_1 <- 1/(1+exp(0.5-c_a*0.1 + c_s*0.2)+exp(0.5-c_a*0.1 + c_s*0.2))
    pu1_2 <- exp(0.5-c_a*0.1 + c_s*0.2)*pu1_1
    pu1_3 <- exp(0.5-c_a*0.1 + c_s*0.2)*pu1_1
    vProb <- cbind(pu1_1,pu1_2,pu1_3)
    u1 <- t(apply(vProb, 1, function(x){rmultinom(1, size = 1, prob = x)}))
    u1 <- apply(u1, 1, function(x){which.max(x)})
    # 0.44:0.32:0.24

    #simulating latent class indicators;
    # high quality indicators are predicted by the latent variable to have a probability near 1 or 0.
    #such indicator is desirable, stablize estimation by increasing information available to estimate latent class;
    # let's consider 4 levels of indicators, 0.9=expit(2.2), 0.8=expit(1.4), 0.75=expit(1.1)
    indic1<-matrix(NA, nobs, nindic)
    if ( qindic==1) {
      ph1 <- expit(-2.2+4.4*(u1==1))
      ph2 <- expit(-2.2+4.4*(u1==2))
      ph3 <- expit(-2.2+4.4*(u1==3))
      for (j1 in 1:floor((nindic)/3)){
        indic1[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic1[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic1[,j3]<- rbinom(nobs,1,prob=ph3)
      }
    } else if (qindic==2 & nindic==5) {

      ph1 <- expit(-2.2+4.4*(u1==1))
      ph2 <- expit(-0.4+0.8*(u1==2))
      ph3 <- expit(-0.4+0.8*(u1==3))


      indic1[,1]<- rbinom(nobs,1,prob=ph1)

      for (j2 in 2:3){
        indic1[,j2]<- rbinom(nobs,1,prob=ph2)
      }

      for (j3 in 4:5){
        indic1[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    }  else if (qindic==2 & nindic==10) {

      ph11 <- expit(-2.2+4.4*(u1==1))
      ph1 <- expit(-0.4+0.8*(u1==1))
      ph2 <- expit(-0.4+0.8*(u1==2))
      ph3 <- expit(-0.4+0.8*(u1==3))


      indic1[,1]<- rbinom(nobs,1,prob=ph11)

      for (j2 in 2:3){
        indic1[,j2]<- rbinom(nobs,1,prob=ph1)
      }

      for (j3 in 4:6){
        indic1[,j3]<- rbinom(nobs,1,prob=ph2)
      }

      for (j4 in 7:10){
        indic1[,j4]<- rbinom(nobs,1,prob=ph3)
      }

    } else if (qindic==3) {

      ph1 <- expit(-1.1+2.2*(u1==1))
      ph2 <- expit(-1.1+2.2*(u1==2))
      ph3 <- expit(-1.1+2.2*(u1==3))

      for (j1 in 1:floor((nindic)/3)){
        indic1[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic1[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic1[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    } else if (qindic==4) {

      ph1 <- expit(-0.4+0.8*(u1==1))
      ph2 <- expit(-0.4+0.8*(u1==2))
      ph3 <- expit(-0.4+0.8*(u1==3))

      for (j1 in 1:floor((nindic)/3)){
        indic1[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic1[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic1[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    }

    pz1 <- expit(-1+1*(u1==2)+1.5*(u1==3))
    z1 <- rbinom(nobs,1,prob=pz1)

    #second visit;
    pu2_1 <- 1/(1+exp(0.5-c_a*0.1 + c_s*0.2 + 1*(u1==2) + 0.5*(u1==3) - 1*z1 )+exp(0.5-c_a*0.1 + c_s*0.2 + 0.5*(u1==2) + 1*(u1==3) - 1*z1))
    pu2_2 <- exp(0.5-c_a*0.1 + c_s*0.2 + 1*(u1==2) + 0.5*(u1==3) - 1*z1)*pu2_1
    pu2_3 <- exp(0.5-c_a*0.1 + c_s*0.2 + 0.5*(u1==2) + 1*(u1==3) - 1*z1)*pu2_1
    vProb2 <- cbind(pu2_1,pu2_2,pu2_3)
    u2 <- t(apply(vProb2, 1, function(x){rmultinom(1, size = 1, prob = x)}))
    u2 <- apply(u2, 1, function(x){which.max(x)})
    # 0.44:0.31:0.24

    indic2<-matrix(NA, nobs, nindic)

    if ( qindic==1) {

      ph1 <- expit(-2.2+4.4*(u2==1))
      ph2 <- expit(-2.2+4.4*(u2==2))
      ph3 <- expit(-2.2+4.4*(u2==3))

      for (j1 in 1:floor((nindic)/3)){
        indic2[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic2[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic2[,j3]<- rbinom(nobs,1,prob=ph3)
      }


    } else if (qindic==2 & nindic ==5) {

      ph1 <- expit(-2.2+4.4*(u2==1))
      ph2 <- expit(-0.4+0.8*(u2==2))
      ph3 <- expit(-0.4+0.8*(u2==3))


      indic2[,1]<- rbinom(nobs,1,prob=ph1)

      for (j2 in 2:3){
        indic2[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in 4:5){
        indic2[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    } else if (qindic==2 & nindic ==10) {

      ph11 <- expit(-2.2+4.4*(u2==1))
      ph1 <- expit(-0.4+0.8*(u2==1))
      ph2 <- expit(-0.4+0.8*(u2==2))
      ph3 <- expit(-0.4+0.8*(u2==3))


      indic2[,1]<- rbinom(nobs,1,prob=ph11)

      for (j2 in 2:3){
        indic2[,j2]<- rbinom(nobs,1,prob=ph1)
      }
      for (j3 in 4:6){
        indic2[,j3]<- rbinom(nobs,1,prob=ph2)
      }

      for (j4 in 7:10){
        indic2[,j4]<- rbinom(nobs,1,prob=ph3)
      }

    } else if (qindic==3) {

      ph1 <- expit(-1.1+2.2*(u2==1))
      ph2 <- expit(-1.1+2.2*(u2==2))
      ph3 <- expit(-1.1+2.2*(u2==3))

      for (j1 in 1:floor((nindic)/3)){
        indic2[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic2[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic2[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    } else if (qindic==4) {

      ph1 <- expit(-0.4+0.8*(u2==1))
      ph2 <- expit(-0.4+0.8*(u2==2))
      ph3 <- expit(-0.4+0.8*(u2==3))

      for (j1 in 1:floor((nindic)/3)){
        indic2[,j1]<- rbinom(nobs,1,prob=ph1)
      }
      for (j2 in (floor((nindic)/3)+1):(floor((nindic)/3*2))){
        indic2[,j2]<- rbinom(nobs,1,prob=ph2)
      }
      for (j3 in (floor((nindic)/3*2)+1):nindic){
        indic2[,j3]<- rbinom(nobs,1,prob=ph3)
      }

    }


    pz2 <- expit(-1+1*(u2==2)+1.5*(u2==3)-1*z1)
    z2 <- rbinom(nobs,1,prob=pz2)

    if ( confounding=="h") {
      my<-0.1+0.5*z1+1*z2-0.5*(u1==2)-1*(u1==3)-1*(u2==2)-1.5*(u2==3) #high confounding effect
      y<-rnorm(nobs, mean=my, sd = 1)}
    else if ( confounding=="m") {
      my<-0.1+0.5*z1+1*z2-0.2*(u1==2)-0.5*(u1==3)-0.5*(u2==2)-1*(u2==3) #high confounding effect
      y<-rnorm(nobs, mean=my, sd = 1)
    }


    dat<-data.frame(y, c_a,c_s,indic1, indic2, z1,z2, u1,u2)

    mytrue<-function(myz1,myz2){

      c_a_sim <- rnorm(samplesize, mean=10, sd=3)
      c_s_sim <- rbinom(samplesize,1,prob=0.6)

      pu1_1 <- 1/(1+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2)+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2))
      pu1_2 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2)*pu1_1
      pu1_3 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2)*pu1_1

      pu2_1_u11 <- 1/(1+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 - 1*myz1 )+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2- 1*myz1))
      pu2_2_u11 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 - 1*myz1)*pu2_1_u11
      pu2_3_u11 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 - 1*myz1)*pu2_1_u11

      pu2_1_u12 <- 1/(1+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 1 - 1*myz1 )+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 0.5- 1*myz1))
      pu2_2_u12 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 1- 1*myz1)*pu2_1_u12
      pu2_3_u12 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 0.5 - 1*myz1)*pu2_1_u12

      pu2_1_u13 <- 1/(1+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 0.5- 1*myz1 )+exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 1 - 1*myz1))
      pu2_2_u13 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 + 0.5 - 1*myz1)*pu2_1_u13
      pu2_3_u13 <- exp(0.5-c_a_sim*0.1 + c_s_sim*0.2 +1 - 1*myz1)*pu2_1_u13


      #u1=1,u2=1;
      my_11<-(0.1+0.5*myz1+1*myz2)*pu1_1*pu2_1_u11
      #u1=1,u2=2;
      my_12<-(0.1+0.5*myz1+1*myz2-1)*pu1_1*pu2_2_u11
      #u1=1,u2=3;
      my_13<-(0.1+0.5*myz1+1*myz2-1.5)*pu1_1*pu2_3_u11

      #u1=2,u2=1;
      my_21<-(0.1+0.5*myz1+1*myz2-0.5)*pu1_2*pu2_1_u12
      #u1=2,u2=2;
      my_22<-(0.1+0.5*myz1+1*myz2-0.5-1)*pu1_2*pu2_2_u12
      #u1=2,u2=3;
      my_23<-(0.1+0.5*myz1+1*myz2-0.5-1.5)*pu1_2*pu2_3_u12

      #u1=3,u2=1;
      my_31<-(0.1+0.5*myz1+1*myz2-1)*pu1_3*pu2_1_u13
      #u1=3,u2=2;
      my_32<-(0.1+0.5*myz1+1*myz2-1-1)*pu1_3*pu2_2_u13
      #u1=3,u2=3;
      my_33<-(0.1+0.5*myz1+1*myz2-1-1.5)*pu1_3*pu2_3_u13

      my<- rowSums(cbind(my_11,my_12,my_13,my_21,my_22,my_23,my_31,my_32,my_33))

      return(mean(my))

    }

    results.it[1,1]<-mytrue(1,1)-mytrue(0,0)
    results.it[1,2]<-mytrue(0,0)
    results.it[1,3]<-mytrue(1,1)

    #naive approach, apply boostrap for s.e. estimation;
    dat_naive<-data.frame(y,z1,z2)
    fit1<-lm(y~z1+z2, data = dat_naive)

    results.it[1,4]<-predict(fit1,newdata = data.frame(z1=1,z2=1), type="response") - predict(fit1,newdata = data.frame(z1=0,z2=0), type="response")
    results.it[1,5]<-predict(fit1,newdata = data.frame(z1=0,z2=0), type="response")
    results.it[1,6]<-predict(fit1,newdata = data.frame(z1=1,z2=1), type="response")

    results.it[1,7]<-sqrt(vcov(fit1)[2,2]+vcov(fit1)[3,3]+2*(vcov(fit1)[2,3]))
    results.it[1,8]<-sqrt(vcov(fit1)[1,1])
    results.it[1,9]<-sqrt(vcov(fit1)[1,1]+vcov(fit1)[2,2]+vcov(fit1)[3,3]+2*(vcov(fit1)[1,2])+2*(vcov(fit1)[1,3])+2*(vcov(fit1)[2,3]))


    #naive regressions, adjusting covariates;
    fit2<-lm(y~z1+z2+c_a+c_s+indic1+indic2, data = dat)

    results.it[1,10]<-mean(predict(fit2,newdata = data.frame(z1=1,z2=1,c_a,c_s,indic1,indic2), type="response"))-mean(predict(fit2,newdata = data.frame(z1=0,z2=0,c_a,c_s,indic1,indic2), type="response"))
    results.it[1,11]<-mean(predict(fit2,newdata = data.frame(z1=0,z2=0,c_a,c_s,indic1,indic2), type="response"))
    results.it[1,12]<-mean(predict(fit2,newdata = data.frame(z1=1,z2=1,c_a,c_s,indic1,indic2), type="response"))

    delta<-rep(NA,nrep)
    delta00<-rep(NA,nrep)
    delta11<-rep(NA,nrep)

    for (k in 1:nrep){

      set.seed(k)
      bootidx<-sample(1:nobs, nobs, replace=TRUE)

      yb<-y[bootidx]
      z1b<-z1[bootidx]
      z2b<-z2[bootidx]
      c_ab<-c_a[bootidx]
      c_sb<-c_s[bootidx]
      indic1b<-indic1[bootidx]
      indic2b<-indic2[bootidx]

      datb<-  data.frame(yb,z1b,z2b, c_ab,c_sb,indic1b, indic2b)
      fit2b<-lm(yb~z1b+z2b+c_ab+c_sb+indic1b+indic2b, data = datb)
      delta[k]<-mean(predict(fit2b,  newdata = data.frame(z1b=1,z2b=1,c_ab,c_sb,indic1b,indic2b), type="response"))-mean(predict(fit2b,newdata = data.frame(z1b=0,z2b=0,c_ab,c_sb,indic1b,indic2b), type="response"))
      delta00[k]<-mean(predict(fit2b,newdata = data.frame(z1b=0,z2b=0,c_ab,c_sb,indic1b,indic2b), type="response"))
      delta11[k]<-mean(predict(fit2b,newdata = data.frame(z1b=1,z2b=1,c_ab,c_sb,indic1b,indic2b), type="response"))
    }

    results.it[1,13]<-sd(delta)
    results.it[1,14]<-sd(delta00)
    results.it[1,15]<-sd(delta11)

    #competitive methods 1, MSM;
    #calculating stabilized weights;

    # IPT weights:
    tmodel1 <- glm(z1 ~ c_a+c_s+indic1, family=binomial(link=logit), data = dat)
    tmodel2 <- glm(z2 ~ c_a+c_s+indic2 + z1, family=binomial(link=logit), data = dat)

    smodel1 <- glm(z1 ~ 1, family=binomial(link=logit), data = dat)
    smodel2 <- glm(z2 ~ z1, family=binomial(link=logit), data = dat)

    tlp1 <- as.matrix(cbind(1.0, dat[,2:(3+nindic)])) %*% as.matrix(coef(tmodel1))
    tlp2 <- as.matrix(cbind(1.0, cbind(dat[,2:3],dat[,(4+nindic):(3+2*nindic)]), dat$z1)) %*% as.matrix(coef(tmodel2))

    slp1 <- as.matrix(cbind(rep(1.0, nobs))) %*% as.matrix(coef(smodel1))
    slp2 <- as.matrix(cbind(1.0, dat$z1)) %*% as.matrix(coef(smodel2))

    pt <- (exp(dat$z1 * tlp1)/(1+exp(tlp1))) * (exp(dat$z2 * tlp2)/(1+exp(tlp2)))
    sc <- (exp(dat$z1 * slp1)/(1+exp(slp1))) * (exp(dat$z2 * slp2)/(1+exp(slp2)))

    iptw2 <- 1.0/pt
    iptws2 <- sc * iptw2

    #stabilized weights: sandwich variance;
    id<-seq(1, nobs, by=1)
    dat2<-data.frame(dat,id,iptws2)
    wGEE_stable<-geese(y~z1+z2, id = id, data =dat2, weights = iptws2 ,family = gaussian, corstr = "independence")

    results.it[1,16]<-wGEE_stable$beta[2]+wGEE_stable$beta[3]
    results.it[1,17]<-wGEE_stable$beta[1]
    results.it[1,18]<-sum(wGEE_stable$beta)

    #sandwich variance;
    results.it[1,19]<-sqrt(wGEE_stable$vbeta[2,2]+wGEE_stable$vbeta[3,3]+2*wGEE_stable$vbeta[2,3])
    results.it[1,20]<-sqrt(wGEE_stable$vbeta[1,1])
    results.it[1,21]<-sqrt(wGEE_stable$vbeta[1,1]+wGEE_stable$vbeta[2,2]+wGEE_stable$vbeta[3,3]+2*wGEE_stable$vbeta[2,3]+2*wGEE_stable$vbeta[1,2]+2*wGEE_stable$vbeta[2,3])


    #competitive method 2, G-computation 1;
    # use gfoRmula package
    id<-rep(1:nobs, 2)
    t0<-c(rep(0, nobs), rep(1,nobs))
    Y<-c(rep(NA,nobs),y)
    indic<-rbind(indic1,indic2)
    A<-c(z1,z2)
    L1<-c(c_a, c_a)
    L2<-c(c_s, c_s)
    datg<-data.frame(id, t0, Y, L1, L2, indic, A)
    datg <- datg[order(id, t0),]


    if (nindic==5){

      covnames<-c("X1","X2","X3","X4","X5","A")
      outcome_type <- 'continuous_eof'
      covtypes <- c(c(rep('binary',6)))
      histories <- c(lagged)
      histvars <- list(c( 'X1','X2','X3','X4','X5','A'))

      covparams <- list(covmodels = c(X1 ~ lag1_A + + L1 + L2 ,
                                      X2 ~ lag1_A + X1 + L1 + L2,
                                      X3 ~ lag1_A + X2 + X1 + L1 + L2 ,
                                      X4 ~ lag1_A + X3 + X2 + X1 + L1 + L2,
                                      X5 ~ lag1_A + X4 + X3 + X2 + X1 +L1 + L2 ,
                                      A ~ lag1_A + X1 + X2 + X3 + X4 + X5 + L1 + L2))


      ymodel <- Y ~ A + L1 + L2 + X1 + X2 + X3 + X4 + X5 + lag1_A + lag1_X1 + lag1_X2 + lag1_X3 + lag1_X4 + lag1_X5

      intvars <- list('A', 'A')
      interventions <- list(list(c(static, rep(0, 2))),
                            list(c(static, rep(1, 2))))
      int_descript <- c('Never treat', 'Always treat')


      gform_cont_eof <- gformula(obs_data = datg,
                                 id = 'id', time_name = 't0',
                                 covnames = covnames, outcome_name = 'Y',
                                 outcome_type = outcome_type, covtypes = covtypes,
                                 covparams = covparams, ymodel = ymodel,
                                 intvars = intvars, interventions = interventions,
                                 int_descript = int_descript, ref_int = 0,
                                 histories = histories, histvars = histvars,
                                 basecovs = c("L1","L2"), nsimul = nrep, seed = 1234, nsamples = nobs, parallel = TRUE, ncores = 2)

      results.it[1,22]<-gform_cont_eof$result$`Mean difference`[3]
      results.it[1,23]<-gform_cont_eof$result$`g-form mean`[2]
      results.it[1,24]<-gform_cont_eof$result$`g-form mean`[3]

      results.it[1,25]<- gform_cont_eof$result$`MD SE`[3]
      results.it[1,26]<- gform_cont_eof$result$`Mean SE`[2]
      results.it[1,27]<- gform_cont_eof$result$`Mean SE`[3]

    } else if (nindic==10){

      covnames<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","A")
      outcome_type <- 'continuous_eof'
      covtypes <- c(c(rep('binary',11)))
      histories <- c(lagged)
      histvars <- list(c( 'X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','A'))

      covparams <- list(covmodels = c(X1 ~ lag1_A + + L1 + L2 ,
                                      X2 ~ lag1_A + X1 + L1 + L2,
                                      X3 ~ lag1_A + X2 + X1 + L1 + L2 ,
                                      X4 ~ lag1_A + X3 + X2 + X1 + L1 + L2,
                                      X5 ~ lag1_A + X4 + X3 + X2 + X1 +L1 + L2 ,
                                      X6 ~ lag1_A + X5 + X4 + X3 + X2 + X1 + L1 + L2,
                                      X7 ~ lag1_A + X6 + X5 + X4 + X3 + X2 + X1 + L1 + L2 ,
                                      X8 ~ lag1_A + X7 + X6 + X5 + X4 + X3 + X2 + X1 + L1 + L2,
                                      X9 ~ lag1_A + X8 + X7 + X6 + X5 + X4 + X3 + X2 + X1 + L1 + L2 ,
                                      X10 ~ lag1_A + X9 + X8 + X7 + X6 + X5 + X4 + X3 + X2 + X1 + L1 + L2,
                                      A ~ lag1_A + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + L1 + L2))

      ymodel <- Y ~ A + L1 + L2 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +
        + lag1_A + lag1_X1 + lag1_X2 + lag1_X3 + lag1_X4 + lag1_X5 + lag1_X6 + lag1_X7 + lag1_X8 + lag1_X9 + lag1_X10

      intvars <- list('A', 'A')
      interventions <- list(list(c(static, rep(0, 2))),
                            list(c(static, rep(1, 2))))
      int_descript <- c('Never treat', 'Always treat')


      gform_cont_eof <- gformula(obs_data = datg,
                                 id = 'id', time_name = 't0',
                                 covnames = covnames, outcome_name = 'Y',
                                 outcome_type = outcome_type, covtypes = covtypes,
                                 covparams = covparams, ymodel = ymodel,
                                 intvars = intvars, interventions = interventions,
                                 int_descript = int_descript, ref_int = 0,
                                 histories = histories, histvars = histvars,
                                 basecovs = c("L1","L2"), nsimul = nrep, seed = 1234, nsamples = nobs, parallel = TRUE, ncores = 2)

      results.it[1,22]<-gform_cont_eof$result$`Mean difference`[3]
      results.it[1,23]<-gform_cont_eof$result$`g-form mean`[2]
      results.it[1,24]<-gform_cont_eof$result$`g-form mean`[3]

      results.it[1,25]<- gform_cont_eof$result$`MD SE`[3]
      results.it[1,26]<- gform_cont_eof$result$`Mean SE`[2]
      results.it[1,27]<- gform_cont_eof$result$`Mean SE`[3]

    }


    #ltmle;
    if (nindic==5){

      datt<- dat[,c(2:8, 14, 9:13 ,15, 1)]
      colnames(datt)[1:2]<-c("W1","W2")
      colnames(datt)[3:7] <- paste("L", 1:5, ".1", sep = "")
      colnames(datt)[8]<-c("A1")
      colnames(datt)[9:13] <- paste("L", 1:5, ".2", sep = "")
      colnames(datt)[14]<-c("A2")
      colnames(datt)[15]<-"Y"

      tmle_model <- ltmle (datt ,Anodes = c ("A1" , "A2") ,
                           Lnodes = c ("L1.1" , "L2.1","L3.1" , "L4.1","L5.1" ,
                                       "L1.2" , "L2.2","L3.2" , "L4.2","L5.2") , Ynodes = "Y",
                           survivalOutcome=FALSE,
                           gform=c("A1 ~ W1 + W2 + L1.1 + L2.1 + L3.1 + L4.1 + L5.1",
                                   "A2 ~ W1 + W2 + A1 + L1.2 + L2.2 + L3.2 + L4.2 + L5.2"),
                           abar = list(c(1,1),c(0,0)))

      results.it[1,28]<-summary(tmle_model)$effect.measures$ATE$estimate
      results.it[1,29]<-summary(tmle_model)$effect.measures$control$estimate
      results.it[1,30]<-summary(tmle_model)$effect.measures$treatment$estimate

      results.it[1,31]<-summary(tmle_model)$effect.measures$ATE$std.dev
      results.it[1,32]<-summary(tmle_model)$effect.measures$control$std.dev
      results.it[1,33]<-summary(tmle_model)$effect.measures$treatment$std.dev

    } else if (nindic==10){
      datt<- dat[,c(2:13,24, 14:23 ,25,1 )]
      colnames(datt)[1:2]<-c("W1","W2")
      colnames(datt)[3:12] <- paste("L", 1:10,".1",sep = "")
      colnames(datt)[13]<-c("A1")
      colnames(datt)[14:23] <- paste("L", 1:10,".2", sep = "")
      colnames(datt)[24]<-c("A2")
      colnames(datt)[25]<-"Y"

      tmle_model <- ltmle (datt ,Anodes = c ("A1" , "A2") ,
                           Lnodes = c ("L1.1","L2.1","L3.1" , "L4.1","L5.1" ,
                                       "L6.1","L7.1","L8.1","L9.1" , "L10.1",
                                       "L1.2","L2.2","L3.2","L4.2","L5.2" ,
                                       "L6.2","L7.2","L8.2","L9.2","L10.2") , Ynodes = "Y", survivalOutcome=FALSE,
                           gform=c("A1 ~ W1 + W2 + L1.1 + L2.1 + L3.1 + L4.1 + L5.1 + L6.1 + L7.1 + L8.1 + L9.1 + L10.1",
                                   "A2 ~ W1 + W2 + A1 + L1.2 + L2.2 + L3.2 + L4.2 + L5.2 + L6.2 + L7.2 + L8.2 + L9.2 + L10.2"),
                           abar = list(c(1,1),c(0,0)))

      results.it[1,28]<-summary(tmle_model)$effect.measures$ATE$estimate
      results.it[1,29]<-summary(tmle_model)$effect.measures$control$estimate
      results.it[1,30]<-summary(tmle_model)$effect.measures$treatment$estimate

      results.it[1,31]<-summary(tmle_model)$effect.measures$ATE$std.dev
      results.it[1,32]<-summary(tmle_model)$effect.measures$control$std.dev
      results.it[1,33]<-summary(tmle_model)$effect.measures$treatment$std.dev

    }


    #Bayesian specification version 1;
    jags.data<-list(c_a= c_a, c_s=c_s,indic1=indic1, z1= z1,indic2=indic2,z2= z2, y= y, N = nobs, nindic=nindic)

    jags.params<-c("a20", "a21", "a22","a30", "a31", "a32","c0", "c1","c2",
                   "d20", "d21", "d22", "d23","d24","d25","d30", "d31", "d32", "d33","d34","d35",
                   "e0","e1", "e2","e3","t0","t1", "t2","t3","t4", "t5","t6", "tau")

    jags.inits<-function(){list(a20=0.5, a21=-0.1, a22=0.1, a30=0.5, a31=-0.1, a32=0.1,
                                c0=-1,c1=1,c2=1.5,
                                d20=0.5, d21=-0.1, d22=0.1, d23 =1, d24= 0.5, d25=-1,
                                d30=0.5, d31=-0.1, d32=0.1, d33 =0.5, d34= 1, d35=-1,
                                e0=-1, e1 = -1,e2 = 1,e3=1.5,
                                t0=0.1,t1=0.5, t2=1,t3=-0.2,t4=-0.5, t5=-0.5,t6=-1,tau=1)}

    jagsfit<- jags(data = jags.data, inits = jags.inits, jags.params, n.iter = 50000,
                   model.file = "model_norm_indicators3.txt", n.chains = 1, n.burnin = 40000, n.thin = 5)

    jags.mcmc <- as.mcmc(jagsfit)
    out.mcmc <- as.matrix(jags.mcmc)

    pred_mcmc <- rep(NA,dim(out.mcmc)[1])

    mymcmc<-function(myz1, myz2){

      for (i2 in 1:(dim(out.mcmc)[1])){

        pu1_1 <- 1/(1+exp(out.mcmc[i2,1]+c_a*out.mcmc[i2,2] + c_s*out.mcmc[i2,3])+exp(out.mcmc[i2,4]+c_a*out.mcmc[i2,5] + c_s*out.mcmc[i2,6]))
        pu1_2 <- exp(out.mcmc[i2,1]+c_a*out.mcmc[i2,2] + c_s*out.mcmc[i2,3])*pu1_1
        pu1_3 <- exp(out.mcmc[i2,4]+c_a*out.mcmc[i2,5] + c_s*out.mcmc[i2,6])*pu1_1

        pu2_1_u11 <- 1/(1+exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11] + c_s*out.mcmc[i2,12] + out.mcmc[i2,15]*myz1)+exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+ c_s*out.mcmc[i2,18]+out.mcmc[i2,21]*myz1))
        pu2_2_u11 <- exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11]+c_s*out.mcmc[i2,12]+out.mcmc[i2,15]*myz1)*pu2_1_u11
        pu2_3_u11 <- exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+c_s*out.mcmc[i2,18]+out.mcmc[i2,21]*myz1)*pu2_1_u11

        pu2_1_u12 <- 1/(1+exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11] + c_s*out.mcmc[i2,12]+out.mcmc[i2,13]+ out.mcmc[i2,15]*myz1)+exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+ c_s*out.mcmc[i2,18]+out.mcmc[i2,19]+out.mcmc[i2,21]*myz1))
        pu2_2_u12 <- exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11]+c_s*out.mcmc[i2,12]+out.mcmc[i2,13]+out.mcmc[i2,15]*myz1)*pu2_1_u12
        pu2_3_u12 <- exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+c_s*out.mcmc[i2,18]+out.mcmc[i2,19]+out.mcmc[i2,21]*myz1)*pu2_1_u12

        pu2_1_u13 <- 1/(1+exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11] + c_s*out.mcmc[i2,12]+out.mcmc[i2,14]+ out.mcmc[i2,15]*myz1)+exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+ c_s*out.mcmc[i2,18]+out.mcmc[i2,20]+out.mcmc[i2,21]*myz1))
        pu2_2_u13 <- exp(out.mcmc[i2,10]+c_a*out.mcmc[i2,11]+c_s*out.mcmc[i2,12]+out.mcmc[i2,14]+out.mcmc[i2,15]*myz1)*pu2_1_u13
        pu2_3_u13 <- exp(out.mcmc[i2,16]+c_a*out.mcmc[i2,17]+c_s*out.mcmc[i2,18]+out.mcmc[i2,20]+out.mcmc[i2,21]*myz1)*pu2_1_u13

        #u1=1,u2=1;
        my_11<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2)*pu1_1*pu2_1_u11
        #u1=1,u2=2;
        my_12<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,32])*pu1_1*pu2_2_u11
        #u1=1,u2=3;
        my_13<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,33])*pu1_1*pu2_3_u11

        #u1=2,u2=1;
        my_21<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,30])*pu1_2*pu2_1_u12
        #u1=2,u2=2;
        my_22<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,30]+out.mcmc[i2,32])*pu1_2*pu2_2_u12
        #u1=2,u2=3;
        my_23<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,30]+out.mcmc[i2,33])*pu1_2*pu2_3_u12

        #u1=3,u2=1;
        my_31<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,31])*pu1_3*pu2_1_u13
        #u1=3,u2=2;
        my_32<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,31]+out.mcmc[i2,32])*pu1_3*pu2_2_u13
        #u1=3,u2=3;
        my_33<-(out.mcmc[i2,27]+out.mcmc[i2,28]*myz1+out.mcmc[i2,29]*myz2+out.mcmc[i2,31]+out.mcmc[i2,33])*pu1_3*pu2_3_u13

        pred_mcmc[i2]<- mean(rowSums(cbind(my_11,my_12,my_13,my_21,my_22,my_23,my_31,my_32,my_33)))

      }
      return(pred_mcmc)
    }

    u00_mcmc<-mymcmc(myz1=0, myz2=0)
    u11_mcmc<-mymcmc(myz1=1, myz2=1)

    results.it[1,34]<-mean(u11_mcmc-u00_mcmc)
    results.it[1,35]<-mean(u00_mcmc)
    results.it[1,36]<-mean(u11_mcmc)

    results.it[1,37]<-sd(u11_mcmc-u00_mcmc)
    results.it[1,38]<-sd(u00_mcmc)
    results.it[1,39]<-sd(u11_mcmc)

    results.it[1,40:41] <- quantile((u11_mcmc-u00_mcmc), probs=c(0.025,0.975))

    cbind(i,results.it)


    # print(i)
  }

  # save(results,file=paste0(outfile,".RData"))
  write.table(results, paste0(outfile,".txt"), row.names = FALSE,col.names = FALSE)

}


mysim(paste0("paper2may17_125","_nindic10","_H1"), from = 1, to = 1000, nobs = 125, nindic = 10, qindic=1, confounding="h", samplesize=10000)

