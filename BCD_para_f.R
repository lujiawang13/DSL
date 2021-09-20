parameter_f<-function(ThetaInv,H,Z,pos_c,lamda,mu,Theta,indnoi,i,p,
                      betamethod){
  ThetaInv11 <- ThetaInv[indnoi,indnoi]; ThetaInv12 = ThetaInv[indnoi,i];
  H11 <- H[indnoi,indnoi]; H12 = H[indnoi,i];
  
  Theta11 <- Theta[indnoi,indnoi];
  
  invTheta11 <- ThetaInv11-ThetaInv12%*%t(ThetaInv12)/ThetaInv[i,i];
  invTheta11H12 <- invTheta11%*%H12;
  
  H11invTheta11 <- H11%*%invTheta11;
  W1 <- invTheta11%*%H11invTheta11;

  beta <- Theta[indnoi,i];   #beta0  
  beta0<-beta
  gam=0;

  iter=0
  repeat{
    iter=iter+1
    betaold <- beta;
    gam.old<-gam
    
    ####### Update gam
   phi <- t(beta)%*%W1%*%beta-2*t(beta)%*%invTheta11H12+H[i,i];
    if(mu==0 | pos_c[i]==0){
      if(lamda==0){
        gam <- phi;
      }else{
        gam <- (-1+sqrt(1+4*lamda*phi))/(2*lamda);
      }
    }else{
      delta=-(t(beta)%*%invTheta11%*%beta-Z[i,i])
      if(lamda==mu){
        gam_1 <- phi;
      }else{
        gam_1 <- (-1+sqrt(1+4*(lamda-mu)*phi))/(2*(lamda-mu));
      }
      gam_2 <- (-1+sqrt(1+4*(lamda+mu)*phi))/(2*(lamda+mu));
      if(gam_1>delta){gam=gam_1}
      if(gam_2<delta){gam=gam_2}
    }
    ####### Update beta          
    V <- W1/as.numeric(gam)+lamda*invTheta11; # A in paper
    d = invTheta11H12/as.numeric(gam);
    
    
    beta<-beta_solution_f(beta,V,d,p,lamda,mu,epsilon_beta=median(abs(beta0)/100),
                          gam,invTheta11,Z,pos_c,indnoi,i,betamethod)
                          
    #if(abs(gam.old-gam)<0.01*abs(gam) & max(abs(betaold-beta))<median(abs(beta0)/100 )| iter>100){
    if(abs(gam.old-gam)<abs(gam)*1e-6 & abs(norm_f(betaold-beta,'2'))<1e-6| iter>100){
      break;
    }

  }
  return(list(gam=gam,beta=beta,invTheta11=invTheta11,iter=iter))
}

beta_solution_f<-function(beta,V,d,p,lamda,mu,epsilon_beta,
                          gam,invTheta11,Z,pos_c,indnoi,i,betamethod){
  if(betamethod=='majorizationminimization'){
    # print(betamethod)
    Z12 = Z[indnoi,i]
    deltah_1=mu*pos_c[i]*(diag(pos_c[-i])%*%sign(beta-Z12))
    deltah_2=mu*pos_c[i]*(invTheta11%*%beta)*as.numeric(sign(t(beta)%*%invTheta11%*%beta+(gam-Z[i,i])))
    deltah=deltah_1+deltah_2
    dnew=d+deltah
    # objective: beta*V*beat-2dnew*beta+2*lambda*||beta||_1
    # solver
    q=length(beta)
    diag2=cbind(diag(q),-diag(q))
    Dmat=t(diag2)%*%V%*%diag2
    dvec=t(dnew)%*%diag2-lamda*rep(1,q*2)
    Amat=diag(2*q)
    bvec=rep(0,2*q)
    # meq=0
    # solve_info=solve.QP(Dmat, dvec, Amat, bvec, meq=meq, factorized=FALSE)
    # alphas=solve_info$solution
    # beta=alphas[1:q]-alphas[(q+1):(2*q)]
    ################
    solve_info=ipop(c=-dvec, H=Dmat, A=Amat, b=bvec, l=bvec, u=bvec+999999, r=bvec+999999, sigf = 7, maxiter = 40, margin = 0.05,
                    bound = 10, verb = 0)
    alphas=attributes(solve_info)$primal
    beta=alphas[1:q]-alphas[(q+1):(2*q)]
  }else{
    if(mu==0 | pos_c[i]==0){
      beta<-beta_solution_fwithoutmu(beta,V,d,p,lamda,epsilon_beta)
    }else{
      Z11 <- Z[indnoi,indnoi]; Z12 = Z[indnoi,i]
      pos_c=pos_c[-i]
      iter=0
      repeat{
        iter=iter+1
        betaold <- beta;
        
        for(j in 1:(p-1)){
          #betaj.old <- beta[j];
          #
          A_1=V-mu*invTheta11
          x_10 = (d[j]+mu*pos_c[j])-(A_1%*%beta)[j]+A_1[j,j]*beta[j];
          x_1 <- max(0,abs(x_10)-lamda)*sign(x_10)/(V[j,j]-mu*invTheta11[j,j]);
          beta_1=beta
          beta_1[j]=x_1
          #
          A_2=V+mu*invTheta11
          x_20 = (d[j]+mu*pos_c[j])-(A_2%*%beta)[j]+A_2[j,j]*beta[j];
          x_2 <- max(0,abs(x_20)-lamda)*sign(x_20)/(V[j,j]+mu*invTheta11[j,j]);
          beta_2=beta
          beta_2[j]=x_2
          #
          A_3=V-mu*invTheta11
          x_30 = (d[j]-mu*pos_c[j])-(A_3%*%beta)[j]+A_3[j,j]*beta[j];
          x_3 <- max(0,abs(x_30)-lamda)*sign(x_30)/(V[j,j]-mu*invTheta11[j,j]);
          beta_3=beta
          beta_3[j]=x_3
          #
          A_4=V+mu*invTheta11
          x_40 = (d[j]-mu*pos_c[j])-(A_4%*%beta)[j]+A_4[j,j]*beta[j];
          x_4 <- max(0,abs(x_40)-lamda)*sign(x_40)/(V[j,j]+mu*invTheta11[j,j]);
          beta_4=beta
          beta_4[j]=x_4
          
          # print("Hi")
          # print(x_1)
          # print(Z12[j])
          # print(t(beta_1)%*%invTheta11%*%beta_1)
          # print(Z[j,j]-gam)
          
          x_candidate<-c()
          if(x_1>=Z12[j] & (t(beta_1)%*%invTheta11%*%beta_1)>=(Z[j,j]-gam)){
            x_candidate<-c(x_candidate,x_1)
          }
          if(x_2>=Z12[j] & (t(beta_2)%*%invTheta11%*%beta_2)<(Z[j,j]-gam)){
            x_candidate<-c(x_candidate,x_2)
          }
          if(x_3<Z12[j] & (t(beta_3)%*%invTheta11%*%beta_3)>=(Z[j,j]-gam)){
            x_candidate<-c(x_candidate,x_3)
          }
          if(x_4<Z12[j] & (t(beta_4)%*%invTheta11%*%beta_4)<(Z[j,j]-gam)){
            x_candidate<-c(x_candidate,x_4)
          }
          
          if(length(x_candidate)==1){
            x=x_candidate
          }else if(length(x_candidate)>1){
            # L=BetaOpt_f(x,beta,V,d,lamda,mu,pos_c,Z12,invTheta11,gam,Z,j)
            betaL=sapply(x_candidate,BetaOpt_f,beta,V,d,lamda,mu,pos_c,Z12,invTheta11,gam,Z,j)
            x=x_candidate[which(betaL==min(betaL))]
          }else{
            print("Error in solving beta")
          }
          
          beta[j] <- x
          
          
          # beta.change = beta[j]-betaj.old;
          # Vbeta = V%*%beta;
          # if(beta.change!=0)
          # {
          #   Vbeta = Vbeta+beta.change*V[,j];
          # }
        }
        
        if(max(abs(betaold-beta))<epsilon_beta | iter>100){
          break;
        }
        
      }# end iteration
    }
  }
  return(beta)
}

beta_solution_fwithoutmu<-function(beta,V,d,p,lamda,epsilon_beta){
  Vbeta = V%*%beta;
  
  iter=0
  repeat{
    iter=iter+1
    
    betaold <- beta;
    
    for(j in 1:(p-1)){
      betaj.old <- beta[j];
      x = d[j]-Vbeta[j]+V[j,j]*beta[j]; # Vbeta[j]-V[j,j]*beta[j] from theta*A*theta
      
      beta[j] <- max(0,abs(x)-lamda)*sign(x)/V[j,j];
      
      beta.change = beta[j]-betaj.old;
      if(beta.change!=0)
      {
        Vbeta = Vbeta+beta.change*V[,j];
      }
    }
    
    if(max(abs(betaold-beta))<epsilon_beta | iter>100){
      break;
    }
    
  }
  return(beta)
}

BetaOpt_f<-function(x,beta,V,d,lamda,mu,pos_c,Z12,invTheta11,gam,Z,j){
  beta[j] <- x
  L=t(beta)%*%V%*%beta-2*t(d)%*%beta+2*lamda*sum(abs(beta))-
    2*mu*sum(abs(diag(pos_c)%*%(beta-Z12)))-mu*abs(t(beta)%*%invTheta11%*%beta+gam-Z[j,j])
  return(L)
}