BCD_f<-function(H,Theta0,Z,pos_c,lamda,mu,epsilon,
                betamethod){
  p = nrow(H);
  Theta=Theta0;
  ThetaInv = solve(Theta0); 

   
  indnoiall = mat.or.vec(p-1,p);
  for( i in 1:p){
    indnoi = seq(1,p);
    indnoi = indnoi[-i];
    indnoiall[,i] = indnoi; 
  }

  ob = -determinant(ThetaInv,logarithm = TRUE)$modulus+sum(ThetaInv*H)+
    lamda*sum(abs(Theta))-mu*sum(abs(diag(pos_c)%*%(Theta-Z)%*%diag(pos_c)));

  obj <- ob;

  iter=0
  repeat{
    iter=iter+1
    ob.old = ob;
    Theta.old = Theta;
    for(i in 1:p){
      indnoi <- indnoiall[,i];

      ##### gam and beta otput
      parameter_result<-parameter_f(ThetaInv,H,Z,pos_c,lamda,mu,Theta,indnoi,i,p,
                                    betamethod)
      gam=parameter_result$gam
      beta=parameter_result$beta
      invTheta11=parameter_result$invTheta11
      inner_iter=parameter_result$iter

      ####
      Theta[indnoi,i] <- beta;
      Theta[i,indnoi] <- beta;
      Theta[i,i] <- gam+t(beta)%*%invTheta11%*%beta;

      
      invTheta11beta <-  invTheta11%*%beta;
      ThetaInv[indnoi,indnoi] <- invTheta11+invTheta11beta%*%t(invTheta11beta)/as.numeric(gam);
     
      ThetaInv12 <- -invTheta11beta/as.numeric(gam);

      ThetaInv[indnoi,i]=ThetaInv12;
      ThetaInv[i,indnoi]=ThetaInv12;
      ThetaInv[i,i] = 1/gam
    } # Theta
    
    ob <- -determinant(ThetaInv,logarithm = TRUE)$modulus+sum(ThetaInv*H)+
      lamda*sum(abs(Theta))-mu*sum(abs(diag(pos_c)%*%(Theta-Z)%*%diag(pos_c)));
    obj <- c(obj,ob);
    
    if(abs(ob-ob.old)/abs(ob)<epsilon | abs(norm_f(Theta.old-Theta,'2'))<1e-6 | iter>500){
      break;
    }
  }
  
  return(list(Theta = Theta, ThetaInv = ThetaInv, obj = obj, niter = iter))

}


