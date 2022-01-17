solve3FINMF=function(Xs, rank=3, lam1=0.1, lam2=0.1, lam3=0.1, V0s, opts=list(init=0, maxIter=40, tol=0.01, ter=2)){
  uptade_S = function(Xs, Us, Vs, Ss, S, lam3){
    add <- function(x) Reduce("+", x)
    nSites = length(Xs)
    
    part1 = add(lapply(1:nSites, function(x) t(Us[[x]]) %*% Xs[[x]] %*% Vs[[x]])) 
    part2 = add(lapply(1:nSites, function(x) t(Us[[x]]) %*% Us[[x]] %*% (Ss[[x]] + S) %*% t(Vs[[x]]) %*% Vs[[x]])) + lam3 * S 
    
    S_new  = S * part1 / part2    
    S_new[is.na(S_new)] = 0
    
    return(S_new)
  }
  
  uptade_Ss = function(Xs, Us, Vs, Ss, S, lam2){
    nSites = length(Xs)
    Ss_new=list()
    for (x in 1:nSites) {
      part1 = t(Us[[x]]) %*% Xs[[x]] %*% Vs[[x]] 
      part2 = t(Us[[x]]) %*% Us[[x]] %*% (Ss[[x]] + S + lam2 * Ss[[x]]) %*% t(Vs[[x]]) %*% Vs[[x]]
      
      s_new  = Ss[[x]] * part1 / part2    
      s_new[is.na(s_new)] = 0
      Ss_new[[x]]=s_new
    }
    
    return(Ss_new)
  }
  
  uptade_Us = function(Xs, Us, Vs, Ss, S, lam2){
    S_sum=lapply(Ss, function(s)return(s+S))
    nSites = length(Xs)
    Us_new=list()
    for (x in 1:nSites) {
      SV=(S_sum[[x]]) %*% t(Vs[[x]])
      part1 = Xs[[x]] %*% t(SV)
      part2 = Us[[x]] %*% (SV %*% t(SV))
      SV2=Ss[[x]] %*% t(Vs[[x]])
      part2 = part2 + lam2 * Us[[x]] %*% (SV2 %*% t(SV2))
      
      u_new  = Us[[x]] * part1 / part2    
      u_new[is.na(u_new)] = 0
      Us_new[[x]]=u_new
    }
    
    return(Us_new)
  }
  
  calcObj = function(Xs, Us, Vs, Ss, S, lam1, lam2, lam3) {
    loss = vector()
    pen1 = vector()
    pen2 = vector()
    
    for (i in 1:length(Xs)) {
      loss[i] = norm(Xs[[i]] - Us[[i]] %*% (S + Ss[[i]]) %*% t(Vs[[i]]), type = "F")^2
      pen1[i] = lam1 * norm(Vs[[i]], type = "F")^2 
      pen2[i] = lam2 * norm(Us[[i]] %*% Ss[[i]] %*% t(Vs[[i]]), type = "F")^2 
      
    }
    pen3 = lam3 * norm(S, type = "F")^2
    
    res=sum(loss) 
    o=res + sum(pen1) + sum(pen2) + pen3
    return(c(o, res))
  }
  
  #######################
  #start the optimization
  #######################
  #get dimensions
  nSites = length(Xs)
  nSubs = sapply(Xs, dim)[1, ]
  nGene = sapply(Xs, dim)[2, ]
  nPath = ncol(V0s[[1]])
  Vs=V0s
  
  #initialize a starting point
  if(opts$init==0){
    Us = lapply(1:nSites, function(x) matrix(data = runif(n=rank*nSubs[x], min = 0, max = 1), nrow = nSubs[x], ncol = rank))
    Ss = lapply(1:nSites, function(x) matrix(data = runif(n=rank*nPath, min = 0, max = 1), nrow = rank, ncol = nPath))
    S = matrix(data = runif(n=nPath*rank, min = 0, max = 1), nrow = rank, ncol = nPath)
  }else if(opts$init==1){
    Ss <- opts$Ss
    S <- opts$S
    Us <- opts$Us
  }    
  
  obj=vector()
  res=vector()
  for (iter in 1:opts$maxIter){
    print(paste0("iteration: ", iter))
    
    # Start of inner loop
    S = uptade_S(Xs=Xs, Us=Us, Vs=Vs, Ss=Ss, S=S, lam3=lam3)
    Ss = uptade_Ss(Xs=Xs, Us=Us, Vs=Vs, Ss=Ss, S=S, lam2=lam2)
    Us = uptade_Us(Xs=Xs, Us=Us, Vs=Vs, Ss=Ss, S=S, lam2=lam2)
    info=calcObj(Xs=Xs, Us=Us, Vs=Vs, Ss=Ss, S=S, lam2=lam2, lam1=lam1, lam3=lam3)
    obj=c(obj, info[1])
    res=c(res, info[2])
    
    
    #check the rule for termination
    if (iter>1){
      convergence = abs(obj[length(obj)] - obj[length(obj)-1])
      if (opts$ter==1 & convergence <= opts$tol){
        break
      } else if(opts$ter==2 & convergence <= opts$tol*obj[length(obj)-1]){
        break
      }  else if(opts$ter==3 & iter==opts$maxIter){
        break
      } 
    }
  }
  return(list(Us=Us, S=S, Ss=Ss, Vs=Vs, iter_to_conv=iter, objList=obj, resList=res))
}



#X=U x S x V
#X: gene or methylation data
#U: sample variation
#S: pathway signature
#V: mapping matrix

#Two simulation dataset
#X1 (200x200) = U1 (200x2) x S (2x20) x V1 (20x200)
#X2 (150x400) = U2 (150x2) x S (2x20) x V2 (20x450)


V0_1=kronecker(diag(20), rep(1,10))
V0_2=kronecker(diag(20), rep(1,20))
Vs0=list(V0_1, V0_2)
Ss=matrix(0, 2,20); Ss[1, 1:3]=1; Ss[2, 4:6]=1; Ss[1, 7:9]=1; Ss[2, 10:12]=1
U1=matrix(0, 200,2); U1[1:100, 1]=1; U1[101:200, 2]=1; 
X1=U1 %*% Ss %*% t(V0_1)
image(t(apply(X1, 2, rev)))
Ss=matrix(0, 2,20); Ss[1, 1:3]=1; Ss[2, 4:6]=1; Ss[1, 13:15]=1; Ss[2, 16:18]=1
U2=matrix(0, 150,2); U2[1:75, 1]=1; U2[76:150, 2]=1; 
X2=U2 %*% Ss %*% t(V0_2)
Xs=list(X1=X1, X2=X2)
image(t(apply(X2, 2, rev)))


#test
solveOpt=list(init=0, maxIter=5000, tol=0.001, ter=2)
model1=solve3FINMF(Xs, rank=2, lam1=0, lam2=0.01, lam3=0, V0s=Vs0, opts=solveOpt)
image(t(apply(model1$S, 2, rev)))
image(t(apply(model1$Ss[[1]], 2, rev)))
image(t(apply(model1$Ss[[2]], 2, rev)))
image(t(apply(model1$Us[[1]], 2, rev)))
image(t(apply(model1$Us[[2]], 2, rev)))


