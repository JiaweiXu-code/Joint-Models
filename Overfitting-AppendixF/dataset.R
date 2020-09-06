options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1])

##################### general case ######################
JMD = function(seed, v, k = 2.5, q = 1, measures = 9, npieces = 1, knots = c(1.91, 2.43, 3.00, 3.80), 
               sigma = 0.65642581, p = c(0.5,0.5), gamma = c(0.26634693,0,-0.3489892,0.3,-0.03393003), 
               beta = -0.3, alpha = c(-0.2,0.76753173), mSigma = 0.71202084,
               llambda = -2.66, eta = 1, max = 10, censor = 5, 
               censorp = 0.05, interval = 0.001, tpoints = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)){
  
  set.seed(seed)
  ss = floor(k*v)
  r0 = runif(ss,0,eta)
  point = data.frame(seq(0,max+interval,interval))
  
  
  dataS = matrix(NA,0,8)
  dataL = matrix(NA,0,7)
  for (i in 1:ss){
    
    t.acc = r0[i] + tpoints
    z = rbinom(1,1,p[1])
    node = rbinom(1,1,p[2])
    theta = rnorm(1,0,mSigma)
    
    ### trajectory function ###
    mu = theta + gamma %*% rbind(1,z,tpoints,z*tpoints,node)
    ### longitudinal data ###
    e = rnorm(measures,0,sigma)
    y = mu + e
    
    ### time-to-event data ###
    hazards = function(t){
      mu = theta + gamma[1:4] %*% rbind(1,z,t,z*t)
      if (npieces == 1){
        exp(beta*mu + alpha%*%c(z,node)+llambda)
      }else{
        new.knots = c(0,knots)
        indicator = sum(ifelse(t >= new.knots,1,0))
        llambda0 = llambda[indicator]
        exp(beta*mu + alpha%*%c(z,node)+llambda0)
      }
    }
    
    ### discretization ###
    h = as.numeric(apply(point,1,hazards))
    meanh = (h[-length(h)]+h[-1])/2
    H = cumsum(interval*meanh)
    Surv = exp(-H)
    Surv[length(Surv)] = 0
    Surv = c(1,Surv)
    
    ### simulate event time ###    
    u = runif(1,0,1)
    Surv.temp = Surv[-length(Surv)]
    loc = sum(ifelse(u <= Surv.temp,1,0))
    frac = (Surv[loc]-u)/(Surv[loc]-Surv[loc+1])
    t = point[loc,] + frac*(point[loc+1,]-point[loc,])
    
    
    c = runif(1,0,censor)
    cp = rbinom(1,1,censorp)
    if (cp == 0){
      c = 10^10
    }
    if (t>c){
      s = c
      ind = 0
    }else{
      s = t
      ind = 1
    }
    sr = s + r0[i]
    
    ### save dataset ###   
    dataS = rbind(dataS,c(i,s,sr,r0[i],ind,z,node,seed))
    dataL = rbind(dataL,cbind(i,as.numeric(y),tpoints,z,t.acc,node,seed))
    
  }
  nameS = c("ID","s","sr","r0","censorship","trt","nodegrp","sim")
  nameL = c("ID","measure","time","trt","toltime","nodegrp","sim")
  colnames(dataS) = nameS
  colnames(dataL) = nameL
  # JMdata = list("surv" = dataS,"long" = dataL)
  datas <<- dataS
  datal <<- dataL
  
}


com = 50
sim_n = 4000
gap = sim_n/com
event = seq(100,400,25)
seeds = seq(1,sim_n,1)
v.index = ceiling(node.idx/gap)
v = event[v.index]
temp = node.idx %% gap
if (temp == 0){
  temp = sim_n/com
}

all_est <- matrix(NA, nrow=0, ncol=11)
all_means = matrix(NA,0,ncol = 4)
for (ii in ((temp-1)*com+1):(com*temp)){
  
  seed_n = seeds[ii]
  JMD(seed_n,v)
  
  ###### duration and sample size ######
  datas.1 = datas[order(datas[,3]),]
  datas.2 = datas.1[which(datas.1[,5]==1),]
  
  if (nrow(datas.2)<v){
    cutoff = datas.2[nrow(datas.2),3]
  }else{
    cutoff = datas.2[v,3]
  }
  
  datas.3 = datas[order(datas[,4]),]
  min = datas.3[1,4]
  
  datas.new = datas[which(datas[,4]<cutoff),]
  for (j in 1:nrow(datas.new)){
    if (datas.new[j,3] > cutoff){
      datas.new[j,2] = cutoff - datas.new[j,4]
      datas.new[j,5] = 0
      datas.new[j,3] = cutoff
    }
  }
  mean_dur = cutoff-min
  sample = nrow(datas.new)
  event = sum(datas.new[,5])
  means = c(seed_n,mean_dur,sample,event)
  
  data = merge(datal,datas.new[,1:5],by = "ID")
  measures = 9
  
  newD = matrix(NA,0,ncol = 11)
  for (i in 1:(sample*measures)){
    if (data[i,3] <= data[i,8]){
      newD = rbind(newD,data[i,])
    }
  }
  est <- newD
  
  all_est = rbind(all_est,est)
  all_means = rbind(all_means,means)
  
}

fname1 = paste("/pine/scr/j/i/jiawei/rev11/overfitting/data/data",node.idx,".csv",sep = "")
write.table(all_est,file = fname1, sep = ",",col.names = F,row.names = F)
