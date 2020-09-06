options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1])

##################### general case ######################
JMD = function(seed, v, k = 2.95, nlinear = 4, measures = 9, npieces = 5, eta = 1, maxt = 10, censor = 5, 
               censorp = 0.05, interval = 0.001, pknots = c(0.25,0.75,1.25), knots = c(1.91, 2.43, 3.00, 3.80), 
               sigma = 0.65642581, p = c(0.5,0.5), mSigma = 0.71202084, gamma = c(0.26634693,0,-0.03393003),
               plinear = c(-0.31869067,
                           -0.72492568,
                           -0.13653627,
                           -0.21580434
               ), pinter = c(0.0,
                             0.0,
                             0.0,
                             0.0
               ),llambda = c(-3.61388568,
                             -2.21636742,
                             -2.25146058,
                             -2.49844602,
                             -2.70051089
               ), beta = -0.15,
               alpha = c(-0.2,
                         0.76753173
               ),tpoints = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)){
  
  set.seed(seed)
  ss = floor(k*v)
  r0 = runif(ss,0,eta)
  point = data.frame(seq(0,maxt+interval,interval))
  
  
  dataS = matrix(NA,0,8)
  dataL = matrix(NA,0,7)
  for (i in 1:ss){
    
    t.acc = r0[i] + tpoints
    z = rbinom(1,1,p[1])
    node = rbinom(1,1,p[2])
    theta = rnorm(1,0,mSigma)
    g = {}
    p.knots = c(tpoints[1],pknots,tpoints[length(tpoints)])
    for (t in 1:measures){
      temp = {}
      for (l in 1:nlinear){
        temp.t = min(max(tpoints[t]-p.knots[l],0),(p.knots[l+1]-p.knots[l]))
        temp = c(temp,temp.t)
      }
      g = rbind(g,temp)
    }
    
    ### trajectory function ###
    mu = theta + plinear %*% t(g) + pinter %*% t(g*z) + as.numeric(gamma %*% rbind(1,z,node))
    ### longitudinal data ###
    e = rnorm(measures,0,sigma)
    y = mu + e
    
    ### time-to-event data ###
    g.knots = c(tpoints[1],pknots,maxt)
    hazards = function(t){
      gg = {}
      for (b in 1:nlinear){
        temp.gg = min(max(t-g.knots[b],0),(g.knots[b+1]-g.knots[b]))
        gg = c(gg,temp.gg)
      }
      mu = theta + plinear %*% gg + pinter %*% gg*z + gamma[1:2]%*%rbind(1,z)
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


seed_n = node.idx
v = 100
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


fname = paste("/pine/scr/j/i/jiawei/rev11/11/ssd/data/mean",node.idx,".csv",sep = "")
write.table(means,file = fname, sep = ",",col.names = F,row.names = F)

