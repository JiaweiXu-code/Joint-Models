
### load libraries ###
library(survival)
library(plyr)
breaks = c(1.91, 2.43, 3.00, 3.80)

### create container for estimates and powers ###
ests = matrix(NA,13,13)
powers = rep(NA,13)

for (node.idx in 1:13){
  temp.est = matrix(NA,4000,13)
  
  for (ii in (80*(node.idx-1)+1):(80*node.idx)){
    fname = paste("/pine/scr/j/i/jiawei/rev11/11/data/data",ii,".csv",sep = "")
    datas = read.csv(fname, header = F, col.names = c("ID","measure","time","trt","toltime","nodegrp","sim","t","sr","r0","censorship"))
    surv = datas[which(datas$time == 0),]
    ints = survSplit(Surv(t,censorship) ~ ., data = surv, cut = breaks, episode = "interval", 
                     start = "tstart", end = "tstop")
    diffs = mutate(ints, exposure = tstop - tstart, 
                   interval=factor(interval
                                   ,labels=paste("(", c(0,breaks), ",", c(breaks,round(max(datas$t),2)), "]", sep=""))) 
    idx = ii%%80
    if (idx == 0){ idx = 80}
    for (jj in (50*(idx-1)+1):(50*idx)){
      sims = diffs[which(diffs$sim==jj),] 
      est = glm(censorship ~ interval + trt + nodegrp + offset(log(exposure)), data = sims, family = poisson)
      #summary(est)
      b = as.vector(est$coefficients)
      if (length(b) == 6) { 
        bb = c(b[1],b[1]+b[2:4],NA,b[5:6])
        se = summary(est)$coefficients[5,2]
      }else{
        bb = c(b[1],b[1]+b[2:5],b[6:7])
        se = summary(est)$coefficients[6,2]
      }
      c = exp(bb[1:5])
      temp.est[jj,] = t(c(bb,c,se))
    }
  }
  means = colMeans(temp.est,na.rm = T)
  ests[node.idx,] = means
  power = mean(ifelse(pnorm(-temp.est[,6]/temp.est[,13])>0.95,1,0))
  powers[node.idx] = power
}

est.mean = colMeans(ests)

write.table(est.mean,"/pine/scr/j/i/jiawei/rev11/11/programs/PH_est.csv",sep = ",",col.names = F,row.names = F)
write.table(powers,"/pine/scr/j/i/jiawei/rev11/11/programs/PH_power.csv",sep = ",",col.names = F,row.names = F)






