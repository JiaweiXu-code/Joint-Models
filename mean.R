means = matrix(NA,nrow = 3,ncol = 0)
for (j in 1:1000){
  
  fname = paste("/pine/scr/j/i/jiawei/rev11/sigma/0.356/ssd/data/mean",j,".csv",sep = "")
  mean = read.csv(file = fname)
  means = cbind(means,mean)
}
mean = rowMeans(means)
write.table(mean,file = "/pine/scr/j/i/jiawei/rev11/sigma/0.356/ssd/means.csv",sep = ",",col.names = F,row.names = F)
