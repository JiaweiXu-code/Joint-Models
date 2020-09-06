
######## misspecification ########

library(sas7bdat)
est_concave = read.sas7bdat("C:/D/paper1/revision/rev11/mis/zest_concave.sas7bdat")
est_leveloff = read.sas7bdat("C:/D/paper1/revision/rev11/mis/zest_leveloff.sas7bdat")
est_type1_concave = read.sas7bdat("C:/D/paper1/revision/rev11/mis/type1/zest_type1_concave.sas7bdat")
est_type1_leveloff = read.sas7bdat("C:/D/paper1/revision/rev11/mis/type1/zest_type1_leveloff.sas7bdat")

##### 0.5/1 ##### concave ####
t = c(0,0.5,1,2)
a = c(-9.113372e-04,5.720166e-04,3.205325e-04)
b = c(0.3764371,0.3724366,-9.067116e-02)
intercept = 0.2699800
trt = -7.135709e-03

# control #
g1 = a[1]*t[2]
g2 = a[2]*(t[3]-t[2])+g1
g3 = a[3]*(t[4]-t[3])+g2
control = c(0,g1,g2,g3)+intercept
# treated #
t1 = (a[1]+b[1])*t[2]
t2 = (a[2]+b[2])*(t[3]-t[2])+t1
t3 = (a[3]+b[3])*(t[4]-t[3])+t2
treated = c(0,t1,t2,t3)+intercept+trt

plot(x = t, y = treated, type = "l",ylim=c(-0.5,1),ylab = "QOL")
lines(t,control,col = "red")

##### 0.5/1 ##### leveloff ####
t = c(0,0.5,1,2)
a = c(-1.196570e-03,1.156407e-03,3.850760e-04)
b = c(0.4568211,0.2585472,1.786127e-02)
intercept = 2.702200e-01
trt = 2.779611e-03

# control #
g1 = a[1]*t[2]
g2 = a[2]*(t[3]-t[2])+g1
g3 = a[3]*(t[4]-t[3])+g2
control = c(0,g1,g2,g3)+intercept
# treated #
t1 = (a[1]+b[1])*t[2]
t2 = (a[2]+b[2])*(t[3]-t[2])+t1
t3 = (a[3]+b[3])*(t[4]-t[3])+t2
treated = c(0,t1,t2,t3)+intercept+trt

plot(x = t, y = treated, type = "l",ylim=c(-0.5,1),ylab = "QOL")
lines(t,control,col = "red")

##### 0.5/1 ##### type1_concave ####
t = c(0,0.5,1,2)
a = c(3.733784e-01,3.737709e-01,-8.986055e-02)
b = c(3.613853e-03,-1.996344e-03,-5.615411e-04)
intercept = 2.637021e-01
trt = -1.289756e-03

# control #
g1 = a[1]*t[2]
g2 = a[2]*(t[3]-t[2])+g1
g3 = a[3]*(t[4]-t[3])+g2
control = c(0,g1,g2,g3)+intercept
# treated #
t1 = (a[1]+b[1])*t[2]
t2 = (a[2]+b[2])*(t[3]-t[2])+t1
t3 = (a[3]+b[3])*(t[4]-t[3])+t2
treated = c(0,t1,t2,t3)+intercept+trt

plot(x = t, y = treated, type = "l",ylim=c(-0.5,1),ylab = "QOL")
lines(t,control,col = "red")

##### 0.5/1 ##### type1_leveloff ####
t = c(0,0.5,1,2)
a = c(4.536036e-01,2.617312e-01,1.840724e-02)
b = c(1.109399e-03,-1.732212e-03,2.959097e-04)
intercept = 2.741821e-01
trt = -5.782955e-04

# control #
g1 = a[1]*t[2]
g2 = a[2]*(t[3]-t[2])+g1
g3 = a[3]*(t[4]-t[3])+g2
control = c(0,g1,g2,g3)+intercept
# treated #
t1 = (a[1]+b[1])*t[2]
t2 = (a[2]+b[2])*(t[3]-t[2])+t1
t3 = (a[3]+b[3])*(t[4]-t[3])+t2
treated = c(0,t1,t2,t3)+intercept+trt

plot(x = t, y = treated, type = "l",ylim=c(-0.5,1),ylab = "QOL")
lines(t,control,col = "red")






##### 0.25/0.5/0.75/1/1.25 ##### true model ####
t = c(0,0.25,0.5,0.75,1,1.25,2)
a = rep(0,6)
b = seq(0.5,0,-0.1)            #leveloff
b = c(0.3,0.4,0.5,0.2,0,-0.1)  #concave
a = seq(0.5,0,-0.1)            #type1-leveloff
a = c(0.3,0.4,0.5,0.2,0,-0.1)  #type1-concave
b = rep(0,6)
intercept = 0.26954762

# control #
g1 = a[1]*t[2]
g2 = a[2]*(t[3]-t[2])+g1
g3 = a[3]*(t[4]-t[3])+g2
g4 = a[4]*(t[5]-t[4])+g3
g5 = a[5]*(t[6]-t[5])+g4
g6 = a[6]*(t[7]-t[6])+g5
control = c(0,g1,g2,g3,g4,g5,g6)+intercept
# treated #
t1 = (a[1]+b[1])*t[2]
t2 = (a[2]+b[2])*(t[3]-t[2])+t1
t3 = (a[3]+b[3])*(t[4]-t[3])+t2
t4 = (a[4]+b[4])*(t[5]-t[4])+t3
t5 = (a[5]+b[5])*(t[6]-t[5])+t4
t6 = (a[6]+b[6])*(t[7]-t[6])+t5
treated = c(0,t1,t2,t3,t4,t5,t6)+intercept

plot(x = t, y = treated, type = "l",ylim=c(-0.5,1),ylab = "QOL")
lines(t,control,col = "red")