
###### Appendix D ######
t = c(0,0.25,0.75,1.25)
a = c(-0.31869067,-0.72492568,-0.13653627,-0.21580434)
intercept = 0.26634693

s1 = -3.61388568
s2 = -2.21636742
s3 = -2.25146058
s4 = -2.49844602
s5 = -2.70051089
knots = c(1.91, 2.43, 3.00, 3.80, 5)

trt = 0
beta = -0.3

# direct effect #
b = rep(0.0,4)
alpha_trt = -0.2

# indirect effect #
b = rep(0.2,4)
alpha_trt = 0

# both effects #
b = rep(0.2,4)
alpha_trt = -0.2

###### compute the survival curves #####

## control survival ##
treatment = 0
A1 = beta*intercept + (beta*trt + alpha_trt)*treatment
A2 = A1 + beta*(a[1] + b[1]*treatment)*t[1]
A3 = A2 + beta*(a[2] + b[2]*treatment)*(t[2]-t[1])
A4 = A3 + beta*(a[3] + b[3]*treatment)*(t[3]-t[2])
B1 = beta*(a[1] + b[1]*treatment)
B2 = beta*(a[2] + b[2]*treatment)
B3 = beta*(a[3] + b[3]*treatment)
B4 = beta*(a[4] + b[4]*treatment)

h1 = exp(s1+A1)*(exp(B1*t[2])-1)/B1
h2 = exp(s1+A2)*(exp(B2*(t[3]-t[2]))-1)/B2 + h1
h3 = exp(s1+A3)*(exp(B3*(t[4]-t[3]))-1)/B3 + h2
h4 = exp(s1+A4)*(exp(B4*(knots[1]-t[4]))-1)/B4 + h3
h5 = exp(s2+A4)*(exp(B4*(knots[2]-t[4]))-exp(B4*(knots[1]-t[4])))/B4 + h4
h6 = exp(s3+A4)*(exp(B4*(knots[3]-t[4]))-exp(B4*(knots[2]-t[4])))/B4 + h5
h7 = exp(s4+A4)*(exp(B4*(knots[4]-t[4]))-exp(B4*(knots[3]-t[4])))/B4 + h6
h8 = exp(s5+A4)*(exp(B4*(knots[5]-t[4]))-exp(B4*(knots[4]-t[4])))/B4 + h7
surv_control = exp(-c(0, h1, h2, h3, h4, h5, h6, h7, h8))


## treated survival ##
treatment = 1
A1 = beta*intercept + (beta*trt + alpha_trt)*treatment
A2 = A1 + beta*(a[1] + b[1]*treatment)*t[1]
A3 = A2 + beta*(a[2] + b[2]*treatment)*(t[2]-t[1])
A4 = A3 + beta*(a[3] + b[3]*treatment)*(t[3]-t[2])
B1 = beta*(a[1] + b[1]*treatment)
B2 = beta*(a[2] + b[2]*treatment)
B3 = beta*(a[3] + b[3]*treatment)
B4 = beta*(a[4] + b[4]*treatment)

m1 = exp(s1+A1)*(exp(B1*t[2])-1)/B1
m2 = exp(s1+A2)*(exp(B2*(t[3]-t[2]))-1)/B2 + m1
m3 = exp(s1+A3)*(exp(B3*(t[4]-t[3]))-1)/B3 + m2
m4 = exp(s1+A4)*(exp(B4*(knots[1]-t[4]))-1)/B4 + m3
m5 = exp(s2+A4)*(exp(B4*(knots[2]-t[4]))-exp(B4*(knots[1]-t[4])))/B4 + m4
m6 = exp(s3+A4)*(exp(B4*(knots[3]-t[4]))-exp(B4*(knots[2]-t[4])))/B4 + m5
m7 = exp(s4+A4)*(exp(B4*(knots[4]-t[4]))-exp(B4*(knots[3]-t[4])))/B4 + m6
m8 = exp(s5+A4)*(exp(B4*(knots[5]-t[4]))-exp(B4*(knots[4]-t[4])))/B4 + m7
surv_treated = exp(-c(0, m1, m2, m3, m4, m5, m6, m7, m8))

plot(x = c(t,knots), y = surv_treated, type = "l", ylim=c(0.7,1),ylab = "Survival Probability",xlab = "Time(years)")
lines(c(t,knots),surv_control,col = "red")

### survival probabilities ###
direct_trt = c(1, 0.9948599, 0.9842142, 0.9735797, 0.9592752, 0.9133950, 0.8656110, 0.8139107, 0.7509320)
direct_con = c(1, 0.9937255, 0.9807531, 0.9678252, 0.9504853, 0.8952583, 0.8383895, 0.7776393, 0.7047876)

indirect_trt = c(1, 0.9937723, 0.9809942, 0.9684431, 0.9521728, 0.9019832, 0.8516172, 0.7994604, 0.7396027)
indirect_con = c(1, 0.9937255, 0.9807531, 0.9678252, 0.9504853, 0.8952583, 0.8383895, 0.7776393, 0.7047876)

both_trt = c(1, 0.9948983, 0.9844124, 0.9740885, 0.9606694, 0.9190086, 0.8767766, 0.8325626, 0.7811690)
both_con = c(1, 0.9937255, 0.9807531, 0.9678252, 0.9504853, 0.8952583, 0.8383895, 0.7776393, 0.7047876)

### save survival probabilities ###
surv_p = rbind(cbind(c(t,knots),direct_con,direct_trt,1),cbind(c(t,knots),both_con,both_trt,2),cbind(c(t,knots),indirect_con,indirect_trt,3))
write.table(surv_p, file = "C:/D/paper1/revision/rev11/surv_curves/survp.csv", sep = ",",row.names = F,col.names = F)

