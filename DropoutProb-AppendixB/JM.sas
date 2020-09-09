libname lib "/pine/scr/j/i/jiawei/rev11/censor/10/means";

data all&sysparm.;
infile "/pine/scr/j/i/jiawei/rev11/censor/10/data/data&sysparm..csv" delimiter = ',';
input ID measure time trt toltime nodegrp sim t sr r0 censorship;
log_t = log(t);
run;

proc sort data=all&sysparm.;
by sim id;
run;

data dataAll&sysparm.;
set all&sysparm.;
by sim id;
  t1 = min(time,0.25);
  t2 = min(max(time - 0.25, 0),(0.75-0.25));
  t3 = min(max(time - 0.75, 0),(1.25-0.75));
  t4 = max(time - 1.25, 0);
if first.id then surv_like = 1;
else surv_like = 0;
run;


ods output SolutionF = longRegParms&sysparm. CovParms = longCovParms&sysparm.;
proc mixed data = dataAll&sysparm. method=reml plots=(none);
by sim;
 model measure = t1 t2 t3 t4 trt t1*trt t2*trt t3*trt t4*trt nodegrp / solution;
 random intercept / subject=id type=un;
run;


data survall&sysparm.;
set all&sysparm.;
by sim id;
  t1 = min(t,1.91);
  t2 = min(max(t - 1.91, 0),(2.43-1.91));
  t3 = min(max(t - 2.43, 0),(3.00-2.43));
  t4 = min(max(t - 3.00, 0),(3.80-3.0));
  t5 = max(t - 3.80, 0);
if first.id then output;
run;

data long&sysparm.;
set survall&sysparm.;
by sim id;
retain int;
if first.id then do;
int = t1;
output;
end;
int = t2;
output;
int = t3;
output;
int = t4;
output;
int = t5;
output;
drop t1-t5;
run;

data piece&sysparm.;
set long&sysparm.;
if int = 0 then delete;
log_t = log(int);
by sim id;
if first.id then c = 0;
c+1;
output;
run;

data cens&sysparm.;
set piece&sysparm.;
by sim id;
new_cen = 0;
if last.id then new_cen = censorship;
run;

data sss&sysparm.;
set cens&sysparm.;
s1 = 0;
s2 = 0;
s3 = 0;
s4 = 0;
s5 = 0;
if c = 1 then s1 = 1;
if c = 2 then s2 = 1;
if c = 3 then s3 = 1;
if c = 4 then s4 = 1;
if c = 5 then s5 = 1;
run;
ods output ParameterEstimates = survRegParms&sysparm.;
proc genmod data = sss&sysparm.;
by sim;
model new_cen(event='1') = s1 s2 s3 s4 s5 trt nodegrp / offset = log_t dist=poisson link=log noint;
run;


data survRegParms&sysparm.;
 set survRegParms&sysparm.;
 where upcase(parameter) ^= 'SCALE';
run;

   data init&sysparm.;
    length parameter varName varPref expression $100.;
    set longRegParms&sysparm.(in=a) survRegParms&sysparm.(in=b) longCovParms&sysparm.(in=c);
    by sim;

	** set parameter name based on different source datasets;
	if Effect > '' then parameter = Effect;
	
	     if c and upcase(CovParm) = 'UN(1,1)' then parameter = 'intercept';
	else if c                                 then parameter = lowcase(covParm);
	
	if c and estimate = 0 then Estimate = 1e-3;
	else if c then estimate = estimate;
	
	     if a then varPref = 'gamma';
	else if b then varPref = 'alpha';
	else if c then varPref = 'sig';

	varName = parameter;
	varName = lowcase(tranwrd(varName,'*','_'));

	if varPref > '' then variable = strip(varPref)||"_"||strip(varName);
	else variable = strip(parameter);

    if find(variable,'gamma','i') or find(variable,'alpha','i') then do;
	 order = 80;
     x = substr(variable,index(variable,"_")+1);
	 y = count(x,'_')+1;
	 do j = 1 to y;
      expression = catx('*',expression,scan(x,j,"_")); 
	 end;
	 if expression = 'intercept' then expression = '';
	end;

	if c then order = 1000;
	output;

	if c then do;
	 order = 90;
	 if find(variable,'intercept','i')   then do;  variable = 'random_intercept'; expression = 'theta0';         end;
     else                                     do;  variable = '';                 expression = '';          end;
	 estimate = .;
     if variable > '' then output;
	end;

	if first.sim then do;
	 order = 99;
     variable = 'beta';
	 estimate = -0.3;
     expression = ' ';
	 output;
	end;

	keep sim order variable estimate expression  ;
   run; proc sort; by sim order variable; run;

   data init_parms&sysparm.;
    set init&sysparm.;
	where find(variable,'random','i')=0;
    if find(variable,'alpha_intercept','i') then delete;
	keep variable estimate sim;
	rename variable=parameter;

   run;
   



ods output ParameterEstimates = lib.jmParmEst&sysparm. AdditionalEstimates = lib.jmEst&sysparm.;
proc nlmixed  data = dataAll&sysparm. /*GCONV=1e-10*/ NOAD QPOINTS=10 MAXFU=10000 MAXIT=1000 REST=100;

by sim; 

parms /data=init_parms&sysparm. BYDATA;

bounds sig_intercept>1e-3, sig_residual>1e-3; 
v11 = sig_intercept**2;
random random_intercept ~ normal(0,v11) subject=id;

  ** longitudinal component;
    muLong = random_intercept + gamma_intercept + gamma_t1*t1 + gamma_t2*t2 + gamma_t3*t3 + gamma_t4*t4 + gamma_trt*trt + (gamma_t1_trt*t1+gamma_t2_trt*t2+gamma_t3_trt*t3+gamma_t4_trt*t4)*trt + gamma_nodegrp*nodegrp;
    loglike = -0.5*log(sig_residual**2) - 0.5/sig_residual**2*(measure-muLong)**2;

	** survival component;
   if surv_like = 1 then do;

    L1 = 0.25;
	L2 = 0.75;
    L3 = 1.25;

	A1 = beta*(random_intercept+gamma_intercept) + (beta*gamma_trt + alpha_trt)*trt + alpha_nodegrp*nodegrp;
	A2 = A1 + beta*(gamma_t1 + gamma_t1_trt*trt)*L1;
	A3 = A2 + beta*(gamma_t2 + gamma_t2_trt*trt)*(L2-L1);
    A4 = A3 + beta*(gamma_t3 + gamma_t3_trt*trt)*(L3-L2);

    B1 = beta*(gamma_t1 + gamma_t1_trt*trt);
    B2 = beta*(gamma_t2 + gamma_t2_trt*trt);
  	B3 = beta*(gamma_t3 + gamma_t3_trt*trt);
    B4 = beta*(gamma_t4 + gamma_t4_trt*trt);

	alphai = alpha_trt*trt + alpha_nodegrp*nodegrp;

        k1 = 1.91;
        k2 = 2.43;
        k3 = 3.00;
        k4 = 3.80;


	** when beta_joint approximation zero the hazard is constant as a function of time and hence the likelihood contribution is different;	  

    if beta > 1e-3 or beta < -1e-3 then do;
	if t <= k1 then do;
        if      t <= L1 then logLike = logLike + censorship*(alpha_s1+A1+B1*t) - exp(alpha_s1+A1)*(exp(B1*t)-1)/B1;
        else if t <= L2 then logLike = logLike + censorship*(alpha_s1+A2+B2*(t-L1)) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(t-L1))-1)/B2;
	    else if t <= L3 then logLike = logLike + censorship*(alpha_s1+A3+B3*(t-L2)) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(t-L2))-1)/B3;
	    else                 logLike = logLike + censorship*(alpha_s1+A4+B4*(t-L3)) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_s1+A4)*(exp(B4*(t-L3))-1)/B4;
      end;
          
	else if t <= k2 then logLike = logLike + censorship*(alpha_s2+A4+B4*(t-L3)) - exp(alpha_s2+A4)/B4*(exp(B4*(t-L3))-exp(B4*(k1-L3))) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_s1+A4)*(exp(B4*(k1-L3))-1)/B4;
    else if t <= k3 then logLike = logLike + censorship*(alpha_s3+A4+B4*(t-L3)) - exp(alpha_s3+A4)/B4*(exp(B4*(t-L3))-exp(B4*(k2-L3))) - exp(alpha_s2+A4)/B4*(exp(B4*(k2-L3))-exp(B4*(k1-L3))) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_s1+A4)*(exp(B4*(k1-L3))-1)/B4;
	else if t <= k4 then logLike = logLike + censorship*(alpha_s4+A4+B4*(t-L3)) - exp(alpha_s4+A4)/B4*(exp(B4*(t-L3))-exp(B4*(k3-L3))) - exp(alpha_s3+A4)/B4*(exp(B4*(k3-L3))-exp(B4*(k2-L3))) - exp(alpha_s2+A4)/B4*(exp(B4*(k2-L3))-exp(B4*(k1-L3))) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_s1+A4)*(exp(B4*(k1-L3))-1)/B4;
    else                 logLike = logLike + censorship*(alpha_s5+A4+B4*(t-L3)) - exp(alpha_s5+A4)/B4*(exp(B4*(t-L3))-exp(B4*(k4-L3))) - exp(alpha_s4+A4)/B4*(exp(B4*(k4-L3))-exp(B4*(k3-L3))) - exp(alpha_s3+A4)/B4*(exp(B4*(k3-L3))-exp(B4*(k2-L3))) - exp(alpha_s2+A4)/B4*(exp(B4*(k2-L3))-exp(B4*(k1-L3))) - exp(alpha_s1+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_s1+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_s1+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_s1+A4)*(exp(B4*(k1-L3))-1)/B4;
    end;
		   
		   

    if 1e-3 > beta > -1e-3 then do;
	    if      t <= k1 then logLike = logLike + censorship*(alpha_s1+alphai) - exp(alpha_s1+alphai)*t;
 	    else if t <= k2 then logLike = logLike + censorship*(alpha_s2+alphai) - exp(alpha_s2+alphai)*(t-k1) - exp(alpha_s1+alphai)*k1;
        else if t <= k3 then logLike = logLike + censorship*(alpha_s3+alphai) - exp(alpha_s3+alphai)*(t-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
	    else if t <= k4 then logLike = logLike + censorship*(alpha_s4+alphai) - exp(alpha_s4+alphai)*(t-k3) - exp(alpha_s3+alphai)*(k3-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
	    else                 logLike = logLike + censorship*(alpha_s5+alphai) - exp(alpha_s5+alphai)*(t-k4) - exp(alpha_s4+alphai)*(k4-k3) - exp(alpha_s3+alphai)*(k3-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
	   end;

    end;
	
	** average hazard ratio;
	t0 = 5;
	del1 = max( 0, min(-(beta*gamma_trt + alpha_trt)/(beta*gamma_t1_trt), L1) );
	del2 = min( L1, max(-(beta*gamma_trt + alpha_trt)/(beta*gamma_t1_trt), 0) );
	del3 = max( L1, min( (beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t2_trt), L2) );
	del4 = min( L2, max( (beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t2_trt), L1) );
	del5 = max( L2, min( (beta*L2*(gamma_t3_trt-gamma_t2_trt)+beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t3_trt), L3) );
	del6 = min( L3, max( (beta*L2*(gamma_t3_trt-gamma_t2_trt)+beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t3_trt), L2) );
	del7 = max( L3, min( (beta*L3*(gamma_t4_trt-gamma_t3_trt)+beta*L2*(gamma_t3_trt-gamma_t2_trt)+beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t4_trt), t0) );
	del8 = min( t0, max( (beta*L3*(gamma_t4_trt-gamma_t3_trt)+beta*L2*(gamma_t3_trt-gamma_t2_trt)+beta*L1*(gamma_t2_trt-gamma_t1_trt)-(beta*gamma_trt + alpha_trt))/(beta*gamma_t4_trt), L3) );


	    phi1 = sign(beta*gamma_t1_trt)*((beta*gamma_trt+alpha_trt)**3+(beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3-(beta*gamma_t1_trt*del1+(beta*gamma_trt+alpha_trt))**3-(beta*gamma_t1_trt*del2+(beta*gamma_trt+alpha_trt))**3)/(3*beta*gamma_t1_trt); 
        phi2 = sign(beta*gamma_t2_trt)*( (beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 + (beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t2_trt*(del3-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t2_trt*(del4-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 )/(3*beta*gamma_t2_trt);
        phi3 = sign(beta*gamma_t3_trt)*( (beta*gamma_t1_trt*L1+beta*gamma_t2_trt*(L2-L1)+(beta*gamma_trt+alpha_trt))**3 + (beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t3_trt*(del5-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t3_trt*(del6-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 )/(3*beta*gamma_t3_trt);
        phi4 = sign(beta*gamma_t4_trt)*( (beta*gamma_t1_trt*L1+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t3_trt*(L3-L2)+(beta*gamma_trt+alpha_trt))**3 + (beta*gamma_t4_trt*(t0-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t4_trt*(del7-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 - (beta*gamma_t4_trt*(del8-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))**3 )/(3*beta*gamma_t4_trt);

        c1 = sign(beta*gamma_t1_trt)*( beta*gamma_t1_trt*L1**2/2 + (beta*gamma_trt+alpha_trt)*(L1-del1-del2) - beta*gamma_t1_trt*del1**2/2 - beta*gamma_t1_trt*del2**2/2 );	
        c2 = sign(beta*gamma_t2_trt)*( beta*gamma_t2_trt*(L2-L1)**2/2 + (beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))*(L2+L1-del3-del4) - beta*gamma_t2_trt*(del3-L1)**2/2 - beta*gamma_t2_trt*(del4-L1)**2/2 );
        c3 = sign(beta*gamma_t3_trt)*( beta*gamma_t3_trt*(L3-L2)**2/2 - beta*gamma_t3_trt*(del5-L2)**2/2 - beta*gamma_t3_trt*(del6-L2)**2/2 + (beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))*(L2+L3-del5-del6) );
        c4 = sign(beta*gamma_t4_trt)*( beta*gamma_t4_trt*(t0-L3)**2/2 - beta*gamma_t4_trt*(del7-L3)**2/2 - beta*gamma_t4_trt*(del8-L3)**2/2 + (beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt+alpha_trt))*(t0+L3-del7-del8) );

        extra = t0*(beta*gamma_trt+alpha_trt) + beta*( gamma_t1_trt*L1**2/2+gamma_t2_trt*(L2-L1)**2/2+gamma_t3_trt*(L3-L2)**2/2+gamma_t4_trt*(t0-L3)**2/2 + gamma_t3_trt*(L3-L2)*(t0-L3)+gamma_t2_trt*(L2-L1)*(t0-L2)+gamma_t1_trt*L1*(t0-L1) );
        eta = 1e-3;

        estimate "Overall Effect eta" (phi1+phi2+phi3+phi4+extra*eta)/(c1+c2+c3+c4+eta*t0);
   
    model t ~  general(loglike);
run;
quit;


