libname lib "/pine/scr/j/i/jiawei/rev11/overfitting/means_BL";

data all&sysparm.;
infile "/pine/scr/j/i/jiawei/rev11/overfitting/data/data&sysparm..csv" delimiter = ',';
input ID measure time trt toltime nodegrp sim t sr r0 censorship;
log_t = log(t);
run;

proc sort data=all&sysparm.;
by sim id;
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


ods output SolutionF = longRegParms&sysparm. CovParms = longCovParms&sysparm.;
proc mixed data = dataAll&sysparm. method=reml plots=(none);
by sim;
 model measure = time trt time*trt nodegrp / solution;
 random intercept / subject=id type=un;
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
	 if find(variable,'intercept','i')   then do;  variable = 'random_intercept'; expression = 'theta0';    end;
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
    muLong = random_intercept + gamma_intercept + gamma_time*time + gamma_trt*trt + gamma_time_trt*time*trt + gamma_nodegrp*nodegrp;
    loglike = -0.5*log(sig_residual**2) - 0.5/sig_residual**2*(measure-muLong)**2;

   if surv_like = 1 then do;

	k1 = 1.91;
        k2 = 2.43;
        k3 = 3.00;
        k4 = 3.80;

	A = beta*(random_intercept + gamma_intercept) + (beta*gamma_trt+alpha_trt)*trt + alpha_nodegrp*nodegrp;
        B = beta*(gamma_time + gamma_time_trt*trt);

	alphai = alpha_trt*trt + alpha_nodegrp*nodegrp;

    if beta > 1e-3 or beta < -1e-3 then do;
        if t <= k1      then logLike = logLike + censorship*(alpha_s1+A+B*t) - exp(alpha_s1+A)/B*(exp(B*t)-1);
        else if t <= k2 then logLike = logLike + censorship*(alpha_s2+A+B*t) - exp(alpha_s2+A)/B*(exp(B*t)-exp(B*k1)) - exp(alpha_s1+A)/B*(exp(B*k1)-1);
        else if t <= k3 then logLike = logLike + censorship*(alpha_s3+A+B*t) - exp(alpha_s3+A)/B*(exp(B*t)-exp(B*k2)) - exp(alpha_s2+A)/B*(exp(B*k2)-exp(B*k1)) - exp(alpha_s1+A)/B*(exp(B*k1)-1);
	else if t <= k4 then logLike = logLike + censorship*(alpha_s4+A+B*t) - exp(alpha_s4+A)/B*(exp(B*t)-exp(B*k3)) - exp(alpha_s3+A)/B*(exp(B*k3)-exp(B*k2)) - exp(alpha_s2+A)/B*(exp(B*k2)-exp(B*k1)) - exp(alpha_s1+A)/B*(exp(B*k1)-1);
        else                 logLike = logLike + censorship*(alpha_s5+A+B*t) - exp(alpha_s5+A)/B*(exp(B*t)-exp(B*k4)) - exp(alpha_s4+A)/B*(exp(B*k4)-exp(B*k3)) - exp(alpha_s3+A)/B*(exp(B*k3)-exp(B*k2)) - exp(alpha_s2+A)/B*(exp(B*k2)-exp(B*k1)) - exp(alpha_s1+A)/B*(exp(B*k1)-1);
    end;

    if 1e-3 > beta > -1e-3 then do;
	    if t <= k1      then logLike = logLike + censorship*(alpha_s1+alphai) - exp(alpha_s1+alphai)*t;
 	    else if t <= k2 then logLike = logLike + censorship*(alpha_s2+alphai) - exp(alpha_s2+alphai)*(t-k1) - exp(alpha_s1+alphai)*k1;
            else if t <= k3 then logLike = logLike + censorship*(alpha_s3+alphai) - exp(alpha_s3+alphai)*(t-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
	    else if t <= k4 then logLike = logLike + censorship*(alpha_s4+alphai) - exp(alpha_s4+alphai)*(t-k3) - exp(alpha_s3+alphai)*(k3-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
	    else                 logLike = logLike + censorship*(alpha_s5+alphai) - exp(alpha_s5+alphai)*(t-k4) - exp(alpha_s4+alphai)*(k4-k3) - exp(alpha_s3+alphai)*(k3-k2) - exp(alpha_s2+alphai)*(k2-k1) - exp(alpha_s1+alphai)*k1;
    end;

end;

	
	** average hazard ratio;
	t0 = 5;
	del = max( 0, min(-(beta*gamma_trt+alpha_trt)/(beta*gamma_time_trt), t0) );
	phi = sign(beta*gamma_time_trt)*((beta*gamma_trt+alpha_trt)**3 + (beta*gamma_time_trt*t0+(beta*gamma_trt+alpha_trt))**3 - 2*(beta*gamma_time_trt*del+(beta*gamma_trt+alpha_trt))**3)/(3*beta*gamma_time_trt); 
        c = sign(beta*gamma_time_trt)*( beta*gamma_time_trt*t0**2/2 + (beta*gamma_trt+alpha_trt)*(t0-2*del) - beta*gamma_time_trt*del**2 );	

        extra = t0*(beta*gamma_trt+alpha_trt) + beta*gamma_time_trt*t0**2/2;
        eta = 1e-3;

        estimate "trt effect" (phi+extra*eta)/(c+eta*t0);
   
    model t ~  general(loglike);
run;
quit;


