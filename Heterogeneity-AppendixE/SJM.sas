libname lib "/pine/scr/j/i/jiawei/rev11/sigma/0.356/means_nore";

data all&sysparm.;
infile "/pine/scr/j/i/jiawei/rev11/sigma/0.356/data/data&sysparm..csv" delimiter = ',';
input ID measure time trt toltime nodegrp sim t sr r0 censorship;
log_t = log(t);
run;

proc sort data=all&sysparm.;
by sim id;
run;

data dataAll&sysparm.;
set all&sysparm.;
by sim id;
if first.id then surv_like = 1;
else surv_like = 0;
run;


ods output SolutionF = longRegParms&sysparm. CovParms = longCovParms&sysparm.;
proc mixed data = dataAll&sysparm. method=reml plots=(none);
by sim;
 model measure = time trt time*trt nodegrp / solution ddfm=kr;
* random intercept / subject=id type=un;
run;


ods output ParameterEstimates = survRegParms&sysparm.;
proc genmod data = dataAll&sysparm.;
by sim;
where surv_like = 1;
model censorship = trt nodegrp / offset = log_t dist=poisson link=log;
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
	
	*     if c and upcase(CovParm) = 'UN(1,1)' then parameter = 'intercept';
	if c                                 then parameter = lowcase(covParm);
	
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
	*where find(variable,'random','i')=0;
        if find(variable,'sig_intercept','i') then delete;
	keep variable estimate sim;
	rename variable=parameter;

   run;


ods output ParameterEstimates = lib.nore_jmParmEst&sysparm. AdditionalEstimates = lib.nore_jmEst&sysparm.;
proc nlmixed  data = dataAll&sysparm. /*GCONV=1e-10*/ NOAD QPOINTS=10 MAXFU=10000 MAXIT=1000 REST=100;

by sim;

parms /data=init_parms&sysparm. BYDATA;

bounds sig_residual>1e-3; 

  ** longitudinal component;
    muLong = gamma_intercept + gamma_time*time + gamma_trt*trt + gamma_time_trt*time*trt + gamma_nodegrp*nodegrp;
    loglike = -0.5*log(sig_residual**2) - 0.5/sig_residual**2*(measure-muLong)**2;

	** survival component;
   if surv_like = 1 then do;

	  A = beta*gamma_intercept + (beta*gamma_trt+alpha_trt)*trt + alpha_nodegrp*nodegrp;

      B = beta*(gamma_time + gamma_time_trt*trt);

	alphai = alpha_trt*trt + alpha_nodegrp*nodegrp;

    if beta > 1e-3 or beta < -1e-3 then logLike = logLike + censorship*(alpha_intercept+A+B*t) - exp(alpha_intercept+A)/B*(exp(B*t)-1);

    if 1e-3 > beta > -1e-3 then logLike = logLike + censorship*(alpha_intercept+alphai) - exp(alpha_intercept+alphai)*t;

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


