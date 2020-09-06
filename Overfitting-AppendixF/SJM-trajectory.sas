libname lib "/pine/scr/j/i/jiawei/rev11/overfitting/means_trajectory_nore";

data all&sysparm.;
infile "/pine/scr/j/i/jiawei/rev11/overfitting/data/data&sysparm..csv" delimiter = ',';
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
 model measure = t1 t2 t3 t4 trt t1*trt t2*trt t3*trt t4*trt nodegrp / solution ddfm=kr;
 *random intercept / subject=id type=un;
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
	keep variable estimate sim;
	rename variable=parameter;

   run;


ods output ParameterEstimates = lib.nore_jmParmEst&sysparm. AdditionalEstimates = lib.nore_jmEst&sysparm.;;
proc nlmixed  data = dataAll&sysparm. /*GCONV=1e-10*/ NOAD QPOINTS=10 MAXFU=10000 MAXIT=1000 REST=100;

by sim;

parms /data=init_parms&sysparm. BYDATA;

bounds sig_residual>1e-3; 

    ** longitudinal component;
    muLong = gamma_intercept + gamma_t1*t1 + gamma_t2*t2 + gamma_t3*t3 + gamma_t4*t4 + gamma_trt*trt + (gamma_t1_trt*t1+gamma_t2_trt*t2+gamma_t3_trt*t3+gamma_t4_trt*t4)*trt + gamma_nodegrp*nodegrp;
    loglike = -0.5*log(sig_residual**2) - 0.5/sig_residual**2*(measure-muLong)**2;

    L1 = 0.25;
    L2 = 0.75;
    L3 = 1.25;

   ** survival component;
   if surv_like = 1 then do;

	A1 = beta*gamma_intercept + (beta*gamma_trt+alpha_trt)*trt + alpha_nodegrp*nodegrp;
	A2 = A1 + beta*(gamma_t1 + gamma_t1_trt*trt)*L1;
	A3 = A2 + beta*(gamma_t2 + gamma_t2_trt*trt)*(L2-L1);
        A4 = A3 + beta*(gamma_t3 + gamma_t3_trt*trt)*(L3-L2);

    B1 = beta*(gamma_t1 + gamma_t1_trt*trt);
    B2 = beta*(gamma_t2 + gamma_t2_trt*trt);
    B3 = beta*(gamma_t3 + gamma_t3_trt*trt);
    B4 = beta*(gamma_t4 + gamma_t4_trt*trt);

   alphai = alpha_trt*trt + alpha_nodegrp*nodegrp;

   if beta > 1e-3 or beta < -1e-3 then do;

            if      t <= L1 then logLike = logLike + censorship*(alpha_intercept+A1+B1*t) - exp(alpha_intercept+A1)*(exp(B1*t)-1)/B1;
            else if t <= L2 then logLike = logLike + censorship*(alpha_intercept+A2+B2*(t-L1)) - exp(alpha_intercept+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_intercept+A2)*(exp(B2*(t-L1))-1)/B2;
	    else if t <= L3 then logLike = logLike + censorship*(alpha_intercept+A3+B3*(t-L2)) - exp(alpha_intercept+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_intercept+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_intercept+A3)*(exp(B3*(t-L2))-1)/B3;
	    else                 logLike = logLike + censorship*(alpha_intercept+A4+B4*(t-L3)) - exp(alpha_intercept+A1)*(exp(B1*L1)-1)/B1 - exp(alpha_intercept+A2)*(exp(B2*(L2-L1))-1)/B2 - exp(alpha_intercept+A3)*(exp(B3*(L3-L2))-1)/B3 - exp(alpha_intercept+A4)*(exp(B4*(t-L3))-1)/B4;
   end;
    
   if 1e-3 > beta > -1e-3 then do;

            logLike = logLike + censorship*(alpha_intercept+alphai) - exp(alpha_intercept+alphai)*t;
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


	phi1 = sign(beta*gamma_t1_trt)*((beta*gamma_trt + alpha_trt)**3+(beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3-(beta*gamma_t1_trt*del1+(beta*gamma_trt + alpha_trt))**3-(beta*gamma_t1_trt*del2+(beta*gamma_trt + alpha_trt))**3)/(3*beta*gamma_t1_trt); 
        phi2 = sign(beta*gamma_t2_trt)*( (beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 + (beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t2_trt*(del3-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t2_trt*(del4-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 )/(3*beta*gamma_t2_trt);
        phi3 = sign(beta*gamma_t3_trt)*( (beta*gamma_t1_trt*L1+beta*gamma_t2_trt*(L2-L1)+(beta*gamma_trt + alpha_trt))**3 + (beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t3_trt*(del5-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t3_trt*(del6-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 )/(3*beta*gamma_t3_trt);
        phi4 = sign(beta*gamma_t4_trt)*( (beta*gamma_t1_trt*L1+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t3_trt*(L3-L2)+(beta*gamma_trt + alpha_trt))**3 + (beta*gamma_t4_trt*(t0-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t4_trt*(del7-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 - (beta*gamma_t4_trt*(del8-L3)+beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))**3 )/(3*beta*gamma_t4_trt);

        c1 = sign(beta*gamma_t1_trt)*( beta*gamma_t1_trt*L1**2/2 + (beta*gamma_trt + alpha_trt)*(L1-del1-del2) - beta*gamma_t1_trt*del1**2/2 - beta*gamma_t1_trt*del2**2/2 );	
        c2 = sign(beta*gamma_t2_trt)*( beta*gamma_t2_trt*(L2-L1)**2/2 + (beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))*(L2+L1-del3-del4) - beta*gamma_t2_trt*(del3-L1)**2/2 - beta*gamma_t2_trt*(del4-L1)**2/2 );
        c3 = sign(beta*gamma_t3_trt)*( beta*gamma_t3_trt*(L3-L2)**2/2 - beta*gamma_t3_trt*(del5-L2)**2/2 - beta*gamma_t3_trt*(del6-L2)**2/2 + (beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))*(L2+L3-del5-del6) );
        c4 = sign(beta*gamma_t4_trt)*( beta*gamma_t4_trt*(t0-L3)**2/2 - beta*gamma_t4_trt*(del7-L3)**2/2 - beta*gamma_t4_trt*(del8-L3)**2/2 + (beta*gamma_t3_trt*(L3-L2)+beta*gamma_t2_trt*(L2-L1)+beta*gamma_t1_trt*L1+(beta*gamma_trt + alpha_trt))*(t0+L3-del7-del8) );

        extra = t0*(beta*gamma_trt + alpha_trt) + beta*( gamma_t1_trt*L1**2/2+gamma_t2_trt*(L2-L1)**2/2+gamma_t3_trt*(L3-L2)**2/2+gamma_t4_trt*(t0-L3)**2/2 + gamma_t3_trt*(L3-L2)*(t0-L3)+gamma_t2_trt*(L2-L1)*(t0-L2)+gamma_t1_trt*L1*(t0-L1) );
        eta = 1e-3;

        estimate "Overall Effect eta" (phi1+phi2+phi3+phi4+extra*eta)/(c1+c2+c3+c4+eta*t0);



   model t ~  general(loglike);

run;
quit;


