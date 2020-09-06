libname lib "/pine/scr/j/i/jiawei/rev11/11/means";
libname pro "/pine/scr/j/i/jiawei/rev11/11/programs";

%macro read(no);
%do s=1 %to &no.;
data nore_jmest&s.;
set lib.jmest_nore&s.;
run;
%end;
%mend;
%read(1040);

%macro tpend(m);
%do s=2 %to &m.;
proc datasets library=work nolist;
append base=nore_jmest1 data=nore_jmest&s. force;
run;
quit;
%end;
%mend;
%tpend(1040);

data nore_jmest;
set nore_jmest1;
order = ceil(_n_/4000);
run;

   data pro.nore_powerJM;
    set nore_jmEst;
    by order;

	retain power;
	if first.order then power = 0;
	power+(cdf('normal',-estimate/StandardError)>0.95);

	if last.order;
     power=power/4000;
   run; 

