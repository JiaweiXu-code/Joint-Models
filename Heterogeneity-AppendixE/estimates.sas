libname lib "/pine/scr/j/i/jiawei/rev11/sigma/0.356/means";
libname pro "/pine/scr/j/i/jiawei/rev11/sigma/0.356/programs";

%macro read(n);
%do nn = 0 %to &n. %by 4;
%let m = %eval(&nn.*80 + 1);
%let mm = %eval( (&nn.+1)*80 );
%do s=&m. %to &mm.;
proc transpose data=lib.int_jmParmEst&s out=trans&s;
by sim;
id parameter;
var estimate;
run;
quit;
%end;
%end;
%mend;
%read(12);

%macro append(m);
%do s=2 %to &m.;
proc datasets library=work nolist;
append base=trans1 data=trans&s force;
run;
quit;
%end;
%mend;
%append(80);

%macro tpend(n);
%do nn = 4 %to &n. %by 4;
%let m = %eval(&nn.*80 + 1);
%let mm = %eval( (&nn.+1)*80 );
%do s=&m. %to &mm.;
proc datasets library=work nolist;
append base=trans1 data=trans&s force;
run;
quit;
%end;
%end;
%mend;
%tpend(12);

data trans;
set trans1;
order = ceil(_n_/4000);
run;

proc means data = trans;
by order;
var alpha: gamma: beta sig:;
output out=pro.zest_int;
run;

