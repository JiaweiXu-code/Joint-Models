libname pro "/pine/scr/j/i/jiawei/rev11/11/programs";
libname lib "/pine/scr/j/i/jiawei/rev11/11/means";

%macro dat(no);
%do s=1 %to &no.;
data logall&s.;
infile "/pine/scr/j/i/jiawei/rev11/11/data/data&s..csv" delimiter = ',';
input ID measure time trt toltime nodegrp sim t sr r0 censorship;
run;
%end;
%mend dat;
%dat(1040);

%macro pend(n);
%do nn = 80 %to &n. %by 80;
%let m = %eval(&nn. - 80 + 1);
%do s=&nn. - 80 + 2 %to &nn.;
proc datasets library=work nolist;
append base=logall&m. data=logall&s. force;
run;
quit;
%end;

proc sort data=logall&m.;
by sim id;
run;

%let t = %eval(&nn./80);

data lsurv&t.;
set logall&m.;
by sim id;
if first.id then output;
keep sim id t censorship trt nodegrp;
run;


ods output  HomTests=lrank&t.;
proc lifetest data = lsurv&t.;
by sim;
time t*censorship(0);
strata trt / test=logrank;
run;
data logrank&t.;
set lrank&t. end=last;
retain power 0;
power+(.<probChiSq<0.05);
if last;
power=power/_n_;
run;

%end;
%mend;
%pend(1040);


%macro append(n);
%do j = 2 %to &n.;
proc datasets library=work nolist;
append base=logrank1 data=logrank&j. force;
run;
quit;
%end;
%mend;
%append(13);

data pro.logrank;
set logrank1;
run;





